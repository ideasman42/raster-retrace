
///
/// Perform cubic curve fitting
///
/// This module takes a complete polygon and optimizes curve fitting
/// and optionally corner calculation,
/// outputting a bezier curve that fits within an error margin.
///

const USE_REFIT: bool = true;
const USE_REFIT_REMOVE: bool = true;
const CORNER_SCALE: f64 = 2.0;  // this is weak, should be made configurable.

macro_rules! unlikely { ($body:expr) => { $body } }

use ::intern::math_vector::{
    add_vnvn, sub_vnvn,
    copy_vnvn,
    dot_vnvn,
    madd_vnvn_fl,
    normalize_vn,
    normalized_vnvn_with_len,
    sq,
    zero_vn,
};

use ::intern::min_heap;

const DIMS: usize = ::intern::math_vector::DIMS;

use std::collections::LinkedList;

use super::curve_fit_single;

#[derive(Copy, Clone, PartialEq)]
pub enum TraceMode {
    Outline,
    Centerline,
}

mod types {
    use super::{
        DIMS,
    };

    pub struct Knot {
        pub next: usize,
        pub prev: usize,

        /// The index of this knot in the point array.
        ///
        /// Currently the same, access as different for now,
        /// since we may want to support different point/knot indices
        pub index: usize,

        pub no_remove: bool,
        pub is_remove: bool,
        pub is_corner: bool,

        pub handles: [f64; 2],

        /// Store the error value, to see if we can improve on it
        /// (without having to re-calculate each time)
        ///
        /// This is the error between this knot and the next.
        pub fit_error_sq_next: f64,

        /// Initially point to contiguous memory, however we may re-assign.
        pub tan: [usize; 2],
    }

    pub struct PointData<'a> {
        /// note, can't use points.len(),
        /// since this may be doubled for cyclic curves
        pub points: &'a Vec<[f64; DIMS]>,
        pub points_len: usize,

        /// This array may be doubled as well.
        pub points_length_cache: &'a Vec<f64>,

        pub tangents: &'a Vec<[f64; DIMS]>,
    }
}

use self::types::{
    Knot,
    PointData,
};

const INVALID: usize = ::std::usize::MAX;

/// Find the knot furthest from the line between \a knot_l & \a knot_r.
/// This is to be used as a split point.
fn knot_find_split_point_on_axis(
    pd: &PointData,
    knots: &Vec<Knot>,
    k_prev: &Knot,
    k_next: &Knot,
    plane_no: &[f64; DIMS],
) -> usize {
    let mut split_point: usize = INVALID;
    let mut split_point_dist_best: f64 = -::std::f64::MAX;

    let knots_end = knots.len() - 1;
    let mut k_step = k_prev.index;
    loop {
        if k_step != knots_end {
            k_step += 1;
        } else {
            // wrap around
            k_step = 0;
        }

        if k_step != k_next.index {
            let knot = &knots[k_step];
            let split_point_dist_test = dot_vnvn(plane_no, &pd.points[knot.index]);
            if split_point_dist_test > split_point_dist_best {
                split_point_dist_best = split_point_dist_test;
                split_point = knot.index;
            }
        } else {
            break;
        }
    }

    return split_point;
}


fn knot_remove_error_value(
    tan_l: &[f64; DIMS],
    tan_r: &[f64; DIMS],
    points_offset: &[[f64; DIMS]],
    points_offset_length_cache: &[f64],
) -> (f64, usize, [f64; 2]) {
    let ((error_sq, error_index), handle_factor_l, handle_factor_r) =
        curve_fit_single::curve_fit_cubic_to_points_single(
            points_offset, points_offset_length_cache,
            tan_l, tan_r,
            );
    return (
        error_sq, error_index,
        [dot_vnvn(tan_l, &sub_vnvn(&handle_factor_l, &points_offset[0])),
         dot_vnvn(tan_r, &sub_vnvn(&handle_factor_r, &points_offset[points_offset.len() - 1]))],
    );
}

fn knot_calc_curve_error_value_and_index(
    pd: &PointData,
    knot_l: &Knot, knot_r: &Knot,
    tan_l: &[f64; DIMS],
    tan_r: &[f64; DIMS],
) -> (f64, usize, [f64; 2]) {
    let points_offset_len =
        if knot_l.index < knot_r.index {
            knot_r.index - knot_l.index
        } else {
            (knot_r.index + pd.points_len) - knot_l.index
        } + 1;

    if points_offset_len != 2 {
        let points_offset_end = knot_l.index + points_offset_len;
        let mut result = knot_remove_error_value(
            tan_l, tan_r,
            &pd.points[knot_l.index..points_offset_end],
            &pd.points_length_cache[knot_l.index..points_offset_end],
            );

        // Adjust the offset index to the global index & wrap if needed.
        result.1 += knot_l.index;
        if result.1 >= pd.points_len {
            result.1 -= pd.points_len;
        }
        return result;
    } else {
        // No points between, use 1/3 handle length with no error as a fallback.
        debug_assert!(points_offset_len == 2);
        let handle_len = pd.points_length_cache[knot_l.index] / 3.0;
        return (0.0, knot_l.index, [handle_len, handle_len]);
    }
}

fn knot_calc_curve_error_value(
    pd: &PointData,
    knot_l: &Knot, knot_r: &Knot,
    tan_l: &[f64; DIMS],
    tan_r: &[f64; DIMS],
) -> (f64, [f64; 2]) {
    let points_offset_len =
        if knot_l.index < knot_r.index {
            knot_r.index - knot_l.index
        } else {
            (knot_r.index + pd.points_len) - knot_l.index
        } + 1;

    if points_offset_len != 2 {
        let points_offset_end = knot_l.index + points_offset_len;
        let result = knot_remove_error_value(
            tan_l, tan_r,
            &pd.points[knot_l.index..points_offset_end],
            &pd.points_length_cache[knot_l.index..points_offset_end],
            );
        return (result.0, result.2);
    } else {
        // No points between, use 1/3 handle length with no error as a fallback.
        debug_assert!(points_offset_len == 2);
        let handle_len = pd.points_length_cache[knot_l.index] / 3.0;
        return (0.0, [handle_len, handle_len]);
    }
}

mod refine_remove {
    use super::{
        INVALID,
        knot_calc_curve_error_value,
    };
    use super::types::{
        Knot,
        PointData,
    };
    use ::intern::min_heap;

    // Store adjacent handles in the case this is removed
    // could make this part of the knot array but its logically
    // more clear whats going on if its kept separate.
    #[derive(Copy, Clone)]
    struct KnotRemoveState {
        // handles for prev/next knots
        index: usize,
        handles: [f64; 2],
    }

    fn knot_remove_error_recalculate(
        pd: &PointData,
        heap: &mut min_heap::MinHeap<f64, KnotRemoveState>,
        knots: &Vec<Knot>,
        knots_handle: &mut Vec<min_heap::NodeHandle>,
        k_curr: &Knot,
        error_max_sq: f64,
    ) {
        debug_assert!(k_curr.no_remove == false);

        let (fit_error_max_sq, handles) = {
            let k_prev = &knots[k_curr.prev];
            let k_next = &knots[k_curr.next];

            knot_calc_curve_error_value(
                pd, k_prev, k_next,
                &pd.tangents[k_prev.tan[1]],
                &pd.tangents[k_next.tan[0]])
        };

        let k_curr_heap_node = &mut knots_handle[k_curr.index];
        if fit_error_max_sq < error_max_sq {
            if *k_curr_heap_node != min_heap::NodeHandle::INVALID {
                heap.remove(*k_curr_heap_node);
            }
            *k_curr_heap_node = heap.insert(
                fit_error_max_sq,
                KnotRemoveState {
                    index: k_curr.index,
                    handles: handles,
                },
            );
        } else {
            if *k_curr_heap_node != min_heap::NodeHandle::INVALID {
                heap.remove(*k_curr_heap_node);
                *k_curr_heap_node = min_heap::NodeHandle::INVALID;
            }
        }
    }

    pub fn curve_incremental_simplify(
        pd: &PointData,
        knots: &mut Vec<Knot>,
        knots_handle: &mut Vec<min_heap::NodeHandle>,
        knots_len_remaining: &mut usize,
        error_max_sq: f64,
    ) {
        let mut heap = min_heap::MinHeap::<f64, KnotRemoveState>::with_capacity(knots.len());

        for k_index in 0..knots.len() {
            let k_curr = &knots[k_index];
            if (k_curr.no_remove == false) &&
               (k_curr.is_remove == false) &&
               (k_curr.is_corner == false)
            {
                knot_remove_error_recalculate(
                    pd, &mut heap, knots, knots_handle, k_curr, error_max_sq);
            }
        }

        while let Some((error_sq, r)) = heap.pop_min_with_value() {
            knots_handle[r.index] = min_heap::NodeHandle::INVALID;

            let k_next_index;
            let k_prev_index;
            {
                // let r: &mut remove_states[r_index];
                let mut k_curr: &mut Knot = &mut knots[r.index];

                if unlikely!(*knots_len_remaining <= 2) {
                    continue;
                }

                k_next_index = k_curr.next;
                k_prev_index = k_curr.prev;

                k_curr.is_remove = true;

                if cfg!(debug_assertions) {
                    k_curr.next = INVALID;
                    k_curr.prev = INVALID;
                }
            }
            knots[k_prev_index].handles[1] = r.handles[0];
            knots[k_next_index].handles[0] = r.handles[1];

            debug_assert!(error_sq <= error_max_sq);

            knots[k_prev_index].fit_error_sq_next = error_sq;
            // Remove ourselves
            knots[k_next_index].prev = k_prev_index;
            knots[k_prev_index].next = k_next_index;


            for k_iter_index in &[k_prev_index, k_next_index] {
                let k_iter = &knots[*k_iter_index];
                if (k_iter.no_remove == false) &&
                   (k_iter.is_corner == false) &&
                   (k_iter.prev != INVALID) &&
                   (k_iter.next != INVALID)
                {
                    knot_remove_error_recalculate(
                        pd, &mut heap, knots, knots_handle, k_iter, error_max_sq);
                }
            }

            *knots_len_remaining -= 1;
        }
        drop(heap);
    }
}
// end refine_remove


mod refine_refit {

    use super::{
        INVALID,
        USE_REFIT_REMOVE,
        knot_calc_curve_error_value,
        knot_calc_curve_error_value_and_index,
    };
    use super::types::{
        Knot,
        PointData,
    };
    use ::intern::min_heap;

    #[derive(Copy, Clone)]
    struct KnotRefitState {
        index: usize,
        // When INVALID - remove this item
        index_refit: usize,

        // Handles for prev/next knots
        handle_pair: [[f64; 2]; 2],

        fit_error_max_sq_pair: [f64; 2],
    }

    fn knot_refit_error_recalculate(
        pd: &PointData,
        heap: &mut min_heap::MinHeap<f64, KnotRefitState>,
        knots: &Vec<Knot>,
        knots_handle: &mut Vec<min_heap::NodeHandle>,
        k_curr: &Knot,
        error_max_sq: f64,
        use_optimize_exhaustive: bool,
    ) {
        debug_assert!(k_curr.no_remove == false);

        let k_curr_heap_node = &mut knots_handle[k_curr.index];

        let k_prev = &knots[k_curr.prev];
        let k_next = &knots[k_curr.next];

        let mut k_refit_index;

        // Support re-fitting to remove points
        {

            let (fit_error_max_sq, fit_error_index, handles) =
                knot_calc_curve_error_value_and_index(
                    pd, k_prev, k_next,
                    &pd.tangents[k_prev.tan[1]],
                    &pd.tangents[k_next.tan[0]],
                    );

            if USE_REFIT_REMOVE && fit_error_max_sq < error_max_sq {
                if *k_curr_heap_node != min_heap::NodeHandle::INVALID {
                    heap.remove(*k_curr_heap_node);
                }

                // Always perform removal before refitting, (make a negative number)
                *k_curr_heap_node = heap.insert(
                    // Weight for the greatest improvement
                    fit_error_max_sq - error_max_sq,
                    KnotRefitState {
                        index: k_curr.index,
                        // INVALID == remove
                        index_refit: INVALID,
                        // 0.0 == unused
                        handle_pair: [[handles[0], 0.0], [0.0, handles[1]]],
                        fit_error_max_sq_pair: [fit_error_max_sq, fit_error_max_sq],
                    }
                );
                return;
            }

            // Use the largest point of difference when removing
            // as the target to refit to.
            k_refit_index = fit_error_index;
        }

        if !use_optimize_exhaustive {
            if (k_refit_index == INVALID) || (k_refit_index == k_curr.index) {
                if *k_curr_heap_node != min_heap::NodeHandle::INVALID {
                    heap.remove(*k_curr_heap_node);
                    *k_curr_heap_node = min_heap::NodeHandle::INVALID;
                    return;
                }
            }
        }

        let cost_sq_src_max = k_prev.fit_error_sq_next.max(k_curr.fit_error_sq_next);
        debug_assert!(cost_sq_src_max <= error_max_sq);

        // Specialized function to avoid duplicate code
        fn knot_calc_curve_error_value_pair_above_error_or_none(
            pd: &PointData, k_prev: &Knot, k_refit: &Knot, k_next: &Knot, error_max_sq: f64,
        ) -> Option<([f64; 2], f64, [f64; 2], f64)> {
            let (fit_error_prev, handles_prev) =
                knot_calc_curve_error_value(
                    pd, k_prev, k_refit,
                    &pd.tangents[k_prev.tan[1]],
                    &pd.tangents[k_refit.tan[0]],
                );

            if fit_error_prev < error_max_sq {
                let (fit_error_next, handles_next) =
                    knot_calc_curve_error_value(
                        pd, k_refit, k_next,
                        &pd.tangents[k_refit.tan[1]],
                        &pd.tangents[k_next.tan[0]],
                    );
                if fit_error_next < error_max_sq {
                    return Some((
                        handles_prev, fit_error_prev,
                        handles_next, fit_error_next,
                    ));
                }
            }
            return None;
        }

        // Instead of using the highest error value,
        // search for *every* possible split point and test it.
        // This is _not_ meant for typical usage (since its obviously very in-efficient).
        //
        // Nevertheless its interesting to have a way to attempt the best possible result.

        // cache result of 'knot_calc_curve_error_value_pair_above_error_or_none'
        let mut refit_result_or_none: Option<([f64; 2], f64, [f64; 2], f64)> = None;

        if use_optimize_exhaustive {

            // loop over inner knots
            let mut k_test_index = k_prev.index + 1;

            // start with current state
            let mut cost_sq_best = cost_sq_src_max;

            loop {
                if k_test_index == knots.len() {
                    k_test_index = 0;
                }
                if k_test_index == k_next.index {
                    break;
                }

                if k_test_index != k_curr.index {
                    if let Some(fit_result_test) =
                        knot_calc_curve_error_value_pair_above_error_or_none(
                            pd, k_prev, &knots[k_test_index], k_next, cost_sq_best)
                    {
                        let cost_sq_test_prev = fit_result_test.1;
                        let cost_sq_test_next = fit_result_test.3;
                        cost_sq_best = cost_sq_test_prev.max(cost_sq_test_next);
                        k_refit_index = k_test_index;

                        // Result for re-use if this is the best fit.
                        refit_result_or_none = Some(fit_result_test);
                    }
                }
                k_test_index += 1;
            }
        } else {
            refit_result_or_none =
                knot_calc_curve_error_value_pair_above_error_or_none(
                    pd, k_prev, &knots[k_refit_index], k_next, cost_sq_src_max)
        }
        // end exhaustive test

        if let Some((
            handles_prev, fit_error_dst_prev,
            handles_next, fit_error_dst_next,
        )) = refit_result_or_none {

            if *k_curr_heap_node != min_heap::NodeHandle::INVALID {
                heap.remove(*k_curr_heap_node);
            }

            let fit_error_dst_max_sq =
                fit_error_dst_prev.max(fit_error_dst_next);
            debug_assert!(fit_error_dst_max_sq < cost_sq_src_max);

            *k_curr_heap_node = heap.insert(
                // Weight for the greatest improvement.
                cost_sq_src_max - fit_error_dst_max_sq,
                KnotRefitState {
                    index: k_curr.index,
                    index_refit: k_refit_index,
                    handle_pair: [handles_prev, handles_next],
                    fit_error_max_sq_pair: [fit_error_dst_prev, fit_error_dst_next],
                }
            );
            return;
        }

        if *k_curr_heap_node != min_heap::NodeHandle::INVALID {
            heap.remove(*k_curr_heap_node);
            *k_curr_heap_node = min_heap::NodeHandle::INVALID;
        }
    }

    pub fn curve_incremental_simplify_refit(
        pd: &PointData,
        knots: &mut Vec<Knot>,
        knots_handle: &mut Vec<min_heap::NodeHandle>,
        knots_len_remaining: &mut usize,
        error_max_sq: f64,
        use_optimize_exhaustive: bool,
    ) {
        let mut heap =
            min_heap::MinHeap::<f64, KnotRefitState>::with_capacity(*knots_len_remaining);

        for k_index in 0..knots.len() {
            let k_curr = &knots[k_index];
            if (k_curr.no_remove == false) &&
               (k_curr.is_remove == false) &&
               (k_curr.is_corner == false)
            {
                knot_refit_error_recalculate(
                    pd, &mut heap, knots, knots_handle, k_curr,
                    error_max_sq, use_optimize_exhaustive);
            }
        }


        while let Some(r) = heap.pop_min() {
            knots_handle[r.index] = min_heap::NodeHandle::INVALID;

            let k_prev_index;
            let k_next_index;
            {
                {
                    let k_old = &knots[r.index];
                    k_prev_index = k_old.prev;
                    k_next_index = k_old.next;
                }

                if r.index_refit == INVALID {
                    // remove
                } else {
                    let k_refit = &mut knots[r.index_refit];
                    k_refit.handles[0] = r.handle_pair[0][1];
                    k_refit.handles[1] = r.handle_pair[1][0];
                }

                knots[k_prev_index].handles[1] = r.handle_pair[0][0];
                knots[k_next_index].handles[0] = r.handle_pair[1][1];
            }
            // finished with 'r'

            // XXX, check this is OK
            if unlikely!(*knots_len_remaining <= 2) {
                continue;
            }

            {
                let k_old = &mut knots[r.index];
                k_old.next = INVALID;
                k_old.prev = INVALID;
                k_old.is_remove = true;
            }

            if r.index_refit == INVALID {
                knots[k_next_index].prev = k_prev_index;
                knots[k_prev_index].next = k_next_index;

                knots[k_prev_index].fit_error_sq_next = r.fit_error_max_sq_pair[0];

                *knots_len_remaining -= 1;
            } else {
                // Remove ourselves
                knots[k_next_index].prev = r.index_refit;
                knots[k_prev_index].next = r.index_refit;

                knots[k_prev_index].fit_error_sq_next = r.fit_error_max_sq_pair[0];

                let k_refit = &mut knots[r.index_refit];
                k_refit.prev = k_prev_index;
                k_refit.next = k_next_index;

                k_refit.fit_error_sq_next = r.fit_error_max_sq_pair[1];

                k_refit.is_remove = false;
            }

            for k_iter_index in &[k_prev_index, k_next_index] {
                let k_iter = &knots[*k_iter_index];
                if (k_iter.no_remove == false) &&
                   (k_iter.is_corner == false) &&
                   (k_iter.prev != INVALID) &&
                   (k_iter.next != INVALID)
                {
                    knot_refit_error_recalculate(
                        pd, &mut heap, knots, knots_handle, k_iter,
                        error_max_sq, use_optimize_exhaustive);
                }
            }
        }

        drop(heap);
    }
}
// end refine_refit

mod refine_corner {
    use super::{
        INVALID,
        knot_calc_curve_error_value,
        knot_find_split_point_on_axis,
    };
    use super::types::{
        Knot,
        PointData,
    };
    use ::intern::math_vector::{
        dot_vnvn,
        len_squared_vnvn,
        project_vnvn_normalized,
        sub_vnvn,
    };
    use ::intern::min_heap;

    // Result of collapsing a corner.
    #[derive(Copy, Clone)]
    struct KnotCornerState {
        index: usize,
        // Merge adjacent handles into this one (may be shared with the 'index').
        index_pair: [usize; 2],

        // Handles for prev/next knots.
        handle_pair: [[f64; 2]; 2],

        fit_error_max_sq_pair: [f64; 2],
    }

    /// (Re)calculate the error incurred from turning this into a corner.
    fn knot_corner_error_recalculate(
        pd: &PointData,
        heap: &mut min_heap::MinHeap<f64, KnotCornerState>,
        knots_handle: &mut Vec<min_heap::NodeHandle>,
        k_split: &Knot,
        k_prev: &Knot,
        k_next: &Knot,
        error_max_sq: f64,
    ) {
        debug_assert!(
            (k_prev.no_remove == false) &&
            (k_next.no_remove == false)
        );

        let k_split_heap_node = &mut knots_handle[k_split.index];

        // Test skipping 'k_prev' by using points (k_prev.prev to k_split).
        {
            let (fit_error_dst_prev, handles_prev) =
                knot_calc_curve_error_value(
                    pd, k_prev, k_split,
                    &pd.tangents[k_prev.tan[1]],
                    &pd.tangents[k_prev.tan[1]],
                    );
            if fit_error_dst_prev < error_max_sq {
                let (fit_error_dst_next, handles_next) =
                    knot_calc_curve_error_value(
                        pd, k_split, k_next,
                        &pd.tangents[k_next.tan[0]],
                        &pd.tangents[k_next.tan[0]],
                        );
                if fit_error_dst_next < error_max_sq {

                    if *k_split_heap_node != min_heap::NodeHandle::INVALID {
                        heap.remove(*k_split_heap_node);
                    }

                    // _must_ be assigned to k_split, later
                    *k_split_heap_node = heap.insert(
                        // Weight for the greatest improvement.
                        fit_error_dst_prev.max(fit_error_dst_next),
                        KnotCornerState {
                            index: k_split.index,
                            // Need to store handle lengths for both sides
                            index_pair: [k_prev.index, k_next.index],
                            handle_pair: [handles_prev, handles_next],
                            fit_error_max_sq_pair: [fit_error_dst_prev, fit_error_dst_next],
                        }
                    );

                    return;
                }
            }
        }

        if *k_split_heap_node != min_heap::NodeHandle::INVALID {
            heap.remove(*k_split_heap_node);
            *k_split_heap_node = min_heap::NodeHandle::INVALID;
        }
    }

    // Attempt to collapse close knots into corners,
    // as long as they fall below the error threshold.
    pub fn curve_incremental_simplify_corners(
        pd: &PointData,
        knots: &mut Vec<Knot>,
        knots_handle: &mut Vec<min_heap::NodeHandle>,
        knots_len_remaining: &mut usize,
        error_max_sq: f64,
        error_sq_collapse_max: f64,
        corner_angle: f64,
    ) {
        // don't pre-allocate, since its likely there are no corners
        let mut heap = min_heap::MinHeap::<f64, KnotCornerState>::with_capacity(0);

        let corner_angle_cos = corner_angle.cos();

        for k_prev_index in 0..knots.len() {
            if let Some((k_prev, k_next)) = {
                let k_prev: &Knot = &knots[k_prev_index];

                if (k_prev.is_remove == false) &&
                   (k_prev.no_remove == false) &&
                   (k_prev.next != INVALID) &&
                   (knots[k_prev.next].no_remove == false)
                {
                    Some((k_prev, &knots[k_prev.next]))
                } else {
                    None
                }
            }
            {
                // Angle outside threshold
                if dot_vnvn(
                    &pd.tangents[k_prev.tan[0]],
                    &pd.tangents[k_next.tan[1]]) < corner_angle_cos
                {
                    // Measure distance projected onto a plane,
                    //since the points may be offset along their own tangents.
                    let plane_no = sub_vnvn(
                        &pd.tangents[k_next.tan[0]],
                        &pd.tangents[k_prev.tan[1]],
                        );

                    // Compare 2x so as to allow both to be changed
                    // by maximum of `error_sq_collapse_max`.
                    let k_split_index = knot_find_split_point_on_axis(
                        pd,
                        knots,
                        k_prev,
                        k_next,
                        &plane_no,
                        );

                    if k_split_index != INVALID {
                        let co_prev  = &pd.points[k_prev.index];
                        let co_next  = &pd.points[k_next.index];
                        let co_split = &pd.points[k_split_index];

                        let k_proj_ref = project_vnvn_normalized(
                            co_prev, &pd.tangents[k_prev.tan[1]]);
                        let k_proj_split = project_vnvn_normalized(
                            co_split, &pd.tangents[k_prev.tan[1]]);

                        if len_squared_vnvn(
                            &k_proj_ref, &k_proj_split) < error_sq_collapse_max
                        {
                            let k_proj_ref = project_vnvn_normalized(
                                co_next, &pd.tangents[k_next.tan[0]]);
                            let k_proj_split = project_vnvn_normalized(
                                co_split, &pd.tangents[k_next.tan[0]]);

                            if len_squared_vnvn(
                                &k_proj_ref, &k_proj_split) < error_sq_collapse_max
                            {
                                knot_corner_error_recalculate(
                                    pd,
                                    &mut heap,
                                    knots_handle,
                                    &knots[k_split_index],
                                    k_prev,
                                    k_next,
                                    error_max_sq,
                                    );
                            }
                        }
                    }
                }
            }
        }

        while let Some(c) = heap.pop_min() {
            knots_handle[c.index] = min_heap::NodeHandle::INVALID;

            let k_split_index = c.index;
            let k_prev_index = c.index_pair[0];
            let k_next_index = c.index_pair[1];

            let tan_prev;
            let tan_next;


            {
                let k_prev  = &mut knots[k_prev_index];
                k_prev.next = k_split_index;
                k_prev.handles[1]  = c.handle_pair[0][0];
                tan_prev = k_prev.tan[1];

                debug_assert!(c.fit_error_max_sq_pair[0] <= error_max_sq);
                k_prev.fit_error_sq_next = c.fit_error_max_sq_pair[0];
            }

            {
                let k_next  = &mut knots[k_next_index];
                k_next.prev = k_split_index;
                tan_next = k_next.tan[0];

                k_next.handles[0] = c.handle_pair[1][1];
            }

            // Remove while collapsing
            {
                let k_split = &mut knots[k_split_index];

                // Insert
                k_split.is_remove = false;
                k_split.is_corner = true;

                k_split.prev = k_prev_index;
                k_split.next = k_next_index;

                // Update tangents
                k_split.tan[0] = tan_prev; // knots[k_prev_index].tan[1];
                k_split.tan[1] = tan_next; // knots[k_next_index].tan[0];

                // Own handles
                k_split.handles[0] = c.handle_pair[0][1];
                k_split.handles[1] = c.handle_pair[1][0];

                debug_assert!(c.fit_error_max_sq_pair[1] <= error_max_sq);
                k_split.fit_error_sq_next = c.fit_error_max_sq_pair[1];
            }

            *knots_len_remaining += 1;
        }


        drop(heap);
    }
}

// end refine_corner


pub fn fit_poly_single(
    // points_orig: &[[f64; 2]],
    points_orig: &Vec<[f64; DIMS]>,
    is_cyclic: bool,
    error_threshold: f64,
    corner_angle: f64,
    use_optimize_exhaustive: bool,
) -> Vec<[[f64; DIMS]; 3]> {
    use ::intern::math_vector::{
        is_finite_vn,
    };

    // Double size to allow extracting wrapped contiguous slices across start/end boundaries.
    let knots_len = points_orig.len();
    let points_len = points_orig.len();
    let points = if is_cyclic {
        [points_orig.as_slice(), points_orig.as_slice()].concat()
    } else {
        // TODO, we don't need to duplicate here,
        // find a way to use the original array!
        [points_orig.as_slice()].concat()
    };

    // del_var!(points_orig);  // TODO

    let mut knots: Vec<Knot> =
        Vec::with_capacity(knots_len);
    let mut knots_handle: Vec<min_heap::NodeHandle> =
        vec![min_heap::NodeHandle::INVALID; knots_len];

    let use_corner = corner_angle < ::std::f64::consts::PI;

    for i in 0..knots_len {
        assert!(is_finite_vn(&points_orig[i]));
        knots.push(Knot {
            next: i.wrapping_add(1),
            prev: i.wrapping_sub(1),
            index: i,
            no_remove: false,
            is_remove: false,
            is_corner: false,
            handles: [-1.0, -1.0], // dummy
            fit_error_sq_next: 0.0,
            tan: [i * 2, i * 2 + 1],
        });
    }

    if is_cyclic {
        let i_last = knots.len() - 1;
        knots[0].prev = i_last;
        knots[i_last].next = 0;
    } else {
        let i_last = knots.len() - 1;
        knots[0].prev = INVALID;
        knots[i_last].next = INVALID;

        knots[0].no_remove = true;
        knots[i_last].no_remove = true;
    }

    // All values will be written to, simplest to initialize to dummy values for now.
    let mut points_length_cache: Vec<f64> = vec![-1.0; points_len * if is_cyclic { 2 } else { 1 }];
    let mut tangents: Vec<[f64; DIMS]> = vec![[-1.0; DIMS]; knots_len * 2];

    // Initialize tangents,
    // also set the values for knot handles since some may not collapse.

    if knots_len < 2 {
        for (i, k) in (&mut knots).iter_mut().enumerate() {
            zero_vn(&mut tangents[k.tan[0]]);
            zero_vn(&mut tangents[k.tan[1]]);
            k.handles[0] = 0.0;
            k.handles[1] = 0.0;
            points_length_cache[i] = 0.0;
        }
    } else if is_cyclic {
        let (mut tan_prev, mut len_prev) = normalized_vnvn_with_len(
            &points[knots_len - 2], &points[knots_len - 1]);

        let mut i_curr = knots.len() - 1;
        for i_next in 0..knots.len() {
            let k = &mut knots[i_curr];

            let (tan_next, len_next) = normalized_vnvn_with_len(
                &points[i_curr], &points[i_next]);
            points_length_cache[i_next] = len_next;

            let mut t = add_vnvn(&tan_prev, &tan_next);
            normalize_vn(&mut t);
            assert!(is_finite_vn(&t));
            copy_vnvn(&mut tangents[k.tan[0]], &t);
            copy_vnvn(&mut tangents[k.tan[1]], &t);

            k.handles[0] = len_prev /  3.0;
            k.handles[1] = len_next / -3.0;

            copy_vnvn(&mut tan_prev, &tan_next);
            len_prev = len_next;
            i_curr = i_next;
        }
    } else {
        points_length_cache[0] = 0.0;
        let (mut tan_prev, mut len_prev) = normalized_vnvn_with_len(
            &points[0], &points[1]);
        points_length_cache[1] = len_prev;

        copy_vnvn(&mut tangents[knots[0].tan[0]], &tan_prev);
        copy_vnvn(&mut tangents[knots[0].tan[1]], &tan_prev);
        knots[0].handles[0] = len_prev /  3.0;
        knots[0].handles[1] = len_prev / -3.0;

        let mut i_curr = 1;
        for i_next in 2..knots.len() {
            let k = &mut knots[i_curr];
            let (tan_next, len_next) = normalized_vnvn_with_len(
                &points[i_curr], &points[i_next]);
            points_length_cache[i_next] = len_next;

            let mut t = add_vnvn(&tan_prev, &tan_next);
            normalize_vn(&mut t);
            assert!(is_finite_vn(&t));
            copy_vnvn(&mut tangents[k.tan[0]], &t);
            copy_vnvn(&mut tangents[k.tan[1]], &t);

            k.handles[0] = len_prev /  3.0;
            k.handles[1] = len_next / -3.0;

            copy_vnvn(&mut tan_prev, &tan_next);
            len_prev = len_next;
            i_curr = i_next;
        }
        // use prev as next since they're copied above
        copy_vnvn(&mut tangents[knots[knots_len - 1].tan[0]], &tan_prev);
        copy_vnvn(&mut tangents[knots[knots_len - 1].tan[1]], &tan_prev);

        knots[knots_len - 1].handles[0] = len_prev /  3.0;
        knots[knots_len - 1].handles[1] = len_prev / -3.0;
    }

    if is_cyclic {
        // TODO, perhaps this can be done more elegantly?
        for i in 0..points_len {
            points_length_cache[i + points_len] = points_length_cache[i];
        }
    }

    let mut knots_len_remaining = knots.len();
    let pd = PointData {
        points: &points,
        points_len: points_len,
        points_length_cache: &points_length_cache,
        tangents: &tangents,
    };

    // `curve_incremental_simplify_refit` can be called here, but its very slow
    // just remove all within the threshold first.
    refine_remove::curve_incremental_simplify(
        &pd, &mut knots, &mut knots_handle, &mut knots_len_remaining,
        sq(error_threshold));

    if use_corner {
        refine_corner::curve_incremental_simplify_corners(
            &pd, &mut knots, &mut knots_handle, &mut knots_len_remaining,
            sq(error_threshold), sq(error_threshold * CORNER_SCALE),
            corner_angle,
            );
    }

    debug_assert!(knots_len_remaining >= 2);

    if USE_REFIT {
        refine_refit::curve_incremental_simplify_refit(
            &pd, &mut knots, &mut knots_handle, &mut knots_len_remaining,
            sq(error_threshold), use_optimize_exhaustive);
    }

    debug_assert!(knots_len_remaining >= 2);

    let mut cubic_array: Vec<[[f64; DIMS]; 3]> = Vec::with_capacity(knots_len_remaining);

    {
        let k_first_index: usize = {
            let mut i_search = INVALID;
            for (i, k) in knots.iter().enumerate() {
                if k.is_remove == false {
                    i_search = i;
                    break;
                }
            }
            debug_assert!(i_search != INVALID);
            i_search
        };

        let mut k_index = k_first_index;
        for _ in 0..knots_len_remaining {
            let k = &knots[k_index];
            let p = &points[k.index];

            // assert!(k.handles[0].is_finite());
            // assert!(k.handles[1].is_finite());

            cubic_array.push([
                madd_vnvn_fl(p, &tangents[k.tan[0]], k.handles[0]),
                *p,
                madd_vnvn_fl(p, &tangents[k.tan[1]], k.handles[1]),
            ]);

            k_index = k.next;
        }
    }

    return cubic_array;
}


pub fn fit_poly_list(
    poly_list_src: LinkedList<(bool, Vec<[f64; DIMS]>)>,
    error_threshold: f64,
    corner_angle: f64,
    use_optimize_exhaustive: bool,
) -> LinkedList<(bool, Vec<[[f64; DIMS]; 3]>)> {
    let mut curve_list_dst: LinkedList<(bool, Vec<[[f64; DIMS]; 3]>)> = LinkedList::new();

    // Single threaded (we may want to allow users to force this).
    if poly_list_src.len() <= 1 {
        for (is_cyclic, poly_src) in poly_list_src {
            let poly_dst = fit_poly_single(
                &poly_src, is_cyclic, error_threshold,
                corner_angle, use_optimize_exhaustive);
            println!("{} -> {}", poly_src.len(), poly_dst.len());
            curve_list_dst.push_back((is_cyclic, poly_dst));
        }
    } else {
        use std::thread;

        let mut join_handles = Vec::with_capacity(poly_list_src.len());
        let mut poly_vec_src = Vec::with_capacity(poly_list_src.len());

        for poly_src in poly_list_src {
            poly_vec_src.push(poly_src);
        }

        // sort length for more even threading
        // and so larger at the end so they are popped off and handled first,
        // smaller ones can be handled when other processors are free.
        poly_vec_src.sort_by(|a, b| a.1.len().cmp(&b.1.len()));

        while let Some((is_cyclic, poly_src_clone)) = poly_vec_src.pop() {
            join_handles.push(thread::spawn(move || {
                let poly_dst = fit_poly_single(
                    &poly_src_clone, is_cyclic, error_threshold,
                    corner_angle, use_optimize_exhaustive);
                println!("{} -> {}", poly_src_clone.len(), poly_dst.len());
                (is_cyclic, poly_dst)
            }));
        }

        for child in join_handles {
            curve_list_dst.push_back(child.join().unwrap());
        }
    }

    return curve_list_dst;
}
