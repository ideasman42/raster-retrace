///
/// Curve Re-fitting Method
/// =======================
///
/// This is a more processor-intensive method of fitting compared to
/// direct curve fitting, and works as follows:
///
/// - First iteratively remove all points under the error threshold.
///
/// - If corner calculation is enabled:
///   - Find adjacent knots that exceed the angle limit.
///   - Collapse the knots into one, or remove entirely
///     (depending on their error values).
///
/// - Run a re-fit pass, where knots are re-positioned between their adjacent knots
///   when their re-fit position has a lower 'error'.
///   While re-fitting, remove knots that fall below the error threshold.
///
/// This module contains three refinement algorithms:
/// - `refine_remove`: Incremental point removal based on error threshold
/// - `refine_refit`: Re-positioning knots to minimize error
/// - `refine_corner`: Corner detection and collapse
///

use ::intern::math_vector::{
    dot_vnvn,
    len_squared_vn,
    normalize_vn,
    project_plane_vnvn_normalized,
    sub_vnvn,
};

use super::curve_fit_single;

/// Number of dimensions for curve fitting (2D points).
const DIMS: usize = ::intern::math_vector::DIMS;

/// Sentinel value representing an invalid index (used for prev/next in non-cyclic curves).
pub const INVALID: usize = ::std::usize::MAX;

/// Type definitions for curve refinement.
pub mod types {
    use super::DIMS;

    /// A knot in the curve fitting linked list.
    ///
    /// Knots form a doubly-linked list where each knot represents a control point
    /// in the simplified curve. During refinement, knots can be removed, repositioned,
    /// or marked as corners.
    pub struct Knot {
        /// Index of next knot in the linked list (INVALID if none).
        pub next: usize,
        /// Index of previous knot in the linked list (INVALID if none).
        pub prev: usize,

        /// The index of this knot in the point array.
        ///
        /// Currently the same, access as different for now,
        /// since we may want to support different point/knot indices
        pub index: usize,

        /// If true, this knot cannot be removed (e.g., endpoints of non-cyclic curves).
        pub no_remove: bool,
        /// If true, this knot has been removed from the active list.
        pub is_remove: bool,
        /// If true, this knot is a corner (tangents are not continuous).
        pub is_corner: bool,

        /// Handle lengths [left, right] for the bezier curves meeting at this knot.
        ///
        /// These are signed values: positive extends in the tangent direction,
        /// negative extends opposite. Used to compute control points as:
        /// `control_point = knot_position + tangent * handle_length`
        pub handles: [f64; 2],

        /// Store the error value, to see if we can improve on it
        /// (without having to re-calculate each time)
        ///
        /// This is the error between this knot and the next.
        pub fit_error_sq_next: f64,

        /// Indices into the tangent array [incoming, outgoing].
        ///
        /// Initially these point to contiguous memory (knot_index * 2),
        /// but may be reassigned when corners are created.
        pub tan: [usize; 2],
    }

    /// Immutable point data used during curve refinement.
    ///
    /// Contains references to the original points, cached segment lengths,
    /// and tangent vectors. For cyclic curves, the arrays may be doubled
    /// to allow extracting contiguous slices across the start/end boundary.
    pub struct PointData<'a> {
        /// The input points (may be doubled for cyclic curves).
        /// Note: can't use points.len() directly since this may be doubled.
        pub points: &'a Vec<[f64; DIMS]>,
        /// The actual number of unique points.
        pub points_len: usize,

        /// Cached segment lengths between consecutive points.
        /// This array may be doubled as well for cyclic curves.
        pub points_length_cache: &'a Vec<f64>,

        /// Tangent vectors at each knot (2 per knot: incoming and outgoing).
        pub tangents: &'a Vec<[f64; DIMS]>,
    }

    /// Handle lengths and error for the 2 segments between 3 knots.
    #[derive(Copy, Clone)]
    pub struct KnotAdjacentParams {
        /// Handle lengths [left, right] for the segment before this knot.
        pub handles_prev: [f64; 2],
        /// Handle lengths [left, right] for the segment after this knot.
        pub handles_next: [f64; 2],
        /// Squared error for the segment before this knot.
        pub error_sq_prev: f64,
        /// Squared error for the segment after this knot.
        pub error_sq_next: f64,
    }
}

/// Re-export core types from the types module.
pub use self::types::{
    Knot,
    PointData,
    KnotAdjacentParams,
};

/// Advance knot index to next knot in array, wrapping at end.
#[inline]
fn knot_step_next_wrap(k_step: &mut usize, knots_end: usize) {
    if *k_step != knots_end {
        *k_step += 1;
    } else {
        // Wrap around.
        *k_step = 0;
    }
}

/// Find the knot furthest from the line between `k_prev` and `k_next`.
///
/// This is used to find a good split point when a curve segment
/// exceeds the error threshold. The point with maximum perpendicular
/// distance from the chord is chosen as the split point.
pub fn knot_find_split_point(
    pd: &PointData,
    knots: &Vec<Knot>,
    k_prev: &Knot,
    k_next: &Knot,
) -> usize {
    let mut split_point: usize = INVALID;
    let mut split_point_dist_best: f64 = -::std::f64::MAX;

    let offset = &pd.points[k_prev.index];

    let mut v_plane = sub_vnvn(&pd.points[k_prev.index], &pd.points[k_next.index]);
    normalize_vn(&mut v_plane);

    let knots_end = knots.len() - 1;
    let mut k_step = k_prev.index;
    loop {
        knot_step_next_wrap(&mut k_step, knots_end);

        if k_step != k_next.index {
            let knot = &knots[k_step];
            let v_offset = sub_vnvn(&pd.points[knot.index], offset);
            let v_proj = project_plane_vnvn_normalized(&v_offset, &v_plane);
            let split_point_dist_test = len_squared_vn(&v_proj);
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

/// Find the knot furthest from the line between `k_prev` and `k_next` along a given axis.
///
/// Similar to `knot_find_split_point`, but projects points onto the given
/// plane normal instead of perpendicular to the chord. Used for corner
/// detection to find split points that best separate angled segments.
pub fn knot_find_split_point_on_axis(
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
        knot_step_next_wrap(&mut k_step, knots_end);

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


/// Fit a curve segment and return error metrics.
///
/// Returns (error_sq, error_index, [handle_left, handle_right]).
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

/// Calculate the number of points from `index_l` to `index_r` inclusive.
///
/// Handles cyclic wrap-around when `index_l > index_r`.
#[inline]
fn knot_span_length(index_l: usize, index_r: usize, points_len: usize) -> usize {
    (if index_l <= index_r {
        index_r - index_l
    } else {
        (index_r + points_len) - index_l
    }) + 1
}

/// Calculate the curve fit error between two knots, including the error index.
///
/// Returns (error_sq, error_index, [handle_left, handle_right]).
/// The error_index is the global index of the point with maximum error.
pub fn knot_calc_curve_error_value_and_index(
    pd: &PointData,
    knot_l: &Knot, knot_r: &Knot,
    tan_l: &[f64; DIMS],
    tan_r: &[f64; DIMS],
) -> (f64, usize, [f64; 2]) {
    let points_offset_len = knot_span_length(knot_l.index, knot_r.index, pd.points_len);

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

/// Calculate the curve fit error between two knots (without error index).
///
/// Returns (error_sq, [handle_left, handle_right]).
pub fn knot_calc_curve_error_value(
    pd: &PointData,
    knot_l: &Knot, knot_r: &Knot,
    tan_l: &[f64; DIMS],
    tan_r: &[f64; DIMS],
) -> (f64, [f64; 2]) {
    let points_offset_len = knot_span_length(knot_l.index, knot_r.index, pd.points_len);

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

macro_rules! unlikely { ($body:expr) => { $body } }

/// First refinement pass: iteratively remove knots below the error threshold.
pub mod refine_remove {
    use super::{
        INVALID,
        knot_calc_curve_error_value,
        Knot,
        PointData,
    };
    use ::intern::min_heap;

    /// State stored in the heap for potential knot removal.
    ///
    /// Stores the handle lengths that would result from removing this knot.
    #[derive(Copy, Clone)]
    struct KnotRemoveState {
        /// Index of the knot being considered for removal.
        index: usize,
        /// Handle lengths if this knot is removed.
        handles: [f64; 2],
    }

    /// (Re)calculate the error for removing a knot and update the heap.
    ///
    /// Tests the error that would result from removing k_curr and fitting
    /// a single curve from k_prev to k_next. Updates the heap if the error
    /// is below the threshold.
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
            heap.insert_or_update(
                k_curr_heap_node,
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

    /// Iteratively remove all points under the error threshold.
    ///
    /// Uses a min-heap to efficiently find and remove the knot with the
    /// smallest error value. After each removal, the error values of
    /// adjacent knots are recalculated and the heap is updated.
    ///
    /// This is the first pass of the curve refinement algorithm.
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
                let k_curr: &mut Knot = &mut knots[r.index];

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
            // Remove ourselves.
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


/// Second refinement pass: reposition knots to minimize error.
///
/// Tests moving each knot to positions between its neighbors to find
/// optimal placement. Knots may also be removed if they fall below threshold.
pub mod refine_refit {

    use super::{
        INVALID,
        knot_calc_curve_error_value,
        knot_calc_curve_error_value_and_index,
        knot_find_split_point,
        Knot,
        KnotAdjacentParams,
        PointData,
    };
    use ::intern::min_heap;

    /// Result from refining a refit index in one direction.
    struct RefineResult {
        /// The best refit index found during the search.
        index_refit: usize,
        /// Handle lengths and errors for the refit position.
        params: KnotAdjacentParams,
        /// True if an improvement was found over the initial position.
        is_refined: bool,
    }

    /// Refine the refit index by searching neighbors for lower error.
    /// Stops when no further improvement is found.
    /// `dir`: -1 to search toward `k_prev`, 1 to search toward `k_next`.
    fn knot_refit_index_refine(
        pd: &PointData,
        knots: &Vec<Knot>,
        k_prev: &Knot,
        k_next: &Knot,
        index_refit: usize,
        dir: isize,
        mut cost_sq_max: f64,
    ) -> RefineResult {
        // Stop before reaching the adjacent knot.
        let index_end = if dir == -1 { k_prev.index } else { k_next.index };
        let points_len = pd.points_len;
        let mut i = index_refit;
        let mut result = RefineResult {
            index_refit,
            params: KnotAdjacentParams {
                handles_prev: [0.0; 2],
                handles_next: [0.0; 2],
                error_sq_prev: 0.0,
                error_sq_next: 0.0,
            },
            is_refined: false,
        };

        // Step through indices in direction `dir`, with wraparound.
        loop {
            i = ((i as isize + dir + points_len as isize) as usize) % points_len;
            if i == index_end {
                break;
            }

            let k_test = &knots[i];
            let (error_sq_prev, handles_prev_test) = knot_calc_curve_error_value(
                pd, k_prev, k_test,
                &pd.tangents[k_prev.tan[1]],
                &pd.tangents[k_test.tan[0]],
            );
            if error_sq_prev >= cost_sq_max {
                break;
            }

            let (error_sq_next, handles_next_test) = knot_calc_curve_error_value(
                pd, k_test, k_next,
                &pd.tangents[k_test.tan[1]],
                &pd.tangents[k_next.tan[0]],
            );
            if error_sq_next >= cost_sq_max {
                break;
            }

            // Raise the bar: subsequent iterations must beat this.
            cost_sq_max = error_sq_prev.max(error_sq_next);
            result.index_refit = i;
            result.params.handles_prev = handles_prev_test;
            result.params.handles_next = handles_next_test;
            result.params.error_sq_prev = error_sq_prev;
            result.params.error_sq_next = error_sq_next;
            result.is_refined = true;
        }
        return result;
    }

    /// State stored in the heap for potential knot repositioning.
    #[derive(Copy, Clone)]
    struct KnotRefitState {
        /// Index of the knot being considered for refit.
        index: usize,
        /// Target position index for refitting. When INVALID, remove this knot instead.
        index_refit: usize,
        /// Handle lengths and errors for the refit position.
        fit_params: KnotAdjacentParams,
    }

    /// (Re)calculate the error and optimal refit position for a knot.
    ///
    /// This function finds the best position for `k_curr` between its neighbors
    /// `k_prev` and `k_next`. It tests potential positions and calculates the
    /// resulting curve fit error for each.
    ///
    /// When `use_optimize_exhaustive` is true, all positions between neighbors
    /// are tested. Otherwise, a faster heuristic based on the split point is used.
    fn knot_refit_error_recalculate(
        pd: &PointData,
        heap: &mut min_heap::MinHeap<f64, KnotRefitState>,
        knots: &Vec<Knot>,
        knots_handle: &mut Vec<min_heap::NodeHandle>,
        k_curr: &Knot,
        error_max_sq: f64,
        use_optimize_exhaustive: bool,
        use_refit_remove: bool,
    ) {
        debug_assert!(k_curr.no_remove == false);

        let k_curr_heap_node = &mut knots_handle[k_curr.index];

        let k_prev = &knots[k_curr.prev];
        let k_next = &knots[k_curr.next];

        let mut k_refit_index;

        if use_refit_remove {
            // Support re-fitting to remove points.
            let (fit_error_max_sq, fit_error_index, handles) =
                knot_calc_curve_error_value_and_index(
                    pd, k_prev, k_next,
                    &pd.tangents[k_prev.tan[1]],
                    &pd.tangents[k_next.tan[0]],
                    );

            if fit_error_max_sq < error_max_sq {
                // Always perform removal before refitting, (make a negative number)
                heap.insert_or_update(
                    k_curr_heap_node,
                    // Weight for the greatest improvement.
                    fit_error_max_sq - error_max_sq,
                    KnotRefitState {
                        index: k_curr.index,
                        // INVALID == remove
                        index_refit: INVALID,
                        fit_params: KnotAdjacentParams {
                            handles_prev: [handles[0], 0.0],  // [1] unused
                            handles_next: [0.0, handles[1]],  // [0] unused
                            error_sq_prev: fit_error_max_sq,
                            error_sq_next: fit_error_max_sq,
                        },
                    }
                );
                return;
            }

            // Use the largest point of difference when removing
            // as the target to refit to.
            k_refit_index = fit_error_index;
        } else {
            k_refit_index = knot_find_split_point(pd, knots, k_prev, k_next);
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

        /// Calculate curve errors for both segments around a refit position.
        /// Returns None if either segment exceeds the error threshold.
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
                    pd, k_prev, &knots[k_refit_index], k_next, cost_sq_src_max);

            // Local refinement: search neighbors for a better refit index.
            // Search both directions independently to avoid bias.
            // Skip when error is zero (e.g. exactly straight lines).
            if let Some((_handles_prev, fit_error_dst_prev, _handles_next, fit_error_dst_next)) =
                refit_result_or_none
            {
                let cost_sq_dst_max_init = fit_error_dst_prev.max(fit_error_dst_next);
                if cost_sq_dst_max_init > 0.0 {
                    // Search toward k_prev (dir=-1) and toward k_next (dir=1).
                    let scan_prev = knot_refit_index_refine(
                        pd, knots, k_prev, k_next, k_refit_index, -1, cost_sq_dst_max_init);
                    let scan_next = knot_refit_index_refine(
                        pd, knots, k_prev, k_next, k_refit_index, 1, cost_sq_dst_max_init);

                    // Pick the best result from both directions.
                    if scan_prev.is_refined || scan_next.is_refined {
                        let scan = if scan_prev.is_refined && scan_next.is_refined {
                            // Both directions found improvements, pick the best.
                            // In the unlikely event of a tie, minimum error breaks it.
                            let cost_sq_max_prev = scan_prev.params.error_sq_prev.max(scan_prev.params.error_sq_next);
                            let cost_sq_max_next = scan_next.params.error_sq_prev.max(scan_next.params.error_sq_next);
                            if cost_sq_max_prev < cost_sq_max_next {
                                &scan_prev
                            } else if cost_sq_max_next < cost_sq_max_prev {
                                &scan_next
                            } else {
                                let cost_sq_min_prev = scan_prev.params.error_sq_prev.min(scan_prev.params.error_sq_next);
                                let cost_sq_min_next = scan_next.params.error_sq_prev.min(scan_next.params.error_sq_next);
                                if cost_sq_min_prev <= cost_sq_min_next { &scan_prev } else { &scan_next }
                            }
                        } else if scan_prev.is_refined {
                            &scan_prev
                        } else {
                            &scan_next
                        };

                        // Use results from the winning direction.
                        k_refit_index = scan.index_refit;
                        refit_result_or_none = Some((
                            scan.params.handles_prev, scan.params.error_sq_prev,
                            scan.params.handles_next, scan.params.error_sq_next,
                        ));
                    }
                }
            }
        }
        // end exhaustive test

        if let Some((
            handles_prev, fit_error_dst_prev,
            handles_next, fit_error_dst_next,
        )) = refit_result_or_none {
            let fit_error_dst_max_sq =
                fit_error_dst_prev.max(fit_error_dst_next);
            debug_assert!(fit_error_dst_max_sq < cost_sq_src_max);
            heap.insert_or_update(
                k_curr_heap_node,
                // Weight for the greatest improvement.
                cost_sq_src_max - fit_error_dst_max_sq,
                KnotRefitState {
                    index: k_curr.index,
                    index_refit: k_refit_index,
                    fit_params: KnotAdjacentParams {
                        handles_prev: handles_prev,
                        handles_next: handles_next,
                        error_sq_prev: fit_error_dst_prev,
                        error_sq_next: fit_error_dst_next,
                    },
                }
            );
            return;
        }

        if *k_curr_heap_node != min_heap::NodeHandle::INVALID {
            heap.remove(*k_curr_heap_node);
            *k_curr_heap_node = min_heap::NodeHandle::INVALID;
        }
    }

    /// Re-adjust the curves by re-fitting points.
    ///
    /// Test the error from moving each knot to positions between its adjacent knots.
    /// If a better position is found (lower error), the knot is moved there.
    ///
    /// Parameters:
    /// - `use_optimize_exhaustive`: When true, search all positions between adjacent knots
    ///   for the optimal refit location. When false, use a faster heuristic.
    /// - `use_refit_remove`: When true, remove knots that fall below the error threshold
    ///   during the refit process.
    pub fn curve_incremental_simplify_refit(
        pd: &PointData,
        knots: &mut Vec<Knot>,
        knots_handle: &mut Vec<min_heap::NodeHandle>,
        knots_len_remaining: &mut usize,
        error_max_sq: f64,
        use_optimize_exhaustive: bool,
        use_refit_remove: bool,
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
                    error_max_sq, use_optimize_exhaustive, use_refit_remove);
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
                    k_refit.handles[0] = r.fit_params.handles_prev[1];
                    k_refit.handles[1] = r.fit_params.handles_next[0];
                }

                knots[k_prev_index].handles[1] = r.fit_params.handles_prev[0];
                knots[k_next_index].handles[0] = r.fit_params.handles_next[1];
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

                knots[k_prev_index].fit_error_sq_next = r.fit_params.error_sq_prev;

                *knots_len_remaining -= 1;
            } else {
                // Remove ourselves.
                knots[k_next_index].prev = r.index_refit;
                knots[k_prev_index].next = r.index_refit;

                knots[k_prev_index].fit_error_sq_next = r.fit_params.error_sq_prev;

                let k_refit = &mut knots[r.index_refit];
                k_refit.prev = k_prev_index;
                k_refit.next = k_next_index;

                k_refit.fit_error_sq_next = r.fit_params.error_sq_next;

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
                        error_max_sq, use_optimize_exhaustive, use_refit_remove);
                }
            }
        }

        drop(heap);
    }
}
// end refine_refit

/// Corner detection pass: identify and collapse sharp angle transitions.
///
/// Finds adjacent knots where the tangent angle exceeds the threshold and
/// collapses them into corner knots (where tangents are discontinuous).
pub mod refine_corner {
    use super::{
        INVALID,
        knot_calc_curve_error_value,
        knot_find_split_point_on_axis,
        Knot,
        KnotAdjacentParams,
        PointData,
    };
    use ::intern::math_vector::{
        dot_vnvn,
        len_squared_vnvn,
        project_vnvn_normalized,
        sub_vnvn,
    };
    use ::intern::min_heap;

    /// Result of collapsing a corner.
    #[derive(Copy, Clone)]
    struct KnotCornerState {
        /// Index of the knot being considered for corner collapse.
        index: usize,
        /// Indices of adjacent knots [k_prev, k_next] whose tangents will be
        /// collapsed into this corner. The corner inherits tangents from these neighbors.
        index_pair: [usize; 2],
        /// Handle lengths and errors for the collapsed corner.
        fit_params: KnotAdjacentParams,
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
                    // _must_ be assigned to k_split, later
                    heap.insert_or_update(
                        k_split_heap_node,
                        // Weight for the greatest improvement.
                        fit_error_dst_prev.max(fit_error_dst_next),
                        KnotCornerState {
                            index: k_split.index,
                            index_pair: [k_prev.index, k_next.index],
                            fit_params: KnotAdjacentParams {
                                handles_prev: handles_prev,
                                handles_next: handles_next,
                                error_sq_prev: fit_error_dst_prev,
                                error_sq_next: fit_error_dst_next,
                            },
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

    /// Attempt to collapse close knots into corners.
    ///
    /// This function finds adjacent knots that exceed the angle limit and
    /// collapses them into a single corner knot, or removes them entirely
    /// (depending on their error values).
    ///
    /// A corner is created when two adjacent curve segments meet at a sharp angle.
    /// The knot at this junction is marked as a corner, which prevents its
    /// tangent from being shared between the two segments.
    ///
    /// Parameters:
    /// - `error_max_sq`: Maximum allowed squared error for the curve.
    /// - `error_sq_collapse_max`: Maximum squared error for collapsing adjacent knots.
    /// - `corner_angle`: Angle threshold (in radians) above which knots become corners.
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
                k_prev.handles[1]  = c.fit_params.handles_prev[0];
                tan_prev = k_prev.tan[1];

                debug_assert!(c.fit_params.error_sq_prev <= error_max_sq);
                k_prev.fit_error_sq_next = c.fit_params.error_sq_prev;
            }

            {
                let k_next  = &mut knots[k_next_index];
                k_next.prev = k_split_index;
                tan_next = k_next.tan[0];

                k_next.handles[0] = c.fit_params.handles_next[1];
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
                k_split.handles[0] = c.fit_params.handles_prev[1];
                k_split.handles[1] = c.fit_params.handles_next[0];

                debug_assert!(c.fit_params.error_sq_next <= error_max_sq);
                k_split.fit_error_sq_next = c.fit_params.error_sq_next;
            }

            *knots_len_remaining += 1;
        }


        drop(heap);
    }
}

// end refine_corner
