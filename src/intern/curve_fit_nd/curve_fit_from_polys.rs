
///
/// Perform cubic curve fitting
///
/// This module takes a complete polygon and optimizes curve fitting
/// and optionally corner calculation,
/// outputting a bezier curve that fits within an error margin.
///

/// Enable the refit pass for improved curve quality.
const USE_REFIT: bool = true;
/// Allow removing knots during the refit pass.
const USE_REFIT_REMOVE: bool = true;
/// Scale factor for corner detection error threshold.
const CORNER_SCALE: f64 = 2.0;  // this is weak, should be made configurable.

use ::intern::math_vector::{
    add_vnvn,
    copy_vnvn,
    madd_vnvn_fl,
    normalize_vn,
    normalized_vnvn_with_len,
    sq,
    zero_vn,
};

use ::intern::min_heap;

/// Number of dimensions for curve fitting (2D points).
const DIMS: usize = ::intern::math_vector::DIMS;

use std::collections::LinkedList;

// Import refit module types and functions
use super::curve_fit_cubic_refit::{
    INVALID,
    Knot,
    PointData,
    refine_remove,
    refine_refit,
    refine_corner,
};

/// Tracing mode for polygon extraction from raster images.
#[derive(Copy, Clone, PartialEq)]
pub enum TraceMode {
    /// Trace the outline (boundary) of shapes.
    Outline,
    /// Trace the centerline (skeleton) of shapes.
    Centerline,
}


/// Fit cubic bezier curves to a single polygon.
///
/// Takes a list of points and fits bezier curves that approximate the polygon
/// within the specified error threshold. Optionally detects corners where the
/// curve direction changes sharply.
///
/// Returns a list of cubic bezier segments, where each segment is represented
/// as `[handle_in, point, handle_out]`.
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

    // `curve_incremental_simplify_refit` can be called here, but it's very slow,
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
            sq(error_threshold), use_optimize_exhaustive, USE_REFIT_REMOVE);
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


/// Fit cubic bezier curves to a list of polygons.
///
/// Processes multiple polygons in parallel when there are more than one.
/// Each polygon is processed independently using `fit_poly_single`.
///
/// The input is a list of `(is_cyclic, points)` tuples.
/// Returns a list of `(is_cyclic, bezier_segments)` tuples.
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
