//! Curve Fitting Module
//!
//! Provides algorithms for fitting cubic bezier curves to point data.
//! This includes single-segment fitting, polygon fitting, and iterative
//! refinement with corner detection.

use std::collections::LinkedList;
use crate::intern::math_vector::DIMS;

/// Polygon to bezier curve conversion with parallel processing support.
mod curve_fit_from_polys;
/// Single-segment cubic bezier fitting with multiple solver strategies.
mod curve_fit_cubic;
/// Iterative curve refinement: knot removal, repositioning, and corner detection.
mod curve_fit_cubic_refit;

/// Re-export math_vector for external access to vector operations.
#[allow(unused_imports)]
pub use crate::intern::math_vector;

/// Re-export TraceMode from curve_fit_from_polys.
#[allow(unused_imports)]
pub use self::curve_fit_from_polys::TraceMode;

/// Re-export refinement types for external use.
#[allow(unused_imports)]
pub use self::curve_fit_cubic_refit::{
    Knot,
    PointData,
    refine_remove,
    refine_refit,
    refine_corner,
};

// ============================================================================
// Wrapper functions for raster-retrace compatibility
// ============================================================================

/// Fit cubic bezier curves to a single polygon (raster-retrace wrapper).
///
/// Converts from fixed-size array format to flat array format and back.
pub fn fit_poly_single(
    poly_src: &Vec<[f64; DIMS]>,
    is_cyclic: bool,
    error_threshold: f64,
    corner_angle: f64,
    use_optimize_exhaustive: bool,
) -> Vec<[[f64; DIMS]; 3]> {
    // Handle edge cases: curves with <= 2 points can't be fit
    if poly_src.len() < 2 {
        return Vec::new();
    }
    if poly_src.len() == 2 {
        // For 2-point curves, return a simple line segment as a degenerate cubic
        let mut result: Vec<[[f64; DIMS]; 3]> = Vec::with_capacity(2);
        // First point: handles point towards second point
        let mut seg0 = [[0.0; DIMS]; 3];
        let mut seg1 = [[0.0; DIMS]; 3];
        for k in 0..DIMS {
            let p0 = poly_src[0][k];
            let p1 = poly_src[1][k];
            let delta = (p1 - p0) / 3.0;
            // First segment
            seg0[0][k] = p0 - delta;  // handle_in (mirrored)
            seg0[1][k] = p0;          // point
            seg0[2][k] = p0 + delta;  // handle_out
            // Second segment
            seg1[0][k] = p1 - delta;  // handle_in
            seg1[1][k] = p1;          // point
            seg1[2][k] = p1 + delta;  // handle_out (mirrored)
        }
        result.push(seg0);
        result.push(seg1);
        return result;
    }

    // Convert to flat array
    let mut points_flat: Vec<f64> = Vec::with_capacity(poly_src.len() * DIMS);
    for point in poly_src {
        points_flat.extend_from_slice(point);
    }

    // Call the upstream function
    let (cubic_flat, _orig_indices) = curve_fit_from_polys::fit_poly_single(
        &points_flat,
        DIMS,
        is_cyclic,
        error_threshold,
        corner_angle,
        use_optimize_exhaustive,
    );

    // Convert back to fixed-size array format
    // Each segment has 3 points (handle_in, point, handle_out), each with DIMS components
    let num_segments = cubic_flat.len() / (3 * DIMS);
    let mut result: Vec<[[f64; DIMS]; 3]> = Vec::with_capacity(num_segments);

    for i in 0..num_segments {
        let base = i * 3 * DIMS;
        let mut segment = [[0.0; DIMS]; 3];
        for j in 0..3 {
            for k in 0..DIMS {
                segment[j][k] = cubic_flat[base + j * DIMS + k];
            }
        }
        result.push(segment);
    }

    result
}

/// Fit cubic bezier curves to a list of polygons with parallel processing.
pub fn fit_poly_list(
    poly_list_src: LinkedList<(bool, Vec<[f64; DIMS]>)>,
    error_threshold: f64,
    corner_angle: f64,
    use_optimize_exhaustive: bool,
) -> LinkedList<(bool, Vec<[[f64; DIMS]; 3]>)> {
    let mut curve_list_dst: LinkedList<(bool, Vec<[[f64; DIMS]; 3]>)> = LinkedList::new();

    // Single threaded for small lists.
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

        // Sort by length for more even threading.
        // Larger polygons at the end so they are popped and handled first.
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

    curve_list_dst
}

