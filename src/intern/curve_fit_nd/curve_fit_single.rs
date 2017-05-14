
use ::intern::math_vector::{
    len_squared_vnvn,
    len_vnvn,
    sub_vnvn,
    dot_vnvn,
};

// weak?
const DIMS: usize = ::intern::math_vector::DIMS;

mod types {
    use super::{
        DIMS,
    };
    #[derive(Copy, Clone)]
    pub struct Cubic {
        pub p0: [f64; DIMS],
        pub p1: [f64; DIMS],
        pub p2: [f64; DIMS],
        pub p3: [f64; DIMS],
    }
}

mod cubic_solve_fallback {
    use super::{
        types,
        DIMS,
    };
    use ::intern::math_vector::{
        len_vnvn,
        madd_vnvn_fl, msub_vnvn_fl,
    };

    pub fn calc(
        points: &[[f64; DIMS]],
        tan_l: &[f64; DIMS],
        tan_r: &[f64; DIMS],
    ) -> types::Cubic {
        let p0 = &points[0];
        let p3 = &points[points.len() - 1];
        let alpha = len_vnvn(p0, p3) / 3.0;

        return types::Cubic {
            p0: *p0,
            p1: msub_vnvn_fl(p0, tan_l, alpha),
            p2: madd_vnvn_fl(p3, tan_r, alpha),
            p3: *p3,
        };
    }
}

mod cubic_solve_least_square {
    use super::{
        types,
        DIMS,
    };
    use ::intern::math_vector::{
        mul_vn_fl,
        madd_vnvn_fl, msub_vnvn_fl,
        is_almost_zero,
    };


    pub fn calc(
        points: &[[f64; DIMS]],
        tan_l: &[f64; DIMS],
        tan_r: &[f64; DIMS],
        u_prime: &[f64],
    ) -> Option<types::Cubic> {
        let p0 = &points[0];
        let p3 = &points[points.len() - 1];

        let (alpha_l, alpha_r) = {
            let mut x: [f64; 2] = [0.0, 0.0];
            let mut c: [[f64; 2]; 2] = [[0.0, 0.0], [0.0, 0.0]];

            for (pt, u) in points.iter().zip(u_prime) {
                let a: [[f64; DIMS]; 2] = [
                    mul_vn_fl(tan_l, bezier::b1(*u)),
                    mul_vn_fl(tan_r, bezier::b2(*u)),
                ];

                let b0_plus_b1 = bezier::b0_plus_b1(*u);
                let b2_plus_b3 = bezier::b2_plus_b3(*u);

                // inline dot product
                for j in 0..DIMS {
                    let tmp = (pt[j] - (p0[j] * b0_plus_b1)) + (p3[j] * b2_plus_b3);

                    x[0] += a[0][j] * tmp;
                    x[1] += a[1][j] * tmp;

                    c[0][0] += a[0][j] * a[0][j];
                    c[0][1] += a[0][j] * a[1][j];
                    c[1][1] += a[1][j] * a[1][j];
                }

                c[1][0] = c[0][1];
            }

            let det_c0_c1 = {
                let tmp = c[0][0] * c[1][1] - c[0][1] * c[1][0];
                if !is_almost_zero(tmp) {
                    tmp
                } else {
                    c[0][0] * c[1][1] * 10e-12
                }
            };
            let det_c_0x  = x[1]    * c[0][0] - x[0]    * c[0][1];
            let det_x_c1  = x[0]    * c[1][1] - x[1]    * c[0][1];

            // may still divide-by-zero, check below will catch nan values.
            (det_x_c1 / det_c0_c1, det_c_0x / det_c0_c1)
        };

        // flip check to catch nan values.
        if !(alpha_l >= 0.0) ||
           !(alpha_r >= 0.0)
        {
            return None;
        } else {
            return Some(types::Cubic {
                p0: *p0,
                p1: msub_vnvn_fl(p0, tan_l, alpha_l),
                p2: madd_vnvn_fl(p3, tan_r, alpha_r),
                p3: *p3,
            })
        }
    }

    // Bezier multipliers
    mod bezier {
        pub fn b1(u: f64) -> f64 {
            let tmp = 1.0 - u;
            return 3.0 * u * tmp * tmp;
        }

        pub fn b2(u: f64) -> f64 {
            return 3.0 * u * u * (1.0 - u);
        }

        pub fn b0_plus_b1(u: f64) -> f64 {
            let tmp = 1.0 - u;
            return tmp * tmp * (1.0 + 2.0 * u);
        }

        pub fn b2_plus_b3(u: f64) -> f64 {
            return u * u * (3.0 - 2.0 * u);
        }
    }

}

mod cubic_solve_circle {
    use super::{
        types,
        DIMS,
    };
    use ::intern::math_vector::{
        len_vnvn,
        len_negated_vnvn,
        dot_vnvn,
        madd_vnvn_fl, msub_vnvn_fl,
    };

    pub fn calc(
        points: &[[f64; DIMS]],
        tan_l: &[f64; DIMS],
        tan_r: &[f64; DIMS],
        points_coords_length: f64,
    ) -> Option<types::Cubic> {
        let p0 = &points[0];
        let p3 = &points[points.len() - 1];

        if let Some(alpha) = points_calc_cubic_scale(p0, p3, tan_l, tan_r, points_coords_length) {
            return Some(types::Cubic {
                p0: *p0,
                p1: msub_vnvn_fl(p0, tan_l, alpha),
                p2: madd_vnvn_fl(p3, tan_r, alpha),
                p3: *p3,
            });
        } else {
            return None;
        }
    }


    // Return a scale value, used to calculate how much the curve handles should be increased,
    //
    // This works by placing each end-point on an imaginary circle,
    // the placement on the circle is based on the tangent vectors,
    // where larger differences in tangent angle cover a larger part of the circle.
    //
    // Return the scale representing how much larger the distance around the circle is.

    fn points_calc_circumference_factor(
        tan_l: &[f64; DIMS],
        tan_r: &[f64; DIMS],
    ) -> f64 {
        use std::f64;
        let dot = dot_vnvn(tan_l, tan_r);

        let len_tangent = if dot < 0.0 { len_vnvn } else { len_negated_vnvn } (tan_l, tan_r);
        if len_tangent > f64::EPSILON {
            // only clamp to avoid precision error.
            let angle = ((-dot.abs()).max(-1.0)).acos();
            // Angle may be less than the length when the
            // tangents define >180 degrees of the circle,
            // (tangents that point away from each other).
            //
            // We could try support this but will likely cause
            // extreme >1 scales which could cause other issues.

            // assert(angle >= len_tangent);
            let factor = angle / len_tangent;
            debug_assert!(factor < (f64::consts::PI / 2.0) + (f64::EPSILON * 10.0));
            return factor;
        } else {
            // tangents are exactly aligned (think two opposite sides of a circle).
            return f64::consts::PI / 2.0;
        }
    }

    // Return the value which the distance between points will need to be scaled by,
    // to define a handle, given both points are on a perfect circle.
    //
    // Note: the return value will need to be multiplied by 1.3... for correct results.
    fn points_calc_circle_tangent_factor(
        tan_l: &[f64; DIMS],
        tan_r: &[f64; DIMS],
    ) -> Option<f64> {
        let eps = 1e-8;
        let tan_dot = dot_vnvn(tan_l, tan_r);
        if tan_dot > 1.0 - eps {
            // no angle difference (use fallback, length wont make any difference)
            return None;
        } else if tan_dot < -1.0 + eps {
            // parallel tangents (half-circle)
            return Some(1.0 / 2.0);
        } else {
            // non-aligned tangents, calculate handle length
            let angle = tan_dot.acos() / 2.0;

            // could also use 'angle_sin = len_vnvn(tan_l, tan_r) / 2.0'
            let angle_sin = angle.sin();
            let angle_cos = angle.cos();
            return Some(((1.0 - angle_cos) / (angle_sin * 2.0)) / angle_sin);
        }
    }

    // Calculate the scale the handles, which serves as a best-guess
    // used as a fallback when the least-square solution fails.
    fn points_calc_cubic_scale(
        v_l: &[f64; DIMS],
        v_r: &[f64; DIMS],
        tan_l: &[f64; DIMS],
        tan_r: &[f64; DIMS],
        coords_length: f64,
    ) -> Option<f64> {

        if let Some(len_circle_factor) = points_calc_circle_tangent_factor(tan_l, tan_r) {

            let len_direct = len_vnvn(v_l, v_r);

            // if this curve is a circle, this value doesn't need modification
            let len_circle_handle = len_direct * (len_circle_factor / 0.75);

            // scale by the difference from the circumference distance
            let len_circle = len_direct * points_calc_circumference_factor(tan_l, tan_r);
            let mut scale_handle = coords_length / len_circle;

            // Could investigate an accurate calculation here,
            // though this gives close results.
            scale_handle = ((scale_handle - 1.0) * 1.75) + 1.0;

            scale_handle *= len_circle_handle;

            if scale_handle.is_finite() {
                return Some(scale_handle);
            }
        }
        return None;
    }
}

mod cubic_solve_offset {
    use super::{
        types,
        DIMS,
    };
    use ::intern::math_vector::{
        sub_vnvn,
        dot_vnvn,
        madd_vnvn_fl, msub_vnvn_fl,
        negated_vn,
        normalized_vnvn,
        normalized_vn,
        project_plane_vnvn_normalized,
        project_vnvn_normalized,
    };

    pub fn calc(
        points: &[[f64; DIMS]],
        tan_l: &[f64; DIMS],
        tan_r: &[f64; DIMS],
    ) -> Option<types::Cubic> {
        use std::f64;

        let p0 = &points[0];
        let p3 = &points[points.len() - 1];

        let dir_unit = normalized_vnvn(p3, p0);
        // note that normalizing output here is only for better accuracy, not essential.
        let a: [[f64; DIMS]; 2] = [
                        normalized_vn(&project_plane_vnvn_normalized(tan_l, &dir_unit)),
            negated_vn(&normalized_vn(&project_plane_vnvn_normalized(tan_r, &dir_unit))),
        ];

        let mut dists: [f64; 2] = [0.0, 0.0];

        // early exit to avoid unnecessary calculation & divide-by-zero.
        let div_l = dot_vnvn(tan_l, &a[0]).abs();
        let div_r = dot_vnvn(tan_r, &a[1]).abs();

        if (div_l < f64::EPSILON) ||
           (div_r < f64::EPSILON)
        {
            return None;
        }

        for pt in &points[1..(points.len() - 1)] {
            for k in 0..2 {
                let tmp = project_vnvn_normalized(&sub_vnvn(p0, pt), &a[k]);
                dists[k] = dists[k].max(dot_vnvn(&tmp, &a[k]));
            }
        }

        let alpha_l = (dists[0] / 0.75) / div_l;
        let alpha_r = (dists[1] / 0.75) / div_r;

        if !(alpha_l >= 0.0) ||
           !(alpha_r >= 0.0)
        {
            return None;
        } else {
            return Some(types::Cubic {
                p0: *p0,
                p1: msub_vnvn_fl(p0, tan_l, alpha_l),
                p2: madd_vnvn_fl(p3, tan_r, alpha_r),
                p3: *p3,
            });
        }
    }
}


/// Use Newton-Raphson iteration to find better root.
///
/// * `cubic` - Current fitted curve.
/// * `p` - Point to test against.
/// * `u` - Parameter value for `p`.
///
/// Note: return value may be `nan` caller must check for this.
fn cubic_find_root(
    cubic: &types::Cubic,
    p: &[f64; DIMS],
    u: f64,
) -> f64 {
    // Newton-Raphson Method.
    // all vectors
    let q0_u = sub_vnvn(&cubic_calc_point(cubic, u), p);
    let q1_u = cubic_calc_speed(cubic, u);
    let q2_u = cubic_calc_acceleration(cubic, u);

    // may divide-by-zero, caller must check for that case.

    // u - (q0_u * q1_u) / (q1_u.length_squared() + q0_u * q2_u)
    return u - dot_vnvn(&q0_u, &q1_u) / (dot_vnvn(&q1_u, &q1_u) + dot_vnvn(&q0_u, &q2_u));
}

/// Given set of points and their parameterization, try to find a better parameterization.
fn cubic_reparameterize(
    cubic: &types::Cubic,
    points: &[[f64; DIMS]],
    u_prime_src: &[f64],

    u_prime_dst: &mut [f64]
) -> bool {
    debug_assert!(points.len() == u_prime_src.len());
    debug_assert!(points.len() == u_prime_dst.len());

    // Recalculate the values of u[] based on the Newton Raphson method.
    for ((u_src, u_dst), pt) in u_prime_src.iter().zip(&mut *u_prime_dst).zip(points) {
        *u_dst = cubic_find_root(cubic, pt, *u_src);
        if !(*u_dst).is_finite() {
            return false;
        }
    }

    // we can safely unwrap here because nan/inf's are caught above
    u_prime_dst.sort_by(|a, b| a.partial_cmp(b).unwrap());

    if (u_prime_dst[0] < 0.0) ||
       (u_prime_dst[points.len() - 1] > 1.0)
    {
        return false;
    }

    debug_assert!(u_prime_dst[0] >= 0.0);
    debug_assert!(u_prime_dst[u_prime_dst.len() - 1] <= 1.0);
    return true;
}

fn points_calc_coord_length(
    points: &[[f64; DIMS]],
    points_length_cache: &[f64],
) -> (Vec<f64>, f64) {
    let mut u: Vec<f64> = Vec::with_capacity(points.len());
    u.push(0.0);

    let mut pt_prev = &points[0];
    let mut l_prev = 0.0;
    for (pt, l) in points.iter().zip(points_length_cache).skip(1) {
        debug_assert!(len_vnvn(pt, pt_prev) == *l);
        let l_curr = l + l_prev;
        u.push(l_curr);

        pt_prev = pt;
        l_prev = l_curr;
    }

    debug_assert!(u.len() == points.len());

    let w = u[u.len() - 1];
    for u_step in &mut u[1..] {
        *u_step /= w;
    }

    return (u, w);
}

fn cubic_calc_point(
    cubic: &types::Cubic, t: f64,
) -> [f64; DIMS] {
    let p0 = &cubic.p0;
    let p1 = &cubic.p1;
    let p2 = &cubic.p2;
    let p3 = &cubic.p3;
    let s = 1.0 - t;
    let mut v_out = [0.0; DIMS];
    for j in 0..DIMS {
        let p01 = (p0[j] * s) + (p1[j] * t);
        let p12 = (p1[j] * s) + (p2[j] * t);
        let p23 = (p2[j] * s) + (p3[j] * t);
        v_out[j] = ((((p01 * s) + (p12 * t))) * s) +
                   ((((p12 * s) + (p23 * t))) * t);
    }
    return v_out;
}

fn cubic_calc_speed(
    cubic: &types::Cubic, t: f64,
) -> [f64; DIMS] {
    let p0 = &cubic.p0;
    let p1 = &cubic.p1;
    let p2 = &cubic.p2;
    let p3 = &cubic.p3;
    let s = 1.0 - t;
    let mut v_out = [0.0; DIMS];
    for j in 0..DIMS {
        v_out[j] =  3.0 * ((p1[j] - p0[j]) * s * s + 2.0 *
                           (p2[j] - p0[j]) * s * t +
                           (p3[j] - p2[j]) * t * t);
    }
    return v_out;
}

fn cubic_calc_acceleration(
    cubic: &types::Cubic, t: f64,
) -> [f64; DIMS] {
    let p0 = &cubic.p0;
    let p1 = &cubic.p1;
    let p2 = &cubic.p2;
    let p3 = &cubic.p3;
    let s = 1.0 - t;
    let mut v_out = [0.0; DIMS];
    for j in 0..DIMS {
        v_out[j] = 6.0 * ((p2[j] - 2.0 * p1[j] + p0[j]) * s +
                          (p3[j] - 2.0 * p2[j] + p1[j]) * t);
    }
    return v_out;
}

#[derive(Clone, Copy)]
struct FitError {
    pub max_sq: f64,
    pub index: usize,
}

fn cubic_calc_error(
    cubic: &types::Cubic,
    points: &[[f64; DIMS]],
    u: &[f64],
) -> FitError {
    let mut error_max_sq = -1.0;

    // no need to measure first & last points
    let skip_endpoints = 1..(points.len() - 1);
    let mut index = 1;
    let mut error_index = 1;
    for (pt_real, u_step) in
        points[skip_endpoints.clone()].iter().zip(
            &u[skip_endpoints.clone()])
    {
        let pt_eval = cubic_calc_point(cubic, *u_step);
        let err_sq = len_squared_vnvn(pt_real, &pt_eval);
        if err_sq > error_max_sq {
            error_max_sq = err_sq;
            error_index = index;
        }
        index += 1;
    }

    debug_assert!(error_max_sq != -1.0);
    return FitError {
        max_sq: error_max_sq,
        index: error_index,
    };
}

/// Like `cubic_calc_error` but return None
/// in the case we can't improve on `error_max_sq_limit`.
fn cubic_calc_error_limit(
    cubic: &types::Cubic,
    points: &[[f64; DIMS]],
    u: &[f64],
    error_max_sq_limit: f64,
) -> Option<FitError> {
    let mut error_max_sq = -1.0;

    // no need to measure first & last points
    let skip_endpoints = 1..(points.len() - 1);
    let mut index = 1;
    let mut error_index = 1;
    for (pt_real, u_step) in
        points[skip_endpoints.clone()].iter().zip(
            &u[skip_endpoints.clone()])
    {
        let pt_eval = cubic_calc_point(cubic, *u_step);
        let err_sq = len_squared_vnvn(pt_real, &pt_eval);
        if err_sq > error_max_sq {
            if err_sq > error_max_sq_limit {
                return None;
            }
            error_max_sq = err_sq;
            error_index = index;
        }
        index += 1;
    }

    // println!("~~ {}", points.len());
    debug_assert!(error_max_sq != -1.0);
    return Some(FitError {
        max_sq: error_max_sq,
        index: error_index,
    });
}

fn fit_cubic_to_points(
    points: &[[f64; DIMS]],
    points_length_cache: &[f64],
    tan_l: &[f64; DIMS],
    tan_r: &[f64; DIMS],
) -> (types::Cubic, FitError) {
    let iteration_max = 4;

    assert!(points.len() > 2);

    let cubic_fallback = cubic_solve_fallback::calc(points, tan_l, tan_r);

    let (mut u, points_length) = points_calc_coord_length(points, points_length_cache);
    let error_fallback = cubic_calc_error(&cubic_fallback, points, &u);
    let mut error_best = error_fallback;
    let mut cubic_best = cubic_fallback;

    macro_rules! cubic_test_error {
        ($cubic_test:expr) => {
            {
                let error_test = cubic_calc_error(
                    $cubic_test, points, &u);
                if error_best.max_sq > error_test.max_sq {
                    cubic_best = *$cubic_test;
                    error_best = error_test;
                }
                error_test
            }
        }
    }

    macro_rules! cubic_test_error_limit {
        ($cubic_test:expr) => {
            {
                if let Some(error_test) = cubic_calc_error_limit(
                    $cubic_test, points, &u, error_best.max_sq)
                {
                    cubic_best = *$cubic_test;
                    error_best = error_test;
                }
            }
        }
    }

    if let Some(cubic_test) = cubic_solve_circle::calc(points, tan_l, tan_r, points_length) {
        cubic_test_error_limit!(&cubic_test);
    }

    if let Some(cubic_test) = cubic_solve_offset::calc(points, tan_l, tan_r) {
        cubic_test_error_limit!(&cubic_test);
    }

    {
        let mut cubic_least_square;
        let mut error_least_square;

        if let Some(cubic_test) = cubic_solve_least_square::calc(points, tan_l, tan_r, &u) {
            // we want the result so we can refine it (even if its currently not the best)
            error_least_square = cubic_test_error!(&cubic_test);
            cubic_least_square = cubic_test;
        } else {
            error_least_square = error_fallback;
            cubic_least_square = cubic_fallback;
        }

        let mut u_prime: Vec<f64> = vec![0.0; u.len()];
        for _ in 0..iteration_max {
            if !cubic_reparameterize(&cubic_least_square, points, &u, &mut u_prime) {
                break;
            }

            if let Some(cubic_test) =
                cubic_solve_least_square::calc(points, tan_l, tan_r, &u_prime)
            {
                let error_test = cubic_calc_error(&cubic_test, points, &u_prime);

                if error_least_square.max_sq > error_test.max_sq {
                    error_least_square = error_test;
                    cubic_least_square = cubic_test;
                } else {
                    // break if we're getting worse
                    // break;
                }
                ::std::mem::swap(&mut u, &mut u_prime);
            } else {
                break;
            }

        }
        drop(u_prime);
        drop(u);

        if error_best.max_sq > error_least_square.max_sq {
            error_best = error_least_square;
            cubic_best = cubic_least_square;
        }
    }

    return (cubic_best, error_best);
}

//
// Return error squared, and both handle locations
//
pub fn curve_fit_cubic_to_points_single(
    points: &[[f64; DIMS]],
    points_length_cache: &[f64],
    tan_l: &[f64; DIMS],
    tan_r: &[f64; DIMS],
) -> ((f64, usize), [f64; DIMS], [f64; DIMS]) {
    let (cubic, fit_error) = fit_cubic_to_points(
        points,
        points_length_cache,
        tan_l, tan_r);

    return ((fit_error.max_sq, fit_error.index), cubic.p1, cubic.p2);
}
