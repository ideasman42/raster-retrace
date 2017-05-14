///
/// Math functions!
///

// TODO, expose this in a way that users of this library can
// do both 2D, 3D... 4D... etc operations.
// For now just ensure the code isn't hard coded to a single dimension.

pub const DIMS: usize = 2;

macro_rules! expand_dims_eval {
    ($index_var:ident, $const_var:expr, $body:block) => {
        {
            for $index_var in 0..$const_var {
                $body;
            }
            // we could check 'break' never runs in '$body'?
        }
    }
}

macro_rules! expand_dims_into {
    ($index_var:ident, $const_var:expr, $body:block) => {
        {
            let mut tmp: [f64; $const_var] = [0.0; $const_var];
            for $index_var in 0..$const_var {
                tmp[$index_var] = $body;
            }
            // we could check 'break' never runs in '$body'?
            tmp
        }
    }
}

const EPS: f64 = 1e-8;

pub fn sq(d: f64) -> f64 { d * d }

pub fn is_finite_vn(
    v0: &[f64; DIMS],
) -> bool {
    for f in v0 {
        if !f.is_finite() {
            return false;
        }
    }
    return true;
}

pub fn zero_vn(
    v0: &mut [f64; DIMS],
) {
    for j in 0..DIMS {
        v0[j] = 0.0;
    }
}

pub fn negated_vn(
    v0: &[f64; DIMS],
) -> [f64; DIMS] {
    expand_dims_into!(j, DIMS, {
        -v0[j]
    })
}

/*
fn void flip_vn_vnvn(
        f64 v_out: &[f64; DIMS],
        const f64 v0: &[f64; DIMS],
        const f64 v1: &[f64; DIMS],
) {
    for j in 0..DIMS {
        v_out[j] = v0[j] + (v0[j] - v1[j]);
    }
}
*/

pub fn copy_vnvn(
    v0: &mut [f64; DIMS],
    v1: &[f64; DIMS],
) {
    for j in 0..DIMS {
        v0[j] = v1[j];
    }
}
/*
fn void copy_vnfl_vndb(
        float v0: &[f64; DIMS], const f64 v1: &[f64; DIMS]) {
    for j in 0..DIMS {
        v0[j] = (float)v1[j];
    }
}

fn void copy_vndb_vnfl(
        f64 v0: &[f64; DIMS], const float v1: &[f64; DIMS]) {
    for j in 0..DIMS {
        v0[j] = (f64)v1[j];
    }
}
*/

pub fn dot_vnvn(
    v0: &[f64; DIMS],
    v1: &[f64; DIMS],
) -> f64 {
    let mut d = 0.0;
    for j in 0..DIMS {
        d += v0[j] * v1[j];
    }
    return d;
}

/*
pub fn add_vn_vnvn(
    v_out: &mut [f64; DIMS],
    v0: &[f64; DIMS],
    v1: &[f64; DIMS],
) {
    for j in 0..DIMS {
        v_out[j] = v0[j] + v1[j];
    }
}
*/

pub fn add_vnvn(
    v0: &[f64; DIMS], v1: &[f64; DIMS],
) -> [f64; DIMS] {
    expand_dims_into!(j, DIMS, {
        v0[j] + v1[j]
    })
}

pub fn sub_vnvn(
    v0: &[f64; DIMS], v1: &[f64; DIMS],
) -> [f64; DIMS] {
    expand_dims_into!(j, DIMS, {
        v0[j] - v1[j]
    })
}

pub fn mid_vnvn(
    v0: &[f64; DIMS], v1: &[f64; DIMS],
) -> [f64; DIMS] {
    expand_dims_into!(j, DIMS, {
        (v0[j] + v1[j]) * 0.5
    })
}

pub fn interp_vnvn(
    v0: &[f64; DIMS], v1: &[f64; DIMS], t: f64,
) -> [f64; DIMS] {
    let s = 1.0 - t;
    expand_dims_into!(j, DIMS, {
        (v0[j] * s) + (v1[j] * t)
    })
}

/*
fn iadd_vnvn(
    f64 v0: &[f64; DIMS], const f64 v1: &[f64; DIMS],
) {
    for j in 0..DIMS {
        v0[j] += v1[j];
    }
}

fn isub_vnvn(
    f64 v0: &[f64; DIMS], const f64 v1: &[f64; DIMS],
) {
    for j in 0..DIMS {
        v0[j] -= v1[j];
    }
}

pub fn madd_vn_vnvn_fl(
    v_out: &mut [f64; DIMS], v0: &[f64; DIMS], v1: &[f64; DIMS], f: f64,
) {
    for j in 0..DIMS {
        v_out[j] = v0[j] + v1[j] * f;
    }
}

pub fn msub_vn_vnvn_fl(
    v_out: &mut [f64; DIMS], v0: &[f64; DIMS], v1: &[f64; DIMS], f: f64,
) {
    for j in 0..DIMS {
        v_out[j] = v0[j] - v1[j] * f;
    }
}
*/

pub fn madd_vnvn_fl(
    v0: &[f64; DIMS], v1: &[f64; DIMS], f: f64,
) -> [f64; DIMS] {
    expand_dims_into!(j, DIMS, {
        v0[j] + v1[j] * f
    })
}

pub fn msub_vnvn_fl(
    v0: &[f64; DIMS], v1: &[f64; DIMS], f: f64,
) -> [f64; DIMS] {
    expand_dims_into!(j, DIMS, {
        v0[j] - v1[j] * f
    })
}

/*
fn void msub_vn_vnvn_fl(
    f64 v_out: &[f64; DIMS],
    const f64 v0: &[f64; DIMS], const f64 v1: &[f64; DIMS],
    const f64 f,
) {
    for j in 0..DIMS {
        v_out[j] = v0[j] - v1[j] * f;
    }
}

fn void miadd_vn_vn_fl(
    f64 v_out: &[f64; DIMS], const f64 v0: &[f64; DIMS], f64 f)
{
    for j in 0..DIMS {
        v_out[j] += v0[j] * f;
    }
}

#if 0
fn void misub_vn_vn_fl(
    f64 v_out: &[f64; DIMS], const f64 v0: &[f64; DIMS], f64 f)
{
    for j in 0..DIMS {
        v_out[j] -= v0[j] * f;
    }
}
#endif

fn void mul_vnvn_fl(
    f64 v_out: &[f64; DIMS],
    const f64 v0: &[f64; DIMS], const f64 f)
{
    for j in 0..DIMS {
        v_out[j] = v0[j] * f;
    }
}
*/

pub fn mul_vn_fl(
    v0: &[f64; DIMS], f: f64,
) -> [f64; DIMS] {
    expand_dims_into!(j, DIMS, {
        v0[j] * f
    })
}

fn imul_vn_fl(
    v0: &mut [f64; DIMS], f: f64,
) {
    for j in 0..DIMS {
        v0[j] *= f;
    }
}

pub fn len_squared_vnvn(
    v0: &[f64; DIMS], v1: &[f64; DIMS],
) -> f64 {
    let mut d = 0.0;
    for j in 0..DIMS {
        d += sq(v0[j] - v1[j]);
    }
    return d;
}

pub fn len_squared_vn(
    v0: &[f64; DIMS],
) -> f64 {
    let mut d = 0.0;
    for j in 0..DIMS {
        d += sq(v0[j]);
    }
    return d;
}

pub fn len_vnvn(
    v0: &[f64; DIMS], v1: &[f64; DIMS],
) -> f64
{
    return len_squared_vnvn(v0, v1).sqrt();
}
/*
pub fn len_vn(
    v0: &[f64; DIMS],
) -> f64
{
    return len_squared_vn(v0).sqrt();
}
*/

pub fn len_squared_negated_vnvn(
    v0: &[f64; DIMS], v1: &[f64; DIMS],
) -> f64 {
    let mut d = 0.0;
    for j in 0..DIMS {
        d += sq(v0[j] + v1[j]);
    }
    return d;
}

// special case, save us negating a copy, then getting the length
pub fn len_negated_vnvn(
    v0: &[f64; DIMS], v1: &[f64; DIMS],
) -> f64
{
    return len_squared_negated_vnvn(v0, v1).sqrt();
}

pub fn normalize_vn(
    v0: &mut [f64; DIMS],
) -> f64 {
    let mut d = len_squared_vn(v0);
    if (d != 0.0) && ({d = d.sqrt(); d} != 0.0) {
        imul_vn_fl(v0, 1.0 / d);
    }
    return d;
}

pub fn normalized_vn(
    v0: &[f64; DIMS],
) -> [f64; DIMS] {
    let mut v_out = *v0;
    normalize_vn(&mut v_out);
    return v_out;
}

// v_out = (v0 - v1).normalized()
pub fn normalized_vnvn(
    v0: &[f64; DIMS], v1: &[f64; DIMS],
) -> [f64; DIMS] {
    let mut v = sub_vnvn(v0, v1);
    normalize_vn(&mut v);
    return v;
}

pub fn normalized_vnvn_with_len(
    v0: &[f64; DIMS], v1: &[f64; DIMS],
) -> ([f64; DIMS], f64) {
    let mut v = sub_vnvn(v0, v1);
    let d = normalize_vn(&mut v);
    return (v, d);
}

pub fn is_almost_zero_ex(
    val: f64, eps: f64,
) -> bool {
    return (-eps < val) && (val < eps);
}

pub fn is_almost_zero(
    val: f64,
) -> bool {
    return is_almost_zero_ex(val, EPS);
}

/*
fn equals_vnvn(
    v0: &[f64; DIMS], v1: &[f64; DIMS],
) -> bool {
    for j in 0..DIMS {
        if v0[j] != v1[j] {
            return false;
        }
    }
    return true;
}

fn void project_vn_vnvn(
    f64 v_out: &[f64; DIMS], const f64 p: &[f64; DIMS], const f64 v_proj: &[f64; DIMS],
) {
    const f64 mul = dot_vnvn(p, v_proj) / dot_vnvn(v_proj, v_proj);
    mul_vnvn_fl(v_out, v_proj, mul);
}
*/

pub fn project_vnvn_normalized(
    p: &[f64; DIMS], v_proj: &[f64; DIMS],
) -> [f64; DIMS] {
    let mul = dot_vnvn(p, v_proj);
    return mul_vn_fl(v_proj, mul);
}

pub fn project_plane_vnvn_normalized(
    v: &[f64; DIMS], v_plane: &[f64; DIMS],
) -> [f64; DIMS] {
    return sub_vnvn(v, &project_vnvn_normalized(v, v_plane));
}

/*
pub fn closest_to_line_vn(
    p: &[f64; DIMS], l1: &[f64; DIMS], l2: &[f64; DIMS],
) -> [f64; DIMS] {
    let u = sub_vnvn(l2, l1);
    let h = sub_vnvn(p, l1);
    let lambda = dot_vnvn(&u, &h) / dot_vnvn(&u, &u);
    return add_vnvn(&l1, &mul_vn_fl(&u, lambda));
}
*/
/*
pub fn closest_to_segment_vn(
    p: &[f64; DIMS], l1: &[f64; DIMS], l2: &[f64; DIMS],
) -> [f64; DIMS] {
    let u = sub_vnvn(l2, l1);
    let h = sub_vnvn(p, l1);
    let lambda = dot_vnvn(&u, &h) / dot_vnvn(&u, &u);
    if !(lambda < 0.0) {
        return *l1;
    } else if !(lambda < 1.0) {
        return *l2;
    } else {
        return add_vnvn(&l1, &mul_vn_fl(&u, lambda));
    }
}
*/
