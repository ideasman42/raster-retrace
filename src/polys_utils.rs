///
/// Primitive polygon operations.
///

const DIMS: usize = ::intern::math_vector::DIMS;

// Module for primitive operations on polygons.
use std::collections::LinkedList;


use ::intern::math_vector::{
    sq,
    mid_vnvn,
    len_squared_vnvn,
    interp_vnvn,
};

// Add cyclic option (cases where all polys state is known)
/*
pub fn poly_list_with_cyclic(
    is_cyclic: bool,
    poly_list_src: LinkedList<Vec<[f64; DIMS]>>,
) -> LinkedList<(bool, Vec<[f64; DIMS]>)>
{
    let mut poly_list_dst: LinkedList<(bool, Vec<[f64; DIMS]>)> = LinkedList::new();
    for poly in poly_list_src {
        poly_list_dst.push_back((is_cyclic, poly));
    }
    return poly_list_dst;
}

pub fn poly_list_without_cyclic(
    poly_list_src: LinkedList<(bool, Vec<[f64; DIMS]>)>,
) -> LinkedList<Vec<[f64; DIMS]>>
{
    let mut poly_list_dst: LinkedList<Vec<[f64; DIMS]>> = LinkedList::new();
    for (_is_cyclic, poly) in poly_list_src {
        poly_list_dst.push_back(poly);
    }
    return poly_list_dst;
}
*/

// Convert from float
pub fn poly_f64_from_i32(
    poly_int: &Vec<[i32; DIMS]>) -> Vec<[f64; DIMS]>
{
    let mut poly_float: Vec<[f64; DIMS]> = Vec::with_capacity(poly_int.len());
    for v_int in poly_int {
        poly_float.push({
            let mut v_as_float = [0.0; DIMS];
            for j in 0..DIMS {
                v_as_float[j] = v_int[j] as f64;
            }
            v_as_float
        });
    }
    return poly_float;
}
pub fn poly_list_f64_from_i32(
    poly_list_int: &LinkedList<(bool, Vec<[i32; DIMS]>)>,
) -> LinkedList<(bool, Vec<[f64; DIMS]>)>
{
    let mut poly_list_float: LinkedList<(bool, Vec<[f64; DIMS]>)> = LinkedList::new();
    for &(is_cyclic, ref poly_int) in poly_list_int {
        poly_list_float.push_back((is_cyclic, poly_f64_from_i32(&poly_int)));
    }
    return poly_list_float;
}

// Subdivide
pub fn poly_subdivide(
    is_cyclic: bool,
    poly_src: &Vec<[f64; DIMS]>,
) -> Vec<[f64; DIMS]>
{
    let mut poly_dst: Vec<[f64; DIMS]> = Vec::with_capacity(poly_src.len() * 2);
    let mut v_orig_prev = &poly_src[if is_cyclic { poly_src.len() - 1 } else { 0 }];
    if !is_cyclic {
        poly_dst.push(*v_orig_prev);
    }

    for v_orig_curr in &poly_src[(if is_cyclic { 0 } else { 1 })..] {
        // subdivided point
        poly_dst.push(mid_vnvn(v_orig_prev, v_orig_curr));
        // regular point
        poly_dst.push(*v_orig_curr);
        v_orig_prev = v_orig_curr;
    }
    return poly_dst;
}

pub fn poly_list_subdivide(
    poly_list_src: &LinkedList<(bool, Vec<[f64; DIMS]>)>,
) -> LinkedList<(bool, Vec<[f64; DIMS]>)>
{
    let mut poly_list_dst: LinkedList<(bool, Vec<[f64; DIMS]>)> = LinkedList::new();
    for &(is_cyclic, ref poly_src) in poly_list_src {
        poly_list_dst.push_back((is_cyclic, poly_subdivide(is_cyclic, poly_src)));
    }
    return poly_list_dst;
}

// Subdivide until segments are smaller then the limit
pub fn poly_subdivide_to_limit(
    is_cyclic: bool,
    poly_src: &Vec<[f64; DIMS]>,
    limit: f64,
) -> Vec<[f64; DIMS]>
{
    // target size isn't known. but will be at least as big as the source
    let mut poly_dst: Vec<[f64; DIMS]> = Vec::with_capacity(poly_src.len());

    let limit_sq = sq(limit);
    let mut v_orig_prev = &poly_src[if is_cyclic { poly_src.len() - 1 } else { 0 }];
    if !is_cyclic {
        poly_dst.push(*v_orig_prev);
    }

    for v_orig_curr in &poly_src[(if is_cyclic { 0 } else { 1 })..] {
        // subdivided point(s)
        let len_sq = len_squared_vnvn(v_orig_prev, v_orig_curr);
        if len_sq > limit_sq {
            let len = len_sq.sqrt();
            let sub = (len / limit).floor();
            let inc = 1.0 / sub;
            let mut step = inc;
            for _ in 0..((sub as usize) - 1) {
                poly_dst.push(interp_vnvn(v_orig_prev, v_orig_curr, step));
                debug_assert!(step > 0.0 && step < 1.0);
                step += inc;
            }
        }
        // regular point
        poly_dst.push(*v_orig_curr);
        v_orig_prev = v_orig_curr;
    }

    return poly_dst;
}

pub fn poly_list_subdivide_to_limit(
    poly_list_src: &LinkedList<(bool, Vec<[f64; DIMS]>)>, limit: f64,
) -> LinkedList<(bool, Vec<[f64; DIMS]>)>
{
    let mut poly_list_dst: LinkedList<(bool, Vec<[f64; DIMS]>)> = LinkedList::new();
    for &(is_cyclic, ref poly_src) in poly_list_src {
        poly_list_dst.push_back(
            (is_cyclic, poly_subdivide_to_limit(is_cyclic, poly_src, limit)));
    }
    return poly_list_dst;
}
