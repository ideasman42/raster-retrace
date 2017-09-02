
use intern::min_heap;

// 2d quadric
mod quadric {

    #[derive(Default, Clone)]
    pub struct Quadric {
        a2: f64, ab: f64, ac: f64,
        b2: f64, bc: f64,
        c2: f64,
    }

    fn to_tensor_matrix_inverse(
        q: &Quadric,
        epsilon: f64,
    ) -> Option<[f64; 3]> {
        let det: f64 =
            (q.a2 * q.b2) -
            (q.ab * q.ab);
        if det.abs() > epsilon {
            let invdet: f64 = 1.0 / det;
            /* 3 components of a 3x3 matrix,
             * we only use some of them, 4th would be identity (1.0) */
            return Some([
                q.b2 *  invdet,  /* [0][0] */
                q.ab * -invdet,  /* [0][1] */
                q.a2 *  invdet,  /* [1][1] */
            ]);
        } else {
            return None;
        }
    }

    // UNUSED
    /*
    pub fn to_position(
        q: &Quadric,
    ) -> [f64; 2] {
        return [
            q.ac,
            q.bc,
        ];
    }
    */

    pub fn add(
        q_a: &Quadric,
        q_b: &Quadric,
    ) -> Quadric {
        return Quadric {
            a2: q_a.a2 + q_b.a2,
            ab: q_a.ab + q_b.ab,
            ac: q_a.ac + q_b.ac,
            b2: q_a.b2 + q_b.b2,
            bc: q_a.bc + q_b.bc,
            c2: q_a.c2 + q_b.c2,
        };
    }

    pub fn iadd(
        q_a: &mut Quadric,
        q_b: &Quadric,
    ) {
        q_a.a2 += q_b.a2;
        q_a.ab += q_b.ab;
        q_a.ac += q_b.ac;
        q_a.b2 += q_b.b2;
        q_a.bc += q_b.bc;
        q_a.c2 += q_b.c2;
    }

    pub fn from_plane(
        v: &[f64; 3]
    ) -> Quadric {
        return Quadric {
            a2: v[0] * v[0],
            b2: v[1] * v[1],

            ab: v[0] * v[1],
            ac: v[0] * v[2],

            bc: v[1] * v[2],
            c2: v[2] * v[2],
        };
    }

    pub fn evaluate(
        q: &Quadric,
        v: &[f64; 2],
    ) -> f64 {
        return (q.a2 * v[0] * v[0]) + (q.ab * 2.0 * v[0] * v[1]) + (q.ac * 2.0 * v[0]) +
               (q.b2 * v[1] * v[1]) + (q.bc * 2.0 * v[1]) +
               (q.c2);
    }

    pub fn optimize(
        q: &Quadric,
        epsilon: f64,
    ) -> Option<[f64; 2]> {
        if let Some(m) = to_tensor_matrix_inverse(q, epsilon) {
            // 3x3 matrix multiply & negate
            // (ac, bc) == (x, y).
            return Some([
                -(m[0] * q.ac),
                -(m[1] * q.ac + m[2] * q.bc),
            ]);
        } else {
            return None;
        }
    }
}

#[inline(always)]
fn dot(a: &[f64; 2], b: &[f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
#[inline(always)]
fn len_sqr(a: &[f64; 2]) -> f64 {
    a[0] * a[0] + a[1] * a[1]
}
#[inline(always)]
fn len(a: &[f64; 2]) -> f64 {
    len_sqr(a).sqrt()
}
#[inline(always)]
fn normalized(a: &[f64; 2]) -> Option<[f64; 2]> {
    let l = len(a);
    if l != 0.0 {
        Some([a[0] / l, a[1] / l])
    } else {
        None
    }
}
#[inline(always)]
fn plane_from_point_normal(p: &[f64; 2], n: &[f64; 2]) -> [f64; 3] {
    [n[0], n[1], -dot(p, n)]
}

const INVALID: usize = ::std::usize::MAX;

struct Edge {
    v1: usize,
    v2: usize,

    index_prev: usize,
    index_next: usize,
}

#[derive(Copy, Clone)]
struct EdgeRemove {
    edge_index: usize,
    collapse_co: [f64; 2],
}

fn edge_heap_insert(
    poly_edit: &Vec<[f64; 2]>,
    quadrics: &Vec<quadric::Quadric>,
    heap: &mut min_heap::MinHeap<f64, EdgeRemove>,
    e: &Edge,
    e_handle: &mut min_heap::NodeHandle,
    i: usize,
    simplify_threshold_sq: f64,
) {
    use std::f64;

    let q1 = &quadrics[e.v1];
    let q2 = &quadrics[e.v2];
    let optimize_co = {
        if let Some(optimize_co) = quadric::optimize(&quadric::add(q1, q2), f64::EPSILON) {
            optimize_co
        } else {
            let v1 = &poly_edit[e.v1];
            let v2 = &poly_edit[e.v2];
            [
                (v1[0] + v2[0]) / 2.0,
                (v1[1] + v2[1]) / 2.0,
            ]
        }
    };

    let cost =
        (quadric::evaluate(q1, &optimize_co) +
         quadric::evaluate(q2, &optimize_co)).abs();

    *e_handle = {
        if cost < simplify_threshold_sq {
            heap.insert(
                cost,
                EdgeRemove {
                    edge_index: i,
                    collapse_co: optimize_co,
                }
            )
        } else {
            min_heap::NodeHandle::INVALID
        }
    };
}

fn edge_heap_update(
    poly_edit: &Vec<[f64; 2]>,
    quadrics: &Vec<quadric::Quadric>,
    heap: &mut min_heap::MinHeap<f64, EdgeRemove>,
    e: &Edge,
    e_handle: &mut min_heap::NodeHandle,
    i: usize,
    simplify_threshold_sq: f64,
) {
    if *e_handle != min_heap::NodeHandle::INVALID {
        heap.remove(*e_handle);
    }
    edge_heap_insert(
        poly_edit,
        quadrics,
        heap,
        e, e_handle, i,
        simplify_threshold_sq,
    );
}

const INVALID_CO: [f64; 2] = [::std::f64::MAX, ::std::f64::MAX];

fn edge_heap_collapse(
    poly_edit: &mut Vec<[f64; 2]>,
    quadrics: &mut Vec<quadric::Quadric>,
    heap: &mut min_heap::MinHeap<f64, EdgeRemove>,
    edges: &mut Vec<Edge>,
    edges_handle: &mut Vec<min_heap::NodeHandle>,
    i: usize,
    collapse_co: &[f64; 2],
    simplify_threshold_sq: f64,
) {
    let (i_prev, i_next) = {
        let e = &mut edges[i];
        let i_prev = e.index_prev;
        let i_next = e.index_next;

        // not needed, just expose invalid reuse
        e.index_prev = INVALID;
        e.index_next = INVALID;

        assert!(i_prev != INVALID);
        assert!(i_next != INVALID);

        assert!(poly_edit[e.v1] != INVALID_CO);
        assert!(poly_edit[e.v2] != INVALID_CO);

        e.v1 = INVALID;
        e.v2 = INVALID;

        (i_prev, i_next)
    };

    edges[i_prev].index_next = i_next;
    edges[i_next].index_prev = i_prev;

    // drop edges[i].v2 (could drop v1 too, doesn't matter which)

    let i_vert_keep = edges[i_prev].v2;
    let i_vert_drop = edges[i_next].v1;

    edges[i_next].v1 = i_vert_keep;

    poly_edit[i_vert_keep] = *collapse_co;
    poly_edit[i_vert_drop] = INVALID_CO;

    // let q = quadrics[i_vert_drop];
    quadrics[i_vert_keep] = quadric::add(&quadrics[i_vert_keep], &quadrics[i_vert_drop]);

    for i_other in &[
        i_prev, edges[i_prev].index_prev,
        i_next, edges[i_next].index_next,
    ] {
        // INVALID checks are needed for non-cyclic polygons.
        if *i_other != INVALID {
            let e = &mut edges[*i_other];
            if e.index_prev != INVALID && e.index_next != INVALID {
                edge_heap_update(
                    poly_edit,
                    quadrics,
                    heap,
                    e, &mut edges_handle[*i_other], *i_other,
                    simplify_threshold_sq,
                );
            }
        }
    }
}

pub fn poly_simplify(
    is_cyclic: bool,
    poly: &Vec<[f64; 2]>,
    simplify_threshold: f64,
) -> Vec<[f64; 2]> {
    // points we're allowed to adjust
    let mut poly_edit = poly.clone();
    let mut edges: Vec<Edge> = Vec::with_capacity(poly.len()  /* is_cyclic TODO */ );

    if is_cyclic {
        let mut j = poly.len() - 1;
        for i in 0..poly.len() {
            edges.push(Edge {
                v1: j,
                v2: i,
                index_prev: j,
                index_next: i.wrapping_add(1),
            });
            j = i;
        }
        edges.last_mut().unwrap().index_next = 0;
    } else {
        let mut j = 0;
        for i in 1..poly.len() {
            edges.push(Edge {
                v1: j,
                v2: i,
                index_prev: j.wrapping_sub(1),
                index_next: i,
            });
            j = i;
        }

        edges[0].index_next = INVALID;
        edges.last_mut().unwrap().index_next = INVALID;
    }

    let mut quadrics = vec![quadric::Quadric::default(); poly.len()];
    for e in &mut edges {
        // -y, x
        let p1 = &poly_edit[e.v1];
        let p2 = &poly_edit[e.v2];

        if let Some(n) = normalized(&[-(p1[1] - p2[1]), p1[0] - p2[0]]) {
            let p = plane_from_point_normal(p1, &n);
            let q = quadric::from_plane(&p);
            quadric::iadd(&mut quadrics[e.v1], &q);
            quadric::iadd(&mut quadrics[e.v2], &q);
        }
    }

    // Edges are setup, now collapse
    let simplify_threshold_sq = simplify_threshold * simplify_threshold;
    let mut heap = min_heap::MinHeap::<f64, EdgeRemove>::with_capacity(edges.len());
    let mut edges_handle = vec![min_heap::NodeHandle::INVALID; edges.len()];
    for i in {
        if is_cyclic {
            0..edges.len()
        } else {
            1..(edges.len() - 1)
        }
    } {
        edge_heap_insert(
            &mut poly_edit,
            &quadrics,
            &mut heap,
            &edges[i], &mut edges_handle[i], i,
            simplify_threshold_sq,
        );
    }

    let poly_minimum_len = if is_cyclic { 4 } else { 2 };
    let mut poly_remaining_len = poly.len();

    while let Some(r) = heap.pop_min() {
        // will never use again, set invalid for hygiene
        edges_handle[r.edge_index] = min_heap::NodeHandle::INVALID;
        if poly_remaining_len <= poly_minimum_len {
            break;
        }
        poly_remaining_len -= 1;

        edge_heap_collapse(
            &mut poly_edit,
            &mut quadrics,
            &mut heap,
            &mut edges,
            &mut edges_handle,
            r.edge_index,
            &r.collapse_co,
            simplify_threshold_sq,
        );
    }

    // reuse poly_edit and return as new polygon.
    let mut i_dst: usize = 0;
    for i_src in 0..poly_edit.len() {
        if poly_edit[i_src] != INVALID_CO {
            if i_dst != i_src {
                poly_edit[i_dst] = poly_edit[i_src];
            }
            i_dst += 1;
        }
    }
    poly_edit.truncate(i_dst);
    poly_edit.shrink_to_fit();

    return poly_edit;
}


use std::collections::LinkedList;

pub fn poly_list_simplify(
    poly_list_src: &LinkedList<(bool, Vec<[f64; 2]>)>,
    simplify_threshold: f64,
) -> LinkedList<(bool, Vec<[f64; 2]>)> {
    let mut poly_list_dst: LinkedList<(bool, Vec<[f64; 2]>)> = LinkedList::new();
    for &(is_cyclic, ref poly_src) in poly_list_src {
        poly_list_dst.push_back(
            (is_cyclic, poly_simplify(is_cyclic, poly_src, simplify_threshold)));
    }
    return poly_list_dst;
}

