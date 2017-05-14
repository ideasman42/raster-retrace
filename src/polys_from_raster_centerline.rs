///
/// Takes an images and returns multiple curves
/// representing the center line for pixel regions.
///
/// Note, the image needs to be pre-processed
/// to ensure lines are only ever 1 pixel width.
///

use std::collections::LinkedList;

const DIMS: usize = ::intern::math_vector::DIMS;

pub fn extract_centerline(
    image: &[bool],
    size: &[usize; 2],
    use_simplify: bool,
) -> LinkedList<(bool, Vec<[i32; DIMS]>)>
{

    mod dir {
        pub const L: u8 = (1 << 0);
        pub const R: u8 = (1 << 1);
        pub const D: u8 = (1 << 2);
        pub const U: u8 = (1 << 3);

        pub const LD: u8 = (1 << 4);
        pub const LU: u8 = (1 << 5);
        pub const RD: u8 = (1 << 6);
        pub const RU: u8 = (1 << 7);
    }

    macro_rules! xy {
        ($x:expr, $y:expr, $x_span:expr) => {
            $x + ($y * $x_span)
        }
    }

    macro_rules! xy_is_filled_l {
        ($x:expr, $y:expr) => {
            if $x != 0 {
                image[xy!($x - 1, $y, size[0])]
            } else {
                false
            }
        }
    }
    macro_rules! xy_is_filled_r {
        ($x:expr, $y:expr) => {
            if $x + 1 != size[0] {
                image[xy!($x + 1, $y, size[0])]
            } else {
                false
            }
        }
    }
    macro_rules! xy_is_filled_d {
        ($x:expr, $y:expr) => {
            if $y != 0 {
                image[xy!($x, $y - 1, size[0])]
            } else {
                false
            }
        }
    }
    macro_rules! xy_is_filled_u {
        ($x:expr, $y:expr) => {
            if $y + 1 != size[1] {
                image[xy!($x, $y + 1, size[0])]
            } else {
                false
            }
        }
    }

    // diagonals
    macro_rules! xy_is_filled_ld {
        ($x:expr, $y:expr) => {
            if $x != 0 && $y != 0 {
                image[xy!($x - 1, $y - 1, size[0])]
            } else {
                false
            }
        }
    }
    macro_rules! xy_is_filled_lu {
        ($x:expr, $y:expr) => {
            if $x != 0 && $y + 1 != size[1] {
                image[xy!($x - 1, $y + 1, size[0])]
            } else {
                false
            }
        }
    }
    macro_rules! xy_is_filled_rd {
        ($x:expr, $y:expr) => {
            if $x + 1 != size[0] && $y != 0 {
                image[xy!($x + 1, $y - 1, size[0])]
            } else {
                false
            }
        }
    }
    macro_rules! xy_is_filled_ru {
        ($x:expr, $y:expr) => {
            if $x + 1 != size[0] && $y + 1 != size[1] {
                image[xy!($x + 1, $y + 1, size[0])]
            } else {
                false
            }
        }
    }

    let mut pimage: Vec<u8> = vec![0; size[0] * size[1]];

    // note, the borders could have special handling for more efficient checks
    for y in 0..size[1] {
        for x in 0..size[0] {
            let index = xy!(x, y, size[0]);
            if image[index] {
                let mut count = 0;
                let mut pf: u8 = 0;

                if xy_is_filled_l!(x, y) {
                    pf |= dir::L;
                    count += 1;
                }
                if xy_is_filled_r!(x, y) {
                    pf |= dir::R;
                    count += 1;
                }
                if xy_is_filled_d!(x, y) {
                    pf |= dir::D;
                    count += 1;
                }
                if xy_is_filled_u!(x, y) {
                    pf |= dir::U;
                    count += 1;
                }

                // connect diagonals when we _only_ have a diagonal connections.
                if (pf & (dir::L | dir::D)) == 0 && xy_is_filled_ld!(x, y) {
                    pf |= dir::LD;
                    count += 1;
                }
                if (pf & (dir::L | dir::U)) == 0 && xy_is_filled_lu!(x, y) {
                    pf |= dir::LU;
                    count += 1;
                }
                if (pf & (dir::R | dir::D)) == 0 && xy_is_filled_rd!(x, y) {
                    pf |= dir::RD;
                    count += 1;
                }
                if (pf & (dir::R | dir::U)) == 0 && xy_is_filled_ru!(x, y) {
                    pf |= dir::RU;
                    count += 1;
                }

                // only walk _to_ 3+ connections, never from.
                if count > 0 && count < 3 {
                    pimage[index] = pf;
                }
            }
        }
    }

    let mut poly_list: LinkedList<(bool, Vec<[i32; DIMS]>)> = LinkedList::new();
    {
        fn poly_from_direction_mask_half(
            pimage: &mut Vec<u8>,
            x_init: i32,
            y_init: i32,
            x_span: usize,
            use_simplify: bool,
            // direction_init: u8,
        ) -> (bool, Vec<[i32; DIMS]>)
        {
            let mut poly: Vec<[i32; DIMS]> = vec![];
            let mut is_cyclic = false;

            let mut x = x_init;
            let mut y = y_init;

            let mut index = xy!(x_init as usize, y_init as usize, x_span);
            loop {
                debug_assert!(index == xy!(x as usize, y as usize, x_span));

                if use_simplify &&
                   (poly.len() > 1) && {
                        let xy_a = &poly[poly.len() - 2];
                        let xy_b = &poly[poly.len() - 1];
                        {
                            (
                                // axis aligned
                                (x == xy_a[0] && x == xy_b[0]) ||
                                (y == xy_a[1] && y == xy_b[1]) ||
                                // diagonal
                                {
                                    let x_a_delta = xy_a[0] - xy_b[0];
                                    let y_a_delta = xy_a[1] - xy_b[1];
                                    let x_b_delta = xy_b[0] - x;
                                    let y_b_delta = xy_b[1] - y;

                                    (x_a_delta != 0 && y_a_delta != 0 &&
                                     // no need to check 'b', signum accounts for that case.

                                     x_a_delta.abs() == y_a_delta.abs() &&
                                     x_b_delta.abs() == y_b_delta.abs() &&


                                     x_a_delta.signum() == x_b_delta.signum() &&
                                     y_a_delta.signum() == y_b_delta.signum())
                                }
                             )
                        }
                   }
                {
                    let xy = poly.last_mut().unwrap();
                    xy[0] = x;
                    xy[1] = y;
                } else {
                    poly.push({
                        let mut xy: [i32; DIMS] = [0; DIMS];
                        xy[0] = x;
                        xy[1] = y;
                        xy
                    });
                }

                let f = pimage[index];
                pimage[index] = 0;

                if (f & dir::L) != 0 { // axis aligned
                    x -= 1;
                    index = index - 1;
                    pimage[index] &= !dir::R;
                } else if (f & dir::R) != 0 {
                    x += 1;
                    index = index + 1;
                    pimage[index] &= !dir::L;
                } else if (f & dir::D) != 0 {
                    y -= 1;
                    index = index - x_span;
                    pimage[index] &= !dir::U;
                } else if (f & dir::U) != 0 {
                    y += 1;
                    index = index + x_span;
                    pimage[index] &= !dir::D;
                } else if (f & dir::LD) != 0 { // diagonals
                    x -= 1;
                    y -= 1;
                    index = (index - 1) - x_span;
                    pimage[index] &= !dir::RU;
                } else if (f & dir::LU) != 0 {
                    x -= 1;
                    y += 1;
                    index = (index - 1) + x_span;
                    pimage[index] &= !dir::RD;
                } else if (f & dir::RD) != 0 {
                    x += 1;
                    y -= 1;
                    index = (index + 1) - x_span;
                    pimage[index] &= !dir::LU;
                } else if (f & dir::RU) != 0 {
                    x += 1;
                    y += 1;
                    index = (index + 1) + x_span;
                    pimage[index] &= !dir::LD;
                } else {
                    break;
                }

                if x == x_init &&
                   y == y_init
                {
                    is_cyclic = true;
                    break;
                }
            }

            return (is_cyclic, poly);
        }


        fn poly_from_direction_mask(
            pimage: &mut Vec<u8>,
            x_init: i32,
            y_init: i32,
            x_span: usize,
            use_simplify: bool,
        ) -> (bool, Vec<[i32; DIMS]>)
        {
            let index = xy!(x_init as usize, y_init as usize, x_span);

            let mut f = pimage[index];

            let (is_cyclic, mut poly) = poly_from_direction_mask_half(
                pimage, x_init, y_init, x_span, use_simplify);
            if is_cyclic == false {
                // remove the first direction, walk the next
                for i in 0..8 {
                    if (f & (1 << i)) != 0 {
                        f &= !(1 << i);
                        break;
                    }
                }
                pimage[index] = f;
                let (_, poly_half) = poly_from_direction_mask_half(
                    pimage, x_init, y_init, x_span, use_simplify);
                // could be more efficient
                poly.reverse();
                // avoid doubling up
                poly.pop();
                poly.extend(poly_half);
            }

            return (is_cyclic, poly);
        }

        for y in 0..size[1] {
            for x in 0..size[0] {
                let index = xy!(x, y, size[0]);
                if pimage[index] != 0 {
                    // walk in 2 directions!
                    let p = poly_from_direction_mask(
                        &mut pimage, x as i32, y as i32, size[0], use_simplify);
                    poly_list.push_back(p);
                }
            }
        }


        return poly_list;
    }
}

