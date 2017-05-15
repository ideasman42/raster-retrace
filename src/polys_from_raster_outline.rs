///
/// Takes an images and returns multiple curves
/// representing the outline of pixel regions.
///

const DIMS: usize = ::intern::math_vector::DIMS;

macro_rules! ensure_const_expr {
    ($value:expr, $t:ty) => {
        {
            const _IGNORE: $t = $value;
        }
    }
}

macro_rules! elem {
    ($val:expr, $($var:expr), *) => {
        $($val == $var) || *
    }
}

use std::collections::LinkedList;

#[derive(Copy, Clone)]
pub enum TurnPolicy {
    Black,
    White,
    Majority,
    Minority,
}

// TODO, split into own file?
//
///
/// Perform the image to bitmap outline generation.
///
/// * `use_simplify` - don't write intermediate points (one per pixel) between corners.
pub fn extract_outline(
    image: &[bool],
    size: &[usize; 2],
    turn_policy: TurnPolicy,
    use_simplify: bool,
) -> LinkedList<(bool, Vec<[i32; DIMS]>)> {
    mod dir {
        pub const L: u8 = (1 << 0);
        pub const R: u8 = (1 << 1);
        pub const D: u8 = (1 << 2);
        pub const U: u8 = (1 << 3);
    }

    let psize: [usize; 2] = [size[0] + 1, size[1] + 1];
    let mut pimage: Vec<u8> = vec![0; psize[0] * psize[1]];

    // assumed in-range
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

    // note, the borders could have special handling for more efficient checks
    let mut steps_total = 0;
    for y in 0..size[1] {
        for x in 0..size[0] {
            let index = xy!(x, y, size[0]);
            if image[index] {
                if !xy_is_filled_l!(x, y) {
                    pimage[xy!(x + 0, y + 0, psize[0])] |= dir::U;
                    steps_total += 1;
                }
                if !xy_is_filled_r!(x, y) {
                    pimage[xy!(x + 1, y + 1, psize[0])] |= dir::D;
                    steps_total += 1;
                }
                if !xy_is_filled_d!(x, y) {
                    pimage[xy!(x + 1, y + 0, psize[0])] |= dir::L;
                    steps_total += 1;
                }
                if !xy_is_filled_u!(x, y) {
                    pimage[xy!(x + 0, y + 1, psize[0])] |= dir::R;
                    steps_total += 1;
                }
            }
        }
    }

    let mut poly_list = LinkedList::new();
    {
        fn poly_from_direction_mask(
            pimage: &mut Vec<u8>,
            x_init: i32,
            y_init: i32,
            x_span: i32,
            // only needed for checking majority turning
            image_data: &(&[bool], [i32; 2]),
            turn_policy: TurnPolicy,
            use_simplify: bool,
            direction_init_prev: u8,
        ) -> (Vec<[i32; DIMS]>, usize) {
            let mut poly: Vec<[i32; DIMS]> = vec![];
            let mut x = x_init;
            let mut y = y_init;
            let mut d_prev: u8 = direction_init_prev;
            let mut handled: usize = 0;
            loop {
                if use_simplify &&
                   (poly.len() > 1) && {
                        let xy_a = &poly[poly.len() - 2];
                        let xy_b = &poly[poly.len() - 1];
                        {
                            ((x == xy_a[0] && x == xy_b[0]) ||
                             (y == xy_a[1] && y == xy_b[1]))
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

                if handled != 0 {
                    if x == x_init && y == y_init {
                        poly.pop(); // weak
                        break;
                    }
                }

                let index: usize = xy!(x, y, x_span) as usize;
                let d = pimage[index];

                macro_rules! step_move {
                    ($dir:expr) => {
                        match $dir {
                            dir::L => { x -= 1; }
                            dir::R => { x += 1; }
                            dir::D => { y -= 1; }
                            dir::U => { y += 1; }
                            _ => { unreachable!(); }

                        }
                    }
                }

                // step along the first match
                macro_rules! step_first_match {
                    // ensure we're constant so the following expression is a constant
                    ($a:expr, $b:expr, $c:expr) => {
                        {
                            ensure_const_expr!($a, u8);
                            ensure_const_expr!($b, u8);
                            ensure_const_expr!($c, u8);

                            if      (d & $a) != 0 { step_move!($a); $a }
                            else if (d & $b) != 0 { step_move!($b); $b }
                            else if (d & $c) != 0 { step_move!($c); $c }
                            else { unreachable!(); }
                        }
                    }
                }

                fn is_majority(
                    x: i32,
                    y: i32,
                    image_data: &(&[bool], [i32; 2]),
                ) -> bool {

                    macro_rules! xy_or {
                        ($x:expr, $y:expr, $default:expr) => {
                            if ($x >= 0 && $x < image_data.1[0]) &&
                               ($y >= 0 && $y < image_data.1[1])
                            {
                                image_data.0[xy!($x, $y, image_data.1[0]) as usize]
                            } else {
                                $default
                            }
                        }
                    }

                    for i in 2..5 {
                        let mut ct: i32 = 0;
                        for a in (-i + 1)..i {
                            ct += if xy_or!(x + a,     y + i - 1, false) { 1 } else { -1 };
                            ct += if xy_or!(x + i - 1, y + a - 1, false) { 1 } else { -1 };
                            ct += if xy_or!(x + a - 1, y - i,     false) { 1 } else { -1 };
                            ct += if xy_or!(x - i,     y + a,     false) { 1 } else { -1 };
                        }
                        if ct > 0 {
                            return true;
                        } else if ct < 0 {
                            return false;
                        }
                    }
                    return false;
                }

                // From the previous direction,
                // take the nearest next step in a counter-clockwise order.

                let d_next: u8 = {
                    if elem!(d, dir::L, dir::R, dir::D, dir::U) {
                        // non-ambiguous case
                        step_move!(d);
                        d
                    } else {
                        // ambiguous case
                        let turn_ccw: bool = {
                            match turn_policy {
                                TurnPolicy::Black => { true },
                                TurnPolicy::White => { false },
                                TurnPolicy::Majority => {  is_majority(x, y, image_data) },
                                TurnPolicy::Minority => { !is_majority(x, y, image_data) },
                            }
                        };

                        if turn_ccw == false {
                            match d_prev {
                                dir::L => { step_first_match!(dir::D, dir::L, dir::U) },
                                dir::U => { step_first_match!(dir::L, dir::U, dir::R) },
                                dir::R => { step_first_match!(dir::U, dir::R, dir::D) },
                                dir::D => { step_first_match!(dir::R, dir::D, dir::L) },
                                _ => { unreachable!(); }
                            }
                        } else {
                            match d_prev {
                                dir::L => { step_first_match!(dir::U, dir::L, dir::D) },
                                dir::U => { step_first_match!(dir::R, dir::U, dir::L) },
                                dir::R => { step_first_match!(dir::D, dir::R, dir::U) },
                                dir::D => { step_first_match!(dir::L, dir::D, dir::R) },
                                _ => { unreachable!(); }
                            }
                        }
                    }
                };

                // never walk this direction again
                pimage[index] &= !d_next;
                d_prev = d_next;

                handled += 1;
            }

            return (poly, handled);
        }

        let mut steps_handled: usize = 0;

        let image_data = (image, [size[0] as i32, size[1] as i32]);

        'outer:
        for y in 0..psize[1] {
            for x in 0..psize[0] {
                let index = xy!(x, y, psize[0]);
                let d = pimage[index];
                // always start searching for up, since we do clockwise search
                if (d & dir::U) != 0 {
                    let (poly, handled) = poly_from_direction_mask(
                        &mut pimage,
                        x as i32,
                        y as i32,
                        psize[0] as i32,
                        &image_data,
                        turn_policy,
                        use_simplify, dir::L);
                    poly_list.push_back((true, poly));
                    steps_handled += handled;

                    if steps_total == steps_handled {
                        break 'outer;
                    }
                }
            }
        }
    }
    return poly_list;
}
