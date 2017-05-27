
///
/// Module for writing curve data to files.
///

const DIMS: usize = ::intern::math_vector::DIMS;

pub mod svg {

    use super::{
        DIMS,
    };
    use std::collections::LinkedList;
    use std::io::prelude::Write;

    pub fn write_header(
        mut f: &::std::fs::File,
        size: &[usize; 2],
        scale: f64,
    ) -> Result<(), ::std::io::Error> {
        writeln!(f, "<?xml version='1.0' encoding='UTF-8'?>")?;
        writeln!(f, concat!(
            "<svg version='1.1' ",
            "width='{}' height='{}' ",
            "viewBox='0 0 {} {}' ",
            "xmlns='http://www.w3.org/2000/svg' ",
            "xmlns:xlink='http://www.w3.org/1999/xlink' ",
            ">"),
            scale * size[0] as f64,
            scale * size[1] as f64,
            scale * size[0] as f64,
            scale * size[1] as f64,
        )?;

        Ok(())
    }

    pub fn write_poly_list_filled(
        mut f: &::std::fs::File,
        _size: &[usize; 2],
        scale: f64,
        poly_list: &LinkedList<(bool, Vec<[f64; DIMS]>)>,
        pass_scale: f64,
    ) -> Result<(), ::std::io::Error> {
        use std::io::prelude::Write;

        f.write_fmt(format_args!(concat!("  ",
            "<g stroke='white' ",
            "stroke-opacity='0.5' ",
            "stroke-width='{:.2}' ",
            "fill='black' ",
            "fill-opacity='0.5' ",
            ">"),
            0.5 * pass_scale,
        ))?;

        f.write(b"    <path d='")?;
        for &(_is_cyclic, ref p) in poly_list {
            f.write(b"M ")?;
            for v in p {
                f.write_fmt(format_args!(
                    "{:.2},{:.2} ",
                    v[0] * scale,
                    v[1] * scale,
                ))?;
            }
            f.write(b" Z\n")?;
        }
        writeln!(f, "' />")?;

        writeln!(f, "  </g>")?;

        Ok(())
    }

    pub fn write_poly_list_centerline(
        mut f: &::std::fs::File,
        _size: &[usize; 2],
        scale: f64,
        poly_list: &LinkedList<(bool, Vec<[f64; DIMS]>)>,
        pass_scale: f64,
    ) -> Result<(), ::std::io::Error> {
        use std::io::prelude::Write;

        f.write_fmt(format_args!(concat!("  ",
            "<g stroke='grey' ",
            "stroke-opacity='0.75' ",
            "stroke-width='{:.2}' ",
            "fill='none' ",
            ">"),
            0.5 * pass_scale,
        ))?;

        f.write(b"    <path d='")?;
        for &(_is_cyclic, ref p) in poly_list {

            f.write(b"M ")?;
            for v in p {
                f.write_fmt(format_args!(
                    "{:.2},{:.2} ",
                    v[0] * scale,
                    v[1] * scale,
                ))?;
            }
        }
        writeln!(f, "' />")?;

        writeln!(f, "  </g>")?;

        Ok(())
    }

    pub fn write_curve_list_with_tangent_info(
        mut f: &::std::fs::File,
        scale: f64,
        poly_list: &LinkedList<(bool, Vec<[[f64; DIMS]; 3]>)>,
        pass_scale: f64,
    ) -> Result<(), ::std::io::Error> {
        // handle segments
        {
            f.write_fmt(format_args!(concat!("  ",
                "<g stroke='black' ",
                "stroke-opacity='0.5' ",
                "stroke-width='{:.2}' ",
                ">"),
                2.0 * pass_scale,
            ))?;
            for &(_is_cyclic, ref p) in poly_list {
                for v in p {
                    f.write_fmt(format_args!(
                        "<line x1='{:.2}' y1='{:.2}' x2='{:.2}' y2='{:.2}' />",
                        v[0][0] * scale, v[0][1] * scale,
                        v[1][0] * scale, v[1][1] * scale,
                    ))?;
                    f.write_fmt(format_args!(
                        "<line x1='{:.2}' y1='{:.2}' x2='{:.2}' y2='{:.2}' />",
                        v[1][0] * scale, v[1][1] * scale,
                        v[2][0] * scale, v[2][1] * scale,
                    ))?;
                }
            }
            writeln!(f, "  </g>")?;
        }

        // circle's
        {
            f.write_fmt(format_args!(concat!("  ",
                "<g stroke='white' ",
                "stroke-opacity='1.0' ",
                "stroke-width='{:.2}' ",
                "fill='black' ",
                "fill-opacity='0.5' ",
                ">"),
                1.0 * pass_scale,
            ))?;

            for &(_is_cyclic, ref p) in poly_list {
                for v in p {
                    for h in v {
                        f.write_fmt(format_args!(
                            "<circle cx='{:.2}' cy='{:.2}' r='{:.2}'/>",
                            h[0] * scale,
                            h[1] * scale,
                            2.0 * pass_scale,
                        ))?;
                    }

                    f.write_fmt(format_args!(
                        "<line x1='{:.2}' y1='{:.2}' x2='{:.2}' y2='{:.2}' />",
                        v[0][0] * scale, v[0][1] * scale,
                        v[1][0] * scale, v[1][1] * scale,
                    ))?;
                    f.write_fmt(format_args!(
                        "<line x1='{:.2}' y1='{:.2}' x2='{:.2}' y2='{:.2}' />",
                        v[1][0] * scale, v[1][1] * scale,
                        v[2][0] * scale, v[2][1] * scale,
                    ))?;
                }
            }
            writeln!(f, "  </g>")?;
        }

        Ok(())
    }

    pub fn write_curve_list_filled(
        mut f: &::std::fs::File,
        _size: &[usize; 2],
        scale: f64,
        poly_list: &LinkedList<(bool, Vec<[[f64; DIMS]; 3]>)>,
    ) -> Result<(), ::std::io::Error> {
        use std::io::prelude::Write;

        writeln!(f, concat!("  ",
            "<g stroke='black' ",
            "stroke-opacity='0.0' ",
            "stroke-width='0' ",
            "fill='black' ",
            "fill-opacity='1' ",
            ">",
        ))?;

        f.write(b"    <path d='")?;
        for &(_is_cyclic, ref p) in poly_list {
            let mut v_prev = p.last().unwrap();
            let mut is_first = true;
            for v_curr in p {

                use intern::math_vector::{
                    is_finite_vn
                };
                debug_assert!(is_finite_vn(&v_curr[0]));
                debug_assert!(is_finite_vn(&v_curr[1]));
                debug_assert!(is_finite_vn(&v_curr[2]));

                let k0 = v_prev[1];
                let h0 = v_prev[2];

                let h1 = v_curr[0];
                let k1 = v_curr[1];

                // Could optimize this, but keep now for simplicity
                if is_first {
                    f.write_fmt(format_args!(
                        "M {:.2},{:.2} ",
                        k0[0] * scale,
                        k0[1] * scale,
                    ))?;
                }
                f.write_fmt(format_args!(
                    "C {:.2},{:.2} {:.2},{:.2} {:.2},{:.2} ",
                    h0[0] * scale, h0[1] * scale,
                    h1[0] * scale, h1[1] * scale,
                    k1[0] * scale, k1[1] * scale,
                ))?;
                v_prev = v_curr;
                is_first = false;
            }

            f.write(b" Z\n")?;

        }
        writeln!(f, "' />")?;

        writeln!(f, "  </g>")?;

        Ok(())
    }

    pub fn write_curve_list_centerline(
        mut f: &::std::fs::File,
        _size: &[usize; 2],
        scale: f64,
        poly_list: &LinkedList<(bool, Vec<[[f64; DIMS]; 3]>)>,
    ) -> Result<(), ::std::io::Error> {
        use std::io::prelude::Write;

        writeln!(f, concat!("  ",
            "<g stroke='black' ",
            "stroke-opacity='1.0' ",
            "stroke-width='1' ",
            "fill='none' ",
            ">",
        ))?;

        for &(is_cyclic, ref p) in poly_list {
            if is_cyclic {
                f.write(b"    <path d='")?;
                let mut v_prev = p.last().unwrap();
                let mut is_first = true;
                for v_curr in p {

                    use intern::math_vector::{
                        is_finite_vn,
                    };
                    debug_assert!(is_finite_vn(&v_curr[0]));
                    debug_assert!(is_finite_vn(&v_curr[1]));
                    debug_assert!(is_finite_vn(&v_curr[2]));

                    let k0 = v_prev[1];
                    let h0 = v_prev[2];

                    let h1 = v_curr[0];
                    let k1 = v_curr[1];

                    // Could optimize this, but keep now for simplicity
                    if is_first {
                        f.write_fmt(format_args!(
                            "M {:.2},{:.2} ",
                            k0[0] * scale,
                            k0[1] * scale,
                        ))?;
                    }
                    f.write_fmt(format_args!(
                        "C {:.2},{:.2} {:.2},{:.2} {:.2},{:.2} ",
                        h0[0] * scale, h0[1] * scale,
                        h1[0] * scale, h1[1] * scale,
                        k1[0] * scale, k1[1] * scale,
                    ))?;
                    v_prev = v_curr;
                    is_first = false;
                }
                f.write(b" Z\n")?;
                writeln!(f, "' />")?;
            } else {
                f.write(b"    <path d='")?;

                let mut v_prev = &p[0];
                let mut is_first = true;
                for v_curr in &p[1..p.len()] {

                    use intern::math_vector::{
                        is_finite_vn,
                    };
                    debug_assert!(is_finite_vn(&v_curr[0]));
                    debug_assert!(is_finite_vn(&v_curr[1]));
                    debug_assert!(is_finite_vn(&v_curr[2]));

                    let k0 = v_prev[1];
                    let h0 = v_prev[2];

                    let h1 = v_curr[0];
                    let k1 = v_curr[1];

                    // Could optimize this, but keep now for simplicity
                    if is_first {
                        f.write_fmt(format_args!(
                            "M {:.2},{:.2} ",
                            k0[0] * scale,
                            k0[1] * scale,
                        ))?;
                    }
                    f.write_fmt(format_args!(
                        "C {:.2},{:.2} {:.2},{:.2} {:.2},{:.2} ",
                        h0[0] * scale, h0[1] * scale,
                        h1[0] * scale, h1[1] * scale,
                        k1[0] * scale, k1[1] * scale,
                    ))?;
                    v_prev = v_curr;
                    is_first = false;
                }

                writeln!(f, "' />")?;
            }
        }

        writeln!(f, "  </g>")?;

        Ok(())
    }

    pub fn write_footer(
        mut f: &::std::fs::File,
    ) -> Result<(), ::std::io::Error> {
        use std::io::prelude::Write;
        writeln!(f, "</svg>")?;
        Ok(())
    }

/*
    pub fn write_full(
        f: &::std::fs::File,
        size: &[usize; 2],
        scale: f64,
        poly_list: &LinkedList<(bool, Vec<[f64; DIMS]>)>,
    ) -> Result<(), ::std::io::Error> {
        write_header(f, size, scale)?;
        write_poly_list_filled(f, size, scale, poly_list)?;
        write_footer(f)?;
        Ok(())
    }
*/
}

