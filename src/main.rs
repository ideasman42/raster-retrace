///
/// Raster Re-Trace, Main function.
///
/// Handles command line arguments,
/// image loading and calling tracing functionality.
///


mod intern;

mod polys_utils;
mod polys_from_raster_outline;
mod polys_from_raster_centerline;

mod polys_simplify_collapse;

mod image_skeletonize;

use std::collections::LinkedList;

// IO
mod curve_write;

use ::intern::{
    curve_fit_nd,
};


const PRINT_STATISTICS: bool = true;

/// Debug passes:
/// useful when investigating changes to internal behavior.
mod debug_pass {
    const DIMS: usize = ::intern::math_vector::DIMS;
    use std::collections::LinkedList;

    pub mod kind {
        /// polygon as extracted from pixels
        pub const PIXEL: u32 = 1 << 0;
        /// polygon before fitting calculation
        pub const PRE_FIT: u32 = 1 << 1;
        /// bezier handles
        pub const TANGENT: u32 = 1 << 2;
    }
    // passes that write out debug info
    pub struct Item {
        pub poly_list: LinkedList<(bool, Vec<[f64; DIMS]>)>,
    }

    pub fn add_pass(
        pass_items: &mut LinkedList<Item>,
        poly_list: &LinkedList<(bool, Vec<[f64; DIMS]>)>,
    )
    {
        pass_items.push_back(
            Item {
                poly_list: poly_list.clone(),
            }
        );
    }
}

pub fn trace_image(
    output_filepath: &String,
    output_scale: f64,
    image: &[bool],
    size: &[usize; 2],
    error_threshold: f64,
    simplify_threshold: f64,
    corner_angle: f64,
    use_optimize_exhaustive: bool,
    length_threshold: f64,
    mode: curve_fit_nd::TraceMode,
    // only for outline
    turn_policy: polys_from_raster_outline::TurnPolicy,
    debug_passes: u32,
    debug_pass_scale: f64,
) -> Result<(), ::std::io::Error>
{
    debug_assert!(size[0] * size[1] == image.len());

    // TODO, we could split these operations per-polygon
    // so they can be easily threaded.

    let mut pass_items: LinkedList<debug_pass::Item> = LinkedList::new();

    let poly_list_to_fit = {
        let poly_list_int = match mode {
            intern::curve_fit_nd::TraceMode::Outline => {
                polys_from_raster_outline::extract_outline(
                    image, &size,
                    turn_policy,
                    true)
            }
            curve_fit_nd::TraceMode::Centerline => {
                use polys_from_raster_centerline;

                polys_from_raster_centerline::extract_centerline(
                    image, &size, true)
            }
        };

        let poly_list_dst =
            polys_utils::poly_list_f64_from_i32(&poly_list_int);

        if (debug_passes & debug_pass::kind::PIXEL) != 0 {
            debug_pass::add_pass(&mut pass_items, &poly_list_dst);
        }

        // Ensure we always have at least one knot between 'corners'
        // this means theres always a middle tangent, giving us more possible
        // tangents when fitting the curve.
        let poly_list_dst =
            polys_utils::poly_list_subdivide(&poly_list_dst);

        let poly_list_dst =
            polys_simplify_collapse::poly_list_simplify(&poly_list_dst, simplify_threshold);

        if (debug_passes & debug_pass::kind::PRE_FIT) != 0 {
            debug_pass::add_pass(&mut pass_items, &poly_list_dst);
        }

        let poly_list_dst =
            polys_utils::poly_list_subdivide(&poly_list_dst);


        // While a little excessive, setting the `length_threshold` around 1.0
        // helps by ensure the density of the polygon is even
        // (without this diagonals will have many more points).
        let poly_list_dst = polys_utils::poly_list_subdivide_to_limit(
            &poly_list_dst, length_threshold);

        poly_list_dst
    };

    // if (debug_passes & debug_pass::kind::PRE_FIT) != 0 {
        // debug_pass::add_pass(&mut pass_items, &poly_list_to_fit);
    // }

    let curve_list =
        curve_fit_nd::fit_poly_list(
            poly_list_to_fit,
            error_threshold,
            corner_angle,
            use_optimize_exhaustive,
        );

    if PRINT_STATISTICS {
        let mut total_points = 0;
        for poly in &curve_list {
            total_points += poly.1.len();
        }
        println!("Total points: {}\n", total_points);
    }

    let f = ::std::fs::File::create(output_filepath).expect("Create output file");
    {
        curve_write::svg::write_header(&f, &size, output_scale)?;

        match mode {
            curve_fit_nd::TraceMode::Outline => {
                curve_write::svg::write_curve_list_filled(
                    &f, &size, output_scale, &curve_list)?;
            },
            curve_fit_nd::TraceMode::Centerline => {
                curve_write::svg::write_curve_list_centerline(
                    &f, &size, output_scale, &curve_list)?;
            }
        };

        // debug info, for developing mostly
        {
            for item in pass_items {
                match mode {
                    curve_fit_nd::TraceMode::Outline => {
                        curve_write::svg::write_poly_list_filled(
                            &f, &size, output_scale, &item.poly_list, debug_pass_scale)?;
                    },
                    curve_fit_nd::TraceMode::Centerline => {
                        curve_write::svg::write_poly_list_centerline(
                            &f, &size, output_scale, &item.poly_list, debug_pass_scale)?;
                    }
                };

            }
            if (debug_passes & debug_pass::kind::TANGENT) != 0 {
                curve_write::svg::write_curve_list_with_tangent_info(
                    &f, output_scale, &curve_list, debug_pass_scale)?;
            }
        }

        curve_write::svg::write_footer(&f)?;
    }

    Ok(())
}

#[derive(Clone)]
pub struct TraceParams {
    pub error_threshold: f64,
    pub simplify_threshold: f64,
    pub corner_threshold: f64,
    pub use_optimize_exhaustive: bool,
    pub input_filepath: String,
    pub output_filepath: String,
    pub output_scale: f64,
    pub mode: curve_fit_nd::TraceMode,
    pub turn_policy: polys_from_raster_outline::TurnPolicy,

    pub debug_passes: u32,
    pub debug_pass_scale: f64,

    show_help: bool,
}

impl Default for TraceParams {
    fn default(
    ) -> TraceParams
    {
        TraceParams {
            error_threshold: 1.0,
            simplify_threshold: 2.5,
            corner_threshold: 30.0_f64.to_radians(),
            use_optimize_exhaustive: false,
            input_filepath: String::new(),
            output_filepath: String::new(),
            output_scale: 1.0,
            mode: curve_fit_nd::TraceMode::Outline,
            turn_policy: polys_from_raster_outline::TurnPolicy::Majority,
            debug_passes: 0,
            debug_pass_scale: 1.0,

            show_help: false,
        }
    }
}

fn main()
{
    use intern::argparse;
    let mut trace_params = TraceParams::default();

    // -----------------------------------------------------------------------
    // Parse Args
    {
        use std::str::FromStr;

        let mut parser = argparse::new(
            &mut trace_params,
            "Bitmap image tracing utility",
            );

        // File Options
        {
            let parser_group = Some(parser.add_argument_group(
                "File Options",
                ""
            ));
            parser.add_argument(
                "-i", "--input",
                "The file path to use for input",
                "FILEPATH",
                Box::new(|dest_data, my_args| {
                    dest_data.input_filepath = my_args[0].clone();
                    return Ok(1);
                }),
                1, argparse::ARGDEF_DEFAULT | argparse::ARGDEF_REQUIRED,
                parser_group,
            );
            parser.add_argument(
                "-o", "--output",
                "The file path to use for writing",
                "FILEPATH",
                Box::new(|dest_data, my_args| {
                    dest_data.output_filepath = my_args[0].clone();
                    return Ok(1);
                }),
                1, argparse::ARGDEF_DEFAULT | argparse::ARGDEF_REQUIRED,
                parser_group,
            );
        }

        // Tracing Methods
        {
            let parser_group = Some(parser.add_argument_group(
                "Tracing Behavior",
                ""
            ));
            parser.add_argument(
                "-m", "--mode",
                concat!("The method used for tracing the image in [OUTLINE, CENTER], ",
                        "(defaults to OUTLINE)."),
                "MODE",
                Box::new(|dest_data, my_args| {
                    match my_args[0].as_ref() {
                        "OUTLINE" => {
                            dest_data.mode = curve_fit_nd::TraceMode::Outline;
                        },
                        "CENTER" => {
                            dest_data.mode = curve_fit_nd::TraceMode::Centerline;
                        },
                        _ => {
                            return Err(format!(
                                "Expected [OUTLINE, CENTER], not '{}'",
                                my_args[0],
                            ));
                        }
                    }
                    return Ok(1);
                }),
                1, argparse::ARGDEF_DEFAULT,
                parser_group,
            );
            parser.add_argument(
                "-z", "--turnpolicy",
                concat!("Method for extracting outlines [BLACK, WHITE, MAJORITY, MINORITY], ",
                        "(defaults to MAJORITY)."),
                "POLICY",
                Box::new(|dest_data, my_args| {
                    match my_args[0].as_ref() {
                        "BLACK" => {
                            dest_data.turn_policy =
                                polys_from_raster_outline::TurnPolicy::Black;
                        }
                        "WHITE" => {
                            dest_data.turn_policy =
                                polys_from_raster_outline::TurnPolicy::White;
                        }
                        "MAJORITY" => {
                            dest_data.turn_policy =
                                polys_from_raster_outline::TurnPolicy::Majority;
                        }
                        "MINORITY" => {
                            dest_data.turn_policy =
                                polys_from_raster_outline::TurnPolicy::Minority;
                        }
                        _ => {
                            return Err(format!(
                                "Expected [BLACK, WHITE, MAJORITY, MINORITY], not '{}'",
                                my_args[0],
                            ));
                        }
                    }
                    return Ok(1);
                }),
                1, argparse::ARGDEF_DEFAULT,
                parser_group,
            );
        }

        // Curve Evaluation
        {
            let parser_group = Some(parser.add_argument_group(
                "Curve Evaluation Options",
                "Parameters controlling curve evaluation behavior."
            ));
            parser.add_argument(
                "-e", "--error",
                "The error threshold (defaults to 1.0)",
                "PIXELS",
                Box::new(|dest_data, my_args| {
                    match f64::from_str(&my_args[0]) {
                        Ok(v) => {
                            dest_data.error_threshold = v;
                            return Ok(1);
                        },
                        Err(e) => {
                            return Err(e.to_string());
                        },
                    }
                }),
                1, argparse::ARGDEF_DEFAULT,
                parser_group,
            );
            parser.add_argument(
                "-t", "--simplify",
                "Simplify polygon before fitting (defaults to 2.0)",
                "PIXELS",
                Box::new(|dest_data, my_args| {
                    match f64::from_str(&my_args[0]) {
                        Ok(v) => {
                            dest_data.simplify_threshold = v;
                            return Ok(1);
                        },
                        Err(e) => {
                            return Err(e.to_string());
                        },
                    }
                }),
                1, argparse::ARGDEF_DEFAULT,
                parser_group,
            );


            parser.add_argument(
                "-c", "--corner",
                "The corner threshold (`pi` or greater to disable, defaults to 30.0)",
                "DEGREES",
                Box::new(|dest_data, my_args| {
                    match f64::from_str(&my_args[0]) {
                        Ok(v) => {
                            dest_data.corner_threshold = v.to_radians();
                            return Ok(1);
                        },
                        Err(e) => {
                            return Err(e.to_string());
                        },
                    }
                }),
                1, argparse::ARGDEF_DEFAULT,
                parser_group,
            );
            parser.add_argument(
                "", "--optimize-exhaustive",
                "When passed, perform exhaustive curve fitting (can be slow!)",
                "",
                Box::new(|dest_data, _my_args| {
                    dest_data.use_optimize_exhaustive = true;
                    return Ok(0);
                }),
                0, argparse::ARGDEF_DEFAULT,
                parser_group,
            );
        }

        // Output Options
        {
            let parser_group = Some(parser.add_argument_group(
                "Output Options",
                "Generic options for output (format agnostic)."
            ));
            parser.add_argument(
                "-s", "--scale",
                "Scale for output, (defaults to 1).",
                "SCALE",
                Box::new(|dest_data, my_args| {
                    match f64::from_str(&my_args[0]) {
                        Ok(v) => {
                            dest_data.output_scale = v;
                            return Ok(1);
                        },
                        Err(e) => {
                            return Err(e.to_string());
                        },
                    }
                }),
                1, argparse::ARGDEF_DEFAULT,
                parser_group,
            );
            parser.add_argument(
                "-p", "--passes",
                concat!("Write extra debug graphics, comma separated list of passes including ",
                        "[PIXEL, PRE_FIT, TANGENT], ",
                        "(defaults to [])."),
                "PASSES",
                Box::new(|dest_data, my_args| {
                    for pass_string in my_args[0].split(",") {
                        match pass_string.as_ref() {
                            "PIXEL" => {
                                dest_data.debug_passes |= debug_pass::kind::PIXEL;
                            }
                            "PRE_FIT" => {
                                dest_data.debug_passes |= debug_pass::kind::PRE_FIT;
                            }
                            "TANGENT" => {
                                dest_data.debug_passes |= debug_pass::kind::TANGENT;
                            }
                            _ => {
                                return Err(format!(
                                    "Expected [PIXEL, PRE_FIT, TANGENT], not '{}'",
                                    my_args[0],
                                ));
                            }
                        }
                    }
                    return Ok(1);
                }),
                1, argparse::ARGDEF_DEFAULT,
                parser_group,
            );
            parser.add_argument(
                "", "--pass-scale",
                "Scale graphic details used in some debug passes, (defaults to 1).",
                "SCALE",
                Box::new(|dest_data, my_args| {
                    match f64::from_str(&my_args[0]) {
                        Ok(v) => {
                            dest_data.debug_pass_scale = v;
                            return Ok(1);
                        },
                        Err(e) => {
                            return Err(e.to_string());
                        },
                    }
                }),
                1, argparse::ARGDEF_DEFAULT,
                parser_group,
            );
        }

        parser.add_argument(
            "-h", "--help",
            "Print help text",
            "",
            Box::new(|dest_data, _my_args| {
                dest_data.show_help = true;
                return Ok(0);
            }),
            0, argparse::ARGDEF_DEFAULT,
            None,
        );

        let args: Vec<String> = ::std::env::args().collect();
        let result = parser.parse(&args[1..]);

        if parser.dest_data.show_help {
            parser.print_help();
            return;
        }

        match result {
            Ok(()) => {}
            Err(e) => {
                use std::io::Write;
                writeln!(&mut std::io::stderr(), "{}, aborting!", e.to_string()).unwrap();
                std::process::exit(1);
            }
        }
    }

    match ::intern::image_load::from_filepath_any(&trace_params.input_filepath) {
        Ok((size, color_max, pixel_buffer)) => {
            println!("{:?} {}", size, color_max);
            let mut image: Vec<bool> = vec![false; pixel_buffer.len()];
            let color_mid = ((color_max / 2) as u32) * 3;
            for (p_src, p_dst) in pixel_buffer.iter().zip(&mut image) {
                let t = (p_src[0] as u32) +
                        (p_src[1] as u32) +
                        (p_src[2] as u32);
                if t < color_mid {
                    *p_dst = true;
                }
            }

            if trace_params.mode == curve_fit_nd::TraceMode::Centerline {
                use image_skeletonize;
                image_skeletonize::calculate(&mut image, &[size[0], size[1]]);
            }

            match trace_image(
                &trace_params.output_filepath,
                trace_params.output_scale,
                &image.as_slice(),
                &size,
                trace_params.error_threshold,
                trace_params.simplify_threshold,
                trace_params.corner_threshold,
                trace_params.use_optimize_exhaustive,
                0.75,
                trace_params.mode,
                trace_params.turn_policy,
                trace_params.debug_passes,
                trace_params.debug_pass_scale * trace_params.output_scale,
                )
            {
                Ok(()) => {}
                Err(e) => {
                    println!("Error writing output {:?}", e);
                }
            }
        }
        Err(e) => {
            println!("Error reading PPM {:?}", e);
        }
    }
}

#[cfg(test)]
#[path="tests.rs"] mod test;
