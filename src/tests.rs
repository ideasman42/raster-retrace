
use curve_fit_nd;
use polys_from_raster_outline;

macro_rules! test_image {
    ($id:ident, $size:expr, $error:expr, $corner_angle:expr, $length:expr, $image:expr) => {
        #[test]
        fn $id() {
            static IMAGE: &'static [bool] = $image;
            let size = $size;
            debug_assert!(IMAGE.len() == (size[0] * size[1]));
            match ::trace_image(
                &String::from(concat!(stringify!($id), ".svg")),
                1.0, IMAGE, &size, $error, $length, $corner_angle, false,
                0.75,
                curve_fit_nd::TraceMode::Outline,
                polys_from_raster_outline::TurnPolicy::Majority,
                0,
            ) {
                Ok(_) => (),
                Err(e) => println!("Error {:?}", e),
            }
        }
    }
}

test_image!(
    test_image_small,
    [10, 10], 0.75, ::std::f64::consts::PI / 8.0, 0.2, &[
    false, false, false, false, true,  true,  false, false, true,  false,
    false, true,  true,  false, false, true,  false, true,  true,  false,
    false, true,  true,  false, false, false, false, true,  true,  false,
    false, true,  true,  true,  false, false, false, true,  true,  false,
    false, false, true,  true,  true,  true,  true,  true,  true,  false,
    false, false, true,  true,  true,  true,  true,  true,  true,  true,
    false, false, true,  true,  false, false, true,  true,  true,  true,
    false, false, true,  true,  false, false, true,  true,  true,  true,
    false, false, true,  true,  false, false, false, true,  true,  true,
    false, false, true,  true,  false, false, false, true,  true,  false,
    ]);

