
mod curve_fit_from_polys;
mod curve_fit_single;
mod curve_fit_cubic_refit;

// we could make this a separate module
pub use ::intern::math_vector;

pub use self::curve_fit_from_polys::{
    TraceMode,
    fit_poly_single,
    fit_poly_list,
};

// Re-export refit types for external use
pub use self::curve_fit_cubic_refit::{
    Knot,
    PointData,
    refine_remove,
    refine_refit,
    refine_corner,
};

