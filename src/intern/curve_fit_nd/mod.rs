
mod curve_fit_from_polys;
mod curve_fit_single;

// we could make this a separate module
pub use ::intern::math_vector;

pub use self::curve_fit_from_polys::{
    TraceMode,
    fit_poly_single,
    fit_poly_list,
};

