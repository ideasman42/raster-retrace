///
/// Generalizes image loading.
///

mod image_load_ppm;

use ::std::io::{
    Error,
    ErrorKind,
};

#[derive(PartialEq, Debug, Copy, Clone)]
pub enum ImageFormat {
    PPM,
    // PNG,
}

fn format_from_filepath(
    filepath: &String,
) -> Option<ImageFormat> {
    if filepath.ends_with(".ppm") {
        return Some(ImageFormat::PPM);
    // } else if filepath.ends_with(".png") {
    //     return Some(ImageFormat::PNG);
    } else {
        return None;
    }
}

pub fn from_filepath_format(
    filepath: &String,
    format: ImageFormat,
) -> Result<([usize; 2], usize, Vec<[u8; 3]>), Error> {
    if format == ImageFormat::PPM {
        let file = ::std::fs::File::open(filepath).expect("open failed");
        return image_load_ppm::from_file(&file);
    // } else if format == ImageFormat::PNG {
    //     return image_load_png::from_filepath(filepath);
    }
    return Err(Error::new(ErrorKind::Other, "Unknown file format"));
}

pub fn from_filepath_any(
    filepath: &String,
) -> Result<([usize; 2], usize, Vec<[u8; 3]>), Error> {
    if let Some(format) = format_from_filepath(filepath) {
        return from_filepath_format(filepath, format);
    }
    return Err(Error::new(ErrorKind::Other, "Unknown file format"));
}

