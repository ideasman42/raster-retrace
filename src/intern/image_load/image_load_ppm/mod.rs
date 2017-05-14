
///
/// Module for reading image data from files.
///


/// TODO
///
/// - 16bpc PPM files.
///   not really that hard, but also not that interesting.
/// - More efficient vector reading (could be a single operation).


macro_rules! elem {
    ($val:expr, $($var:expr), *) => {
        $($val == $var) || *
    }
}

/// Returns (size, pixel_data), or fail.

use ::std::io::{
    Error,
    ErrorKind,
    SeekFrom,
};

use std::io::prelude::*;
use std::str::FromStr;

pub fn from_file(
    mut f: &::std::fs::File,
) -> Result<([usize; 2], usize, Vec<[u8; 3]>), Error> {

    fn read_until_newline(
        mut f: &::std::fs::File,
    ) -> Result<(), Error> {
        let mut buf: [u8; 1] = [0];
        loop {
            f.read_exact(&mut buf)?;
            if buf[0] == '\n' as u8 {
                break;
            }
        }
        Ok(())
    }

    fn read_peek_byte(
        mut f: &::std::fs::File,
    ) -> Result<u8, Error> {
        let mut buf: [u8; 1] = [0];
        f.read_exact(&mut buf)?;
        f.seek(SeekFrom::Current(-1))?;
        return Ok(buf[0]);
    }

    fn read_as_usize_skip_ws(
        mut f: &::std::fs::File,
    ) -> Result<usize, Error> {
        // note, we could attempt to evaluate this as bytes
        // (atio style). for now it seems Rust's std lib doesn't support this.
        let mut num_str = String::with_capacity(16);
        let mut buf: [u8; 1] = [0];
        loop {
            f.read_exact(&mut buf)?;

            if elem!(buf[0], ' ' as u8, '\t' as u8, '\r' as u8, '\n' as u8) {
                if num_str.len() != 0 {
                    break;
                }
            } else {
                num_str.push(buf[0] as char);
            }

            if num_str.len() == 0 {
                return Err(Error::new(ErrorKind::Other, "No number found"));
            }
        }

        return match usize::from_str(num_str.as_str()) {
            Ok(n) => { Ok(n) }
            Err(e) => { Err(Error::new(ErrorKind::Other, e.to_string())) }
        };
    }

    // Header Magic
    {
        let mut header: [u8; 2] = [0; 2];
        f.read_exact(&mut header)?;
        if !(header[0] == 'P' as u8 && header[1] == '6' as u8) {
            return Err(Error::new(ErrorKind::Other, "Invalid header"));
        }
        read_until_newline(f)?;
    }

    // Header Content
    let mut size: [usize; 2] = [0; 2];
    let color_max;  // range is 1-65535
    loop {
        let byte = read_peek_byte(f)?;
        if elem!(byte, '#' as u8, ' ' as u8, '\t' as u8, '\r' as u8, '\n' as u8) {
            read_until_newline(f)?;
        } else {
            // check if size has been set
            if size[0] == 0 {
                size = [
                    read_as_usize_skip_ws(f)?,
                    read_as_usize_skip_ws(f)?,
                ];
                if !(size[0] > 0 && size[1] > 0) {
                    return Err(Error::new(ErrorKind::Other, "Invalid size"));
                }
            } else {
                color_max = read_as_usize_skip_ws(f)?;
                if !(color_max > 0 && color_max < 65536)  {
                    return Err(Error::new(ErrorKind::Other, "Invalid color range"));
                }
                // Nothing left to read,
                // we have a single whitespace character between this and the real data.
                // which we will have already read, so can jump directly into reading the data.
                break;
            }
        }
    }

    // All header data is read.

    // TODO, support allocation failure
    let pixel_buffer_len = size[0] * size[1];
    let mut pixel_buffer = Vec::<[u8; 3]>::with_capacity(pixel_buffer_len);
    let mut pixel: [u8; 3] = [0; 3];
    for _ in 0..pixel_buffer_len {
        f.read_exact(&mut pixel)?;
        pixel_buffer.push(pixel);
    }
    return Ok((size, color_max, pixel_buffer));
}

