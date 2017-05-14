// 3-D medial surface/axis thinning algorithms.
// Computer Vision, Graphics, and Image Processing, 56(6):462-478, 1994.

struct Bitmap<'a> {
    data: &'a mut Vec<bool>,
    size: [i32; 2],
}

pub fn calculate(
    data: &mut Vec<bool>,
    size: &[usize; 2],
) {
    compute_thin_image(
        &mut Bitmap {
            data: data,
            size: [
                size[0] as i32,
                size[1] as i32,
            ],
        }
    );
}

fn compute_thin_image(image: &mut Bitmap)
{
    let mut simple_border_points: Vec<[i32; 2]> = Vec::new();

    // Loop through the image several times until there is no change.
    let mut unchanged_borders = 0;

    // loop until no change for all the six border types
    while unchanged_borders < 4 {
        unchanged_borders = 0;
        for current_border in 0..4 {
            let mut no_change: bool = true;

            // Loop over each pixel
            for y in 0..image.size[1] {
                for x in 0..image.size[0] {
                    // check if point is foreground
                    if pixel_get_no_check(image, x, y) == false {
                        // current point is already background
                        continue;
                    }

                    // check 4-neighbors if point is a border point of type current_border
                    if match current_border {
                        0 => pixel_get(image, x, y - 1),
                        1 => pixel_get(image, x, y + 1),
                        2 => pixel_get(image, x + 1, y),
                        3 => pixel_get(image, x - 1, y),
                        _ => unreachable!() }
                    {
                        // current point is not deletable
                        continue;
                    }

                    if pixel_is_endpoint(image, x, y) {
                        continue;
                    }

                    let neighborhood = neighborhood_get_no_center(image, x, y);

                    // Check if point is Euler invariant (condition 1 in Lee[94])
                    if !is_euler_invariant(neighborhood) {
                        // current point is not deletable
                        continue;
                    }

                    // Check if point is simple
                    // (deletion does not change connectivity in the 3x3 neighborhood)
                    // (conditions 2 and 3 in Lee[94])
                    if !is_simple_point(neighborhood) {
                        // current point is not deletable
                        continue;
                    }
                    // add all simple border points to a list for sequential re-checking
                    simple_border_points.push([x, y]);
                }
            }

            // sequential re-checking to preserve connectivity when
            // deleting in a parallel way
            for index in &simple_border_points {
                // Check if border points is simple
                if is_simple_point(neighborhood_get_no_center(image, index[0], index[1])) {
                    // we can delete the current point
                    pixel_set(image, index[0], index[1], false);
                    no_change = false;
                }
            }

            if no_change {
                unchanged_borders += 1;
            }

            simple_border_points.clear();
        }
    }
}

/// Check if a point in the given stack is at the end of an arc.
/// return true if the point has exactly one neighbor
fn pixel_is_endpoint(image: &Bitmap, x: i32, y: i32) -> bool
{
    let mut number_of_neighbors: u32 = 0;
    let neighbors = neighborhood_get_no_center(image, x, y);
    for i in 0..DIR_FLAG_NUM_NO_CENTER {
        if neighbors & (1 << i) != 0 {
            number_of_neighbors += 1;
        }
    }
    // 2 and not 1 because the center pixel will be counted as well
    return number_of_neighbors == 1;
}

// const DIR_FLAG_NUM:           u32 = 9;
const DIR_FLAG_NUM_NO_CENTER: u32 = 8;

const DIR_SW: u32 = 1 <<  0;
const DIR_S:  u32 = 1 << 1;
const DIR_SE: u32 = 1 << 2;

const DIR_W:  u32 = 1 << 3;
const DIR_E:  u32 = 1 << 4;

const DIR_NW: u32 = 1 << 5;
const DIR_N:  u32 = 1 << 6;
const DIR_NE: u32 = 1 << 7;

// Currently unused
// const DIR_C: u32 = 1 << 8;

/// Get neighborhood of a pixel in a 2D image (0 border conditions)
///
/// return corresponding DIR_FLAG_NUM-pixels neighborhood (0 if out of image)
/*
fn neighborhood_get(image: &Bitmap, x: i32, y: i32) -> u32
{
    return (
        if pixel_get(image, x - 1, y - 1)     { DIR_SW      } else { 0 } |
        if pixel_get(image, x,     y - 1)     { DIR_S       } else { 0 } |
        if pixel_get(image, x + 1, y - 1)     { DIR_SE      } else { 0 } |

        if pixel_get(image, x - 1, y)         { DIR_W       } else { 0 } |
        if pixel_get(image, x + 1, y)         { DIR_E       } else { 0 } |

        if pixel_get(image, x - 1, y + 1)     { DIR_NW      } else { 0 } |
        if pixel_get(image, x,     y + 1)     { DIR_N       } else { 0 } |
        if pixel_get(image, x + 1, y + 1)     { DIR_NE      } else { 0 } |

        if pixel_get(image, x,     y)         { DIR_C       } else { 0 }
    );
}
*/

fn neighborhood_get_no_center(image: &Bitmap, x: i32, y: i32) -> u32
{
    return
        if pixel_get(image, x - 1, y - 1)     { DIR_SW      } else { 0 } |
        if pixel_get(image, x,     y - 1)     { DIR_S       } else { 0 } |
        if pixel_get(image, x + 1, y - 1)     { DIR_SE      } else { 0 } |

        if pixel_get(image, x - 1, y)         { DIR_W       } else { 0 } |
        if pixel_get(image, x + 1, y)         { DIR_E       } else { 0 } |

        if pixel_get(image, x - 1, y + 1)     { DIR_NW      } else { 0 } |
        if pixel_get(image, x,     y + 1)     { DIR_N       } else { 0 } |
        if pixel_get(image, x + 1, y + 1)     { DIR_NE      } else { 0 }
    ;
}

/// Get pixel in 2D image (0 border conditions)
fn pixel_get(image: &Bitmap, x: i32, y: i32) -> bool
{
    if x >= 0 && x < image.size[0] &&
       y >= 0 && y < image.size[1]
    {
        return image.data[(x + y * image.size[0]) as usize];
    } else {
        return false;
    }
}

/// Get pixel in 2D image (no border checking)
fn pixel_get_no_check(image: &Bitmap, x: i32, y: i32) -> bool
{
    return image.data[(x + y * image.size[0]) as usize];
}

/// Set pixel in 2D image
fn pixel_set(image: &mut Bitmap, x: i32, y: i32, value: bool)
{
    if x >= 0 && x < image.size[0] &&
       y >= 0 && y < image.size[1]
    {
        image.data[(x + y * image.size[0]) as usize] = value;
    }
}

const INDEX_LUT: [i32; 32] = [
    0,  1, 0, -1, 0, -1, 0,  1, 0, -3, 0, -1, 0, -1, 0,  1,
    0, -1, 0,  1, 0,  1, 0, -1, 0,  3, 0,  1, 0,  1, 0, -1,
];

/// Check if a point is Euler invariant
///
/// return true or false if the point is Euler invariant or not
fn is_euler_invariant(neighbors: u32) -> bool
{
    // Calculate Euler characteristic for each quadrant and sum up
    let mut euler_char: i32 = 0;

    euler_char += INDEX_LUT[index_quadrant_sw(neighbors) as usize];
    euler_char += INDEX_LUT[index_quadrant_se(neighbors) as usize];
    euler_char += INDEX_LUT[index_quadrant_nw(neighbors) as usize];
    euler_char += INDEX_LUT[index_quadrant_ne(neighbors) as usize];

    return euler_char == 0;
}

fn index_quadrant_ne(neighbors: u32) -> u8
{
    return
        if neighbors & DIR_S != 0 { (1 << 4) } else { 0 } |
        if neighbors & DIR_E != 0 { (1 << 1) } else { 0 } |
        (1 << 0)
    ;
}

fn index_quadrant_nw(neighbors: u32) -> u8
{
    return
        if neighbors & DIR_W != 0 { (1 << 4) } else { 0 } |
        if neighbors & DIR_S != 0 { (1 << 2) } else { 0 } |
        (1 << 0)
    ;
}

fn index_quadrant_se(neighbors: u32) -> u8
{
    return
        if neighbors & DIR_N != 0 { (1 << 4) } else { 0 } |
        if neighbors & DIR_E != 0 { (1 << 1) } else { 0 } |
        (1 << 0)
    ;
}

fn index_quadrant_sw(neighbors: u32) -> u8
{
    return
        if neighbors & DIR_N != 0 { (1 << 4) } else { 0 } |
        if neighbors & DIR_W != 0 { (1 << 2) } else { 0 } |
        (1 << 0)
    ;
}

/// Check if current point is a Simple Point.
/// This method is named 'N(v)_labeling' in [Lee94].
/// Outputs the number of connected objects in a neighborhood of a point
/// after this point would have been removed.
///
/// * `neighbors` - neighbors neighbor pixels of the point.
///
/// Return true or false if the point is simple or not.
fn is_simple_point(neighbors: u32) -> bool
{
    let mut quad: u32 = neighbors;

    // set initial label
    let mut label: u32 = 0;
    // for all points in the neighborhood
    for i in 0..DIR_FLAG_NUM_NO_CENTER {
        // voxel has not been labeled yet
        let flag = 1 << i;
        if (quad & flag) != 0 {
            // Start recursion with any quadrant that contains the point 'i'.
            //
            // set points in this quadrant to current label
            // and recursive labeling of adjacent quadrants.
            match flag {
                DIR_SW | DIR_S | DIR_W => {
                    quadtree_labeling_sw(&mut quad);
                }
                DIR_SE | DIR_E => {
                    quadtree_labeling_se(&mut quad);
                }
                DIR_NW | DIR_N => {
                    quadtree_labeling_nw(&mut quad);
                }
                DIR_NE => {
                    quadtree_labeling_ne(&mut quad);
                }
                _ => {
                    unreachable!();
                }
            }
            label += 1;
            if label >= 2 {
                return false;
            }
        }
    }
    // return label-2; in [Lee94] if the number of connected components would be needed
    return true;
}

/// This is a recursive method that calculates the number of connected
/// components in the 2D neighborhood after the center pixel would
/// have been removed.

fn quadtree_labeling_sw(quad: &mut u32)
{
    if *quad & DIR_SW != 0 {
        *quad &= !DIR_SW;
    }
    if *quad & DIR_S != 0 {
        *quad &= !DIR_S;
        quadtree_labeling_se(quad);
    }
    if *quad & DIR_W != 0 {
        *quad &= !DIR_W;
        quadtree_labeling_nw(quad);
    }
}

fn quadtree_labeling_se(quad: &mut u32)
{
    if *quad & DIR_S != 0 {
        *quad &= !DIR_S;
        quadtree_labeling_sw(quad);
        quadtree_labeling_se(quad);
    }
    if *quad & DIR_SE != 0 {
        *quad &= !DIR_SE;
        quadtree_labeling_se(quad);
    }
    if *quad & DIR_E != 0 {
        *quad &= !DIR_E;
        quadtree_labeling_ne(quad);
        quadtree_labeling_se(quad);
    }
}

fn quadtree_labeling_nw(quad: &mut u32)
{
    if *quad & DIR_W != 0 {
        *quad &= !DIR_W;
        quadtree_labeling_sw(quad);
        quadtree_labeling_nw(quad);
    }
    if *quad & DIR_NW != 0 {
        *quad &= !DIR_NW;
        quadtree_labeling_nw(quad);
    }
    if *quad & DIR_N != 0 {
        *quad &= !DIR_N;
        quadtree_labeling_ne(quad);
        quadtree_labeling_nw(quad);
    }
}

fn quadtree_labeling_ne(quad: &mut u32)
{
    if *quad & DIR_E != 0 {
        *quad &= !DIR_E;
        quadtree_labeling_se(quad);
        quadtree_labeling_ne(quad);
    }
    if *quad & DIR_N != 0 {
        *quad &= !DIR_N;
        quadtree_labeling_nw(quad);
        quadtree_labeling_ne(quad);
    }
    if *quad & DIR_NE != 0 {
        *quad &= !DIR_NE;
        quadtree_labeling_ne(quad);
    }
}

