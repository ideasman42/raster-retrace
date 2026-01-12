fn main() {
    // Fixed to 2D for raster-retrace
    println!("cargo:rustc-env=DIMS_MAX=2");
}
