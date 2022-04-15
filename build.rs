extern crate bindgen;
extern crate cc;
use std::env;
use std::path::PathBuf;

fn main() {
    cc::Build::new()
    // .files("third_party/fastme/")
    .file("third_party/fastme/BIONJ.c")
.file("third_party/fastme/MVR.c")
.file("third_party/fastme/NNI.c")
.file("third_party/fastme/SPR.c")
.file("third_party/fastme/bNNI.c")
.file("third_party/fastme/bme.c")
.file("third_party/fastme/distance.c")
.file("third_party/fastme/fastme.c")
.file("third_party/fastme/gme.c")
.file("third_party/fastme/graph.c")
.file("third_party/fastme/heap.c")
.file("third_party/fastme/inputs.c")
.file("third_party/fastme/interface_options.c")
.file("third_party/fastme/interface_utilities.c")
.file("third_party/fastme/newick.c")
.file("third_party/fastme/p_bootstrap.c")
.file("third_party/fastme/p_eigen.c")
.file("third_party/fastme/p_lk.c")
.file("third_party/fastme/p_models.c")
.file("third_party/fastme/p_optimiz.c")
.file("third_party/fastme/p_utils.c")
.file("third_party/fastme/random.c")
.file("third_party/fastme/traverse.c")
.file("third_party/fastme/utils.c")
        .compile("fastme");
    println!("cargo:rustc-link-search=native=third_party/fastme/");
    println!("cargo:rustc-link-lib=fastme");
    // Tell cargo to invalidate the built crate whenever the wrapper changes
    println!("cargo:rerun-if-changed=wrapper.h");
    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("wrapper.h")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .derive_default(true)
        .blocklist_item("FP_NAN")
        .blocklist_item("FP_INFINITE")
        .blocklist_item("FP_ZERO")
        .blocklist_item("FP_SUBNORMAL")
        .blocklist_item("FP_NORMAL")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}