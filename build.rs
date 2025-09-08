use std::env;
use std::path::PathBuf;

fn main() {
    // Path to the installed Sundials includes and libs - adjust as needed.
    let include_dir = "/Users/kmaitreys/Documents/software/cvode-7.4.0/install/include";
    let lib_dir = "/Users/kmaitreys/Documents/software/cvode-7.4.0/install/lib";

    println!("cargo:rerun-if-changed={}", include_dir);

    // Link against Sundials libraries used by the example
    println!("cargo:rustc-link-search=native={}", lib_dir);
    println!("cargo:rustc-link-lib=sundials_cvode");
    println!("cargo:rustc-link-lib=sundials_nvecserial");
    println!("cargo:rustc-link-lib=sundials_nvecmanyvector");
    println!("cargo:rustc-link-lib=sundials_core");

    // Configure bindgen
    let bindings = bindgen::Builder::default()
        // Use a header that declares the API you want to bind.
        // Prefer creating a small header `robertson_bind.h` that includes necessary sundials headers
        .header("wrapper.h")
        // Tell clang where to find Sundials headers
        .clang_arg(format!("-I{}", include_dir))
        // If SUNDIALS uses macros to select real type, define appropriate macros:
        // .clang_arg("-DSUN_REALTYPE=double") or whichever is required.
        // .clang_arg("-DSUN_PRECISION=double") // example if needed
        // Limit the generated bindings to what you need:
        .derive_default(true)
        .derive_debug(true)
        .generate_comments(true)
        .generate()
        .expect("bindgen failed");

    // Write the bindings to $OUT_DIR/bindings.rs
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("failed to write bindings");
}
