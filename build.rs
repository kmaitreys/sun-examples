fn main() {
    // Adjust this once if your SUNDIALS install moves
    let sundials_prefix = "/Users/kmaitreys/Documents/software/cvode-7.4.0/install";
    let examples_dir = "/Users/kmaitreys/Documents/probe/sun-examples/examples";

    // Add all C example files here
    let sources = ["cv_roberts_dns.c", "cv_diurnal_kry.c"];

    for src in &sources {
        println!("cargo:rerun-if-changed={}/{}", examples_dir, src);
    }

    let mut build = cc::Build::new();
    for src in &sources {
        build.file(format!("{}/{}", examples_dir, src));
    }
    build
        .include(format!("{}/include", sundials_prefix))
        .flag_if_supported("-std=c11")
        .compile("robertson");

    println!("cargo:rustc-link-search=native={}/lib", sundials_prefix);

    // Required SUNDIALS components
    println!("cargo:rustc-link-lib=sundials_core");
    println!("cargo:rustc-link-lib=sundials_cvode");
    println!("cargo:rustc-link-lib=sundials_nvecserial");
    println!("cargo:rustc-link-lib=sundials_sunmatrixdense");
    println!("cargo:rustc-link-lib=sundials_sunlinsoldense");

    // system math lib
    println!("cargo:rustc-link-lib=m");
}
