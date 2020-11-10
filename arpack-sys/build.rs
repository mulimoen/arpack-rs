fn main() {
    let builder = bindgen::Builder::default();
    let bindings = builder
        .header("wrapper.h")
        .generate()
        .expect("Failed to generate bindigns");

    let out_path = std::path::PathBuf::from(std::env::var_os("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("arpack_bindings.rs"))
        .expect("Failed to write bindings");

    println!("cargo:rustc-link-lib=arpack");
}
