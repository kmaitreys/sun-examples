unsafe extern "C" {
    fn run_diurnal_kry() -> i32;
}

fn main() {
    unsafe { run_diurnal_kry() };
}
