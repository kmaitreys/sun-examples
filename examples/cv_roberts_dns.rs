unsafe extern "C" {
    fn run_robertson(y1: *mut f64, y2: *mut f64, y3: *mut f64) -> i32;
}

fn main() {
    let mut y1 = 0.0;
    let mut y2 = 0.0;
    let mut y3 = 0.0;

    let ret = unsafe { run_robertson(&mut y1, &mut y2, &mut y3) };

    println!("return code = {}", ret);
    println!("y = [{:.6e}, {:.6e}, {:.6e}]", y1, y2, y3);
}
