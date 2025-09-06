use sundials_sys::*;

extern "C" fn f(
    _t: sunrealtype,
    _y: N_Vector,
    _ydot: N_Vector,
    _user_data: *mut std::os::raw::c_void,
) -> i32 {
    0
}

fn main() {
    todo!()
}
