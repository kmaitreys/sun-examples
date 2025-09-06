#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

use std::ffi::{c_int, c_void};
use std::ptr;
use sundials_sys::*;

pub type sunbooleantype = i32;
pub const SUNTRUE: sunbooleantype = 1;
pub const SUNFALSE: sunbooleantype = 0;

const NEQ: sunindextype = 2;
const FORCE: f64 = 2200.0;
const MASSR: f64 = 10.0;
const MASSF0: f64 = 1.0;
const BRATE: f64 = 0.1;
const DRAG: f64 = 0.3;
const GRAV: f64 = 32.0;
const HCUT: f64 = 4000.0;

const Y1: f64 = 0.0;
const Y2: f64 = 0.0;
const RTOL: f64 = 1e-5;
const ATOL1: f64 = 1e-2;
const ATOL2: f64 = 1e-1;
const T0: f64 = 0.0;
const T1: f64 = 1.0;
const TINC: f64 = 1.0;
const NOUT: i32 = 70;

extern "C" fn f(t: sunrealtype, y: N_Vector, ydot: N_Vector, user_data: *mut c_void) -> i32 {
    unsafe {
        let engine_on = *(user_data as *mut sunbooleantype);
        let ydata = N_VGetArrayPointer(y);
        let yddata = N_VGetArrayPointer(ydot);
        let v = *ydata.add(1);
        *yddata.add(0) = v;
        // dbg!(*ydata, *yddata);
        let acc = if engine_on == SUNTRUE {
            FORCE / (MASSR + MASSF0 - BRATE * t)
        } else {
            0.0
        };
        *yddata.add(1) = acc - DRAG * v - GRAV;
    }
    0
}

extern "C" fn jac(
    _t: sunrealtype,
    _y: N_Vector,
    _fy: N_Vector,
    j: SUNMatrix,
    _user_data: *mut c_void,
    _tmp1: N_Vector,
    _tmp2: N_Vector,
    _tmp3: N_Vector,
) -> i32 {
    unsafe {
        let jdata = SUNDenseMatrix_Data(j);
        *jdata.add(1) = 1.0;
        *jdata.add(3) = -DRAG;
    }
    0
}

extern "C" fn g(
    t: sunrealtype,
    y: N_Vector,
    gout: *mut sunrealtype,
    user_data: *mut c_void,
) -> i32 {
    unsafe {
        let engine_on = *(user_data as *mut sunbooleantype);
        let ydata = N_VGetArrayPointer(y);
        if engine_on == SUNTRUE {
            *gout.add(0) = MASSF0 - BRATE * t;
            *gout.add(1) = *ydata.add(0) - HCUT;
        } else {
            *gout.add(0) = *ydata.add(1);
        }
    }
    0
}

// ---- Print helpers ----
unsafe fn print_output(t: f64, y1: f64, y2: f64) {
    println!("At t = {:0.4e}      y ={:14.6e}  {:14.6e}", t, y1, y2);
}

unsafe fn print_root_info(root_f1: i32, root_f2: i32, numroot: i32) {
    if numroot == 2 {
        println!("    rootsfound[] = {:3} {:3}", root_f1, root_f2);
    } else {
        println!("    rootsfound[] = {:3}", root_f1);
    }
}

unsafe fn print_final_stats(cvode_mem: *mut c_void) {
    let mut nst: i64 = 0;
    let mut nfe: i64 = 0;
    let mut nsetups: i64 = 0;
    let mut nje: i64 = 0;
    let mut nfeLS: i64 = 0;
    let mut nni: i64 = 0;
    let mut ncfn: i64 = 0;
    let mut netf: i64 = 0;
    let mut nge: i64 = 0;

    unsafe {
        CVodeGetNumSteps(cvode_mem, &mut nst);
        CVodeGetNumRhsEvals(cvode_mem, &mut nfe);
        CVodeGetNumLinSolvSetups(cvode_mem, &mut nsetups);
        CVodeGetNumErrTestFails(cvode_mem, &mut netf);
        CVodeGetNumNonlinSolvIters(cvode_mem, &mut nni);
        CVodeGetNumNonlinSolvConvFails(cvode_mem, &mut ncfn);
        CVodeGetNumJacEvals(cvode_mem, &mut nje);
        CVodeGetNumLinRhsEvals(cvode_mem, &mut nfeLS);
        CVodeGetNumGEvals(cvode_mem, &mut nge);
    }

    println!("\nFinal Statistics:");
    println!(
        "nst = {:<6} nfe  = {:<6} nsetups = {:<6} nfeLS = {:<6} nje = {:<6}",
        nst, nfe, nsetups, nfeLS, nje
    );
    println!(
        "nni = {:<6} ncfn = {:<6} netf = {:<6} nge = {}",
        nni, ncfn, netf, nge
    );
}

fn main() {
    unsafe {
        let mut sunctx: SUNContext = ptr::null_mut();
        let mut retval = SUNContext_Create(SUN_COMM_NULL, &mut sunctx);
        assert_eq!(retval, 0);

        let y = N_VNew_Serial(NEQ, sunctx);
        let abstol = N_VNew_Serial(NEQ, sunctx);

        let ydata = N_VGetArrayPointer(y);
        *ydata.add(0) = Y1;
        *ydata.add(1) = Y2;

        let adata = N_VGetArrayPointer(abstol);
        *adata.add(0) = ATOL1;
        *adata.add(1) = ATOL2;

        let cvode_mem = CVodeCreate(CV_BDF, sunctx);
        CVodeInit(cvode_mem, Some(f), T0, y);
        CVodeSVtolerances(cvode_mem, RTOL, abstol);

        let mut engine_on: Box<sunbooleantype> = Box::new(SUNTRUE);
        CVodeSetUserData(cvode_mem, &mut *engine_on as *mut _ as *mut c_void);
        CVodeRootInit(cvode_mem, 2, Some(g));

        let a = SUNDenseMatrix(NEQ, NEQ, sunctx);
        let ls = SUNLinSol_Dense(y, a, sunctx);
        CVodeSetLinearSolver(cvode_mem, ls, a);
        let flag = CVodeSetJacFn(cvode_mem, Some(jac));
        assert_eq!(flag, 0, "Jacobian function registration failed");

        println!("\nAccelerating rocket problem\n");

        let mut iout = 0;
        let mut tout = T1;
        let mut t: f64 = 0.0;
        let mut numroot = 2;
        let mut rootsfound: Vec<c_int> = vec![0; numroot as usize];

        loop {
            retval = CVode(cvode_mem, tout, y, &mut t, CV_NORMAL);
            if retval < 0 {
                break;
            }

            let ydata = N_VGetArrayPointer(y);
            print_output(t, *ydata.add(0), *ydata.add(1));

            if *engine_on == SUNTRUE && retval == CV_ROOT_RETURN {
                println!("here");
                CVodeGetRootInfo(cvode_mem, rootsfound.as_mut_ptr());
                print_root_info(rootsfound[0], rootsfound[1], numroot);
                *engine_on = SUNFALSE;
                numroot = 1;
                CVodeRootInit(cvode_mem, 1, Some(g));
                retval = CVodeReInit(cvode_mem, t, y);
            } else if *engine_on == SUNFALSE && retval == CV_ROOT_RETURN {
                retval = CVodeGetRootInfo(cvode_mem, rootsfound.as_mut_ptr());
                print_root_info(rootsfound[0], rootsfound[1], numroot);
            }

            if retval == CV_SUCCESS {
                iout += 1;
                tout += TINC;
            }
            if iout == NOUT {
                break;
            }
            if *ydata.add(0) < 0.0 {
                break;
            }
        }

        print_final_stats(cvode_mem);

        N_VDestroy(y);
        N_VDestroy(abstol);
        CVodeFree(&cvode_mem as *const _ as *mut _);
        SUNLinSolFree(ls);
        SUNMatDestroy(a);
        SUNContext_Free(&mut sunctx);
    }
}
