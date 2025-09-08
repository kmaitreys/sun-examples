#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use std::{ffi::c_void, ptr};
use sundials_sys::*;

pub type sunbooleantype = i32;

/// Get the content of a serial vector.
///
/// # Safety
/// - `v` must be a valid non-null pointer to an `N_Vector`.
/// - `v` must actually contain serial vector content (`N_VectorContent_Serial`).
/// - Caller is responsible for ensuring the pointer is not dangling or invalid.
#[inline]
pub unsafe fn nv_content_s(v: N_Vector) -> N_VectorContent_Serial {
    unsafe { (*v).content as N_VectorContent_Serial }
}

/// Get the length of a serial vector.
///
/// # Safety
/// - `v` must be a valid non-null pointer to an `N_Vector`.
/// - `v` must contain valid serial vector content.
/// - Dereferencing an invalid or mismatched vector will cause undefined behavior.
#[inline]
pub unsafe fn nv_length_s(v: N_Vector) -> sunindextype {
    unsafe { (*nv_content_s(v)).length }
}

/// Check whether the serial vector owns its data array.
///
/// # Safety
/// - `v` must be a valid non-null pointer to an `N_Vector`.
/// - `v` must contain valid serial vector content.
/// - Caller must not assume memory ownership rules beyond what this flag indicates.
#[inline]
pub unsafe fn nv_own_data_s(v: N_Vector) -> sunbooleantype {
    unsafe { (*nv_content_s(v)).own_data }
}

/// Get a raw pointer to the data array of a serial vector.
///
/// # Safety
/// - `v` must be a valid non-null pointer to an `N_Vector`.
/// - Returned pointer must not outlive the vector.
/// - Caller must ensure proper bounds when indexing into this array.
/// - No aliasing or concurrent mutation rules are enforced here.
#[inline]
pub unsafe fn nv_data_s(v: N_Vector) -> *mut sunrealtype {
    unsafe { (*nv_content_s(v)).data }
}

/// Get the i-th element (zero-based) from the serial vector.
///
/// # Safety
/// - `v` must be a valid non-null pointer to an `N_Vector`.
/// - Index `i` must be less than the vector length.
/// - Caller must ensure no out-of-bounds access.
/// - Dereferencing an invalid pointer results in undefined behavior.
#[inline]
pub unsafe fn nv_ith_s(v: N_Vector, i: usize) -> *mut sunrealtype {
    unsafe { nv_data_s(v).add(i) }
}

extern "C" fn f(t: sunrealtype, y: N_Vector, ydot: N_Vector, user_data: *mut c_void) -> i32 {
    unsafe {
        let rdata = user_data as *mut sunrealtype;
        let lambda = *rdata;
        let u = *nv_ith_s(y, 0);

        *nv_ith_s(ydot, 0) = lambda * u + 1.0 / (1.0 + t * t) - lambda * t.atan();
    }
    0
}

fn MatrixEmbeddedLS(cvode_mem: *mut c_void, sunctx: SUNContext) -> SUNLinearSolver {
    unsafe {
        let LS = SUNLinSolNewEmpty(sunctx);
        assert!(!LS.is_null());

        (*(*LS).ops).gettype = Some(MatrixEmbeddedLSType);
        (*(*LS).ops).solve = Some(MatrixEmbeddedLSSolve);
        (*(*LS).ops).free = Some(MatrixEmbeddedLSFree);

        (*LS).content = cvode_mem;
        LS
    }
}

extern "C" fn MatrixEmbeddedLSType(_S: SUNLinearSolver) -> SUNLinearSolver_Type {
    SUNLinearSolver_Type_SUNLINEARSOLVER_MATRIX_EMBEDDED
}

extern "C" fn MatrixEmbeddedLSSolve(
    LS: SUNLinearSolver,
    A: SUNMatrix,
    x: N_Vector,
    b: N_Vector,
    tol: f64,
) -> i32 {
    unsafe {
        let mut tcur: sunrealtype = 0.0;
        let mut gamma: sunrealtype = 0.0;
        let mut rl1: sunrealtype = 0.0;
        let mut y: N_Vector = std::ptr::null_mut();
        let mut ypred: N_Vector = std::ptr::null_mut();
        let mut fn_: N_Vector = std::ptr::null_mut();
        let mut zn1: N_Vector = std::ptr::null_mut();
        let mut user_data: *mut c_void = std::ptr::null_mut();
        let retval = CVodeGetNonlinearSystemData(
            (*LS).content,
            &mut tcur,
            &mut ypred,
            &mut y,
            &mut fn_,
            &mut gamma,
            &mut rl1,
            &mut zn1,
            &mut user_data,
        );
        if retval != SUN_SUCCESS {
            return -1;
        }
        let rdata = user_data as *mut sunrealtype;
        let lambda = *rdata;
        *nv_ith_s(x, 0) = *nv_ith_s(b, 0) / (1.0 - gamma * lambda);

        SUN_SUCCESS
    }
}

unsafe extern "C" fn MatrixEmbeddedLSFree(LS: SUNLinearSolver) -> SUNErrCode {
    unsafe {
        if LS.is_null() {
            return SUN_SUCCESS;
        }
        (*LS).content = std::ptr::null_mut();
        SUNLinSolFreeEmpty(LS);
        return SUN_SUCCESS;
    }
}

fn main() {
    unsafe {
        const T0: f64 = 0.0;
        const Tf: f64 = 10.0;
        const dTout: f64 = 1.0;
        const NEQ: sunindextype = 1;
        const reltol: f64 = 1.0e-6;
        const abstol: f64 = 1.0e-10;
        const lambda: f64 = -100.0;

        let mut sunctx: SUNContext = ptr::null_mut();
        let mut retval = SUNContext_Create(SUN_COMM_NULL, &mut sunctx);
        assert!(retval == 0);

        println!("\nAnalytical ODE test problem:\n");
        println!("   lambda = {}", lambda);
        println!("   reltol = {:.1e}", reltol);
        println!("   abstol = {:.1e}", abstol);

        let y = N_VNew_Serial(NEQ, sunctx);
        assert!(!y.is_null());

        let cvode_mem = CVodeCreate(CV_BDF, sunctx);
        assert!(!cvode_mem.is_null());

        retval = CVodeInit(cvode_mem, Some(f), T0, y);
        assert!(retval == 0);

        retval = CVodeSetUserData(cvode_mem, &lambda as *const _ as *mut c_void);
        assert!(retval == 0);

        retval = CVodeSStolerances(cvode_mem, reltol, abstol);
        assert!(retval == 0);

        let LS = MatrixEmbeddedLS(cvode_mem, sunctx);
        assert!(!LS.is_null());

        retval = CVodeSetLinearSolver(cvode_mem, LS, std::ptr::null_mut());
        assert!(retval == 0);

        let mut t = T0;
        let mut tout = T0 + dTout;
        println!("        t           u\n");
        println!("   ---------------------\n");

        while Tf - t > 1.0e-15 {
            retval = CVode(cvode_mem, tout, y, &mut t, CV_NORMAL);
            if retval < 0 {
                break;
            }

            println!("   {:10.6e}  {:10.6e}", t, *nv_ith_s(y, 0));
            if retval >= 0 {
                tout += dTout;
                tout = if tout > Tf { Tf } else { tout };
            } else {
                println!("Solver failure, stopping integration");
                break;
            }
        }
        println!("   ---------------------\n");
    }
}
