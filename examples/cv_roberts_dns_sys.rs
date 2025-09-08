#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

use libc::{FILE, fopen};
use std::ffi::{CString, c_int, c_void};

use std::ptr;

use cvode_sys::*;

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

/// Get the i-th element (one-based) from the serial vector.
///
/// # Safety
/// - `v` must be a valid non-null pointer to an `N_Vector`.
/// - Index `i` must be in the range `1..=length`.
/// - Converts to zero-based indexing internally.
/// - Caller must ensure no out-of-bounds access.
#[inline]
pub unsafe fn ith(v: N_Vector, i: usize) -> *mut sunrealtype {
    unsafe { nv_ith_s(v, i - 1) }
}

/// Returns the dense matrix content of a SUNMatrix.
///
/// # Safety
/// - `A` must be a valid, non-null `SUNMatrix`.
/// - `A` must actually contain `SUNMatrixContent_Dense`.
#[inline]
pub unsafe fn sm_content_d(A: SUNMatrix) -> SUNMatrixContent_Dense {
    unsafe { (*A).content as SUNMatrixContent_Dense }
}

/// Returns the number of rows in a dense SUNMatrix.
///
/// # Safety
/// - `A` must be a valid dense SUNMatrix pointer.
#[inline]
pub unsafe fn sm_rows_d(A: SUNMatrix) -> sunindextype {
    unsafe { (*sm_content_d(A)).M }
}

/// Returns the number of columns in a dense SUNMatrix.
///
/// # Safety
/// - `A` must be a valid dense SUNMatrix pointer.
#[inline]
pub unsafe fn sm_columns_d(A: SUNMatrix) -> sunindextype {
    unsafe { (*sm_content_d(A)).N }
}

/// Returns the leading dimension (`ldata`) of a dense SUNMatrix.
///
/// # Safety
/// - `A` must be a valid dense SUNMatrix pointer.
#[inline]
pub unsafe fn sm_ldata_d(A: SUNMatrix) -> sunindextype {
    unsafe { (*sm_content_d(A)).ldata }
}

/// Returns a pointer to the contiguous data array of a dense SUNMatrix.
///
/// # Safety
/// - `A` must be a valid dense SUNMatrix pointer.
/// - The returned pointer must not be used after `A` is freed.
#[inline]
pub unsafe fn sm_data_d(A: SUNMatrix) -> *mut sunrealtype {
    unsafe { (*sm_content_d(A)).data }
}

/// Returns a pointer to the array of column pointers (`cols`) of a dense SUNMatrix.
///
/// # Safety
/// - `A` must be a valid dense SUNMatrix pointer.
/// - The returned pointer must not be dereferenced beyond the number of columns (`N`).
#[inline]
pub unsafe fn sm_cols_d(A: SUNMatrix) -> *mut *mut sunrealtype {
    unsafe { (*sm_content_d(A)).cols }
}

/// Returns a pointer to the `j`-th column of a dense SUNMatrix.
///
/// # Safety
/// - `A` must be a valid dense SUNMatrix pointer.
/// - `j` must be less than the number of columns (`N`) in the matrix.
#[inline]
pub unsafe fn sm_column_d(A: SUNMatrix, j: usize) -> *mut sunrealtype {
    unsafe { *(*sm_content_d(A)).cols.add(j) }
}

/// Returns a pointer to the element at row `i` and column `j` of a dense SUNMatrix.
///
/// # Safety
/// - `A` must be a valid dense SUNMatrix pointer.
/// - `i` must be less than the number of rows (`M`).
/// - `j` must be less than the number of columns (`N`).
#[inline]
pub unsafe fn sm_element_d(A: SUNMatrix, i: usize, j: usize) -> *mut sunrealtype {
    unsafe { sm_column_d(A, j).add(i) }
}

/// Get the (i, j)-th element of a dense SUNMatrix with
/// one-based indexing.
///
/// # Safety
/// - `A` must be a valid dense SUNMatrix pointer.
/// - `i` must be less than or equal to the number of rows (`M`).
/// - `j` must be less than or equal to the number of columns (`N`).
#[inline]
pub unsafe fn ijth(A: SUNMatrix, i: usize, j: usize) -> *mut sunrealtype {
    unsafe { sm_element_d(A, i - 1, j - 1) }
}

pub const NEQ: sunindextype = 3; // number of equations
pub const Y1: f64 = 1.0; // initial y components
pub const Y2: f64 = 0.0;
pub const Y3: f64 = 0.0;
pub const RTOL: f64 = 1.0e-8; // scalar relative tolerance
pub const ATOL1: f64 = 1.0e-14; // vector absolute tolerance components
pub const ATOL2: f64 = 1.0e-14;
pub const ATOL3: f64 = 1.0e-14;
pub const T0: f64 = 0.0; // initial time
pub const T1: f64 = 0.4; // first output time
pub const TMULT: f64 = 10.0; // output time factor
pub const NOUT: sunindextype = 12; // number of output times

pub const ZERO: f64 = 0.0;

unsafe fn print_output(t: f64, y1: f64, y2: f64, y3: f64) {
    println!(
        "At t = {:0.4e}      y ={:14.6e}  {:14.6e}  {:14.6e}",
        t, y1, y2, y3
    );
}

unsafe fn print_root_info(root_f1: i32, root_f2: i32) {
    println!("    rootsfound[] = {:3} {:3}", root_f1, root_f2);
}

extern "C" fn f(_t: sunrealtype, y: N_Vector, ydot: N_Vector, _user_data: *mut c_void) -> i32 {
    unsafe {
        let y1 = *ith(y, 1);
        let y2 = *ith(y, 2);
        let y3 = *ith(y, 3);

        // Compute derivatives
        let yd1 = -0.04 * y1 + 1.0e4 * y2 * y3;
        // println!("{:.4e} {:.6e} {:.6e} {:.6e} {:.6e}", t, y1, y2, y3, yd1);
        *ith(ydot, 1) = yd1;

        let yd3 = 3.0e7 * y2 * y2;
        *ith(ydot, 3) = yd3;

        *ith(ydot, 2) = -yd1 - yd3;
        // println!("ydot: [{:.5e}, {:.5e}, {:.5e}]", yd1, -yd1 - yd3, yd3);
    }
    0
}

extern "C" fn jac(
    _t: sunrealtype,
    y: N_Vector,
    _fy: N_Vector,
    J: SUNMatrix,
    _user_data: *mut c_void,
    _tmp1: N_Vector,
    _tmp2: N_Vector,
    _tmp3: N_Vector,
) -> i32 {
    unsafe {
        // println!("In jacobian");
        let y2 = *ith(y, 2);
        let y3 = *ith(y, 3);

        *ijth(J, 1, 1) = -0.04;
        // println!("J[1,1] = {}", *ijth(J, 1, 1));
        *ijth(J, 1, 2) = 1.0e4 * y3;
        *ijth(J, 1, 3) = 1.0e4 * y2;
        *ijth(J, 2, 1) = 0.04;
        *ijth(J, 2, 2) = -1.0e4 * y3 - 6.0e7 * y2;
        *ijth(J, 2, 3) = -1.0e4 * y2;
        *ijth(J, 3, 1) = 0.0;
        *ijth(J, 3, 2) = 6.0e7 * y2;
        *ijth(J, 3, 3) = 0.0;
    }
    0
}

extern "C" fn g(
    _t: sunrealtype,
    y: N_Vector,
    gout: *mut sunrealtype,
    _user_data: *mut c_void,
) -> i32 {
    unsafe {
        let y1 = *ith(y, 1);
        let y3 = *ith(y, 3);

        *gout.add(0) = y1 - 0.0001;
        *gout.add(1) = y3 - 0.01;
    }

    0
}

fn main() {
    unsafe {
        let mut sunctx: SUNContext = ptr::null_mut();
        let mut retval = SUNContext_Create(SUN_COMM_NULL.try_into().unwrap(), &mut sunctx);
        assert!(retval == 0);

        let y = N_VNew_Serial(NEQ, sunctx);
        assert!(!y.is_null());

        // Initialize
        *ith(y, 1) = Y1;
        *ith(y, 2) = Y2;
        *ith(y, 3) = Y3;

        let abstol = N_VNew_Serial(NEQ, sunctx);
        assert!(!abstol.is_null());

        *ith(abstol, 1) = ATOL1;
        *ith(abstol, 2) = ATOL2;
        *ith(abstol, 3) = ATOL3;

        let cvode_mem = CVodeCreate(CV_BDF.try_into().unwrap(), sunctx);

        assert!(!cvode_mem.is_null());

        retval = CVodeInit(cvode_mem, Some(f), T0, y);
        assert!(retval == 0);

        retval = CVodeSVtolerances(cvode_mem, RTOL, abstol);
        assert!(retval == 0);

        retval = CVodeRootInit(cvode_mem, 2, Some(g));
        assert!(retval == 0);

        let A = SUNDenseMatrix(NEQ, NEQ, sunctx);
        assert!(!A.is_null());

        let LS = SUNLinSol_Dense(y, A, sunctx);
        assert!(!LS.is_null());

        retval = CVodeSetLinearSolver(cvode_mem, LS, A);
        assert!(retval == 0);

        retval = CVodeSetJacFn(cvode_mem, Some(jac));
        assert!(retval == 0);

        println!("\n3-species kinetics problem\n\n");

        let mut iout = 0;
        let mut t: f64 = 0.0;
        let mut tout = T1;
        let mut rootsfound: Vec<c_int> = vec![0; 2];

        let path = CString::new(
            "/Users/kmaitreys/Documents/probe/sun-examples/examples/cv_roberts_dns_out.csv",
        )
        .unwrap();
        let mode = CString::new("w").unwrap();
        let fid: *mut FILE = fopen(path.as_ptr(), mode.as_ptr());
        assert!(!fid.is_null());
        let fid_sys: *mut __sFILE = fid.cast();

        loop {
            retval = CVode(cvode_mem, tout, y, &mut t, CV_NORMAL.try_into().unwrap());
            // println!("CVode order: current={}, last={}", order, last_order);
            print_output(t, *ith(y, 1), *ith(y, 2), *ith(y, 3));

            if retval == CV_ROOT_RETURN.try_into().unwrap() {
                let retvalr = CVodeGetRootInfo(cvode_mem, rootsfound.as_mut_ptr());
                assert!(retvalr == 0);
                print_root_info(rootsfound[0], rootsfound[1]);
            }

            if retval < 0 {
                // CVode error codes are negative
                break;
            }

            if retval == CV_SUCCESS.try_into().unwrap() {
                iout += 1;
                tout *= TMULT;
            }

            retval = CVodePrintAllStats(cvode_mem, fid_sys, SUNOutputFormat_SUN_OUTPUTFORMAT_CSV);
            assert!(retval == 0);

            if iout == NOUT {
                break;
            }
        }
        fclose(fid_sys);

        println!("\nFinal statistics:\n");
        // retval = CVodePrintAllStats(cvode_mem, stdout, SUNOutputFormat_SUN_OUTPUTFORMAT_TABLE);
        // assert!(retval == 0);

        N_VDestroy(y);
        N_VDestroy(abstol);
        CVodeFree(&cvode_mem as *const _ as *mut _);
        SUNLinSolFree(LS);
        SUNMatDestroy(A);
        SUNContext_Free(&mut sunctx);
    }
}
