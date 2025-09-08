#![allow(clippy::missing_safety_doc)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(clippy::approx_constant)]

use cvode_sys::*;
use std::{ffi::c_void, ptr};

pub type sunbooleantype = i32;

const NUM_SPECIES: i32 = 2;
const KH: f64 = 4.0e-6;
const VEL: f64 = 0.001;
const KV0: f64 = 1.0e-8;
const Q1: f64 = 1.63e-16;
const Q2: f64 = 4.66e-16;
const C3: f64 = 3.7e16;
const A3: f64 = 22.62;
const A4: f64 = 7.601;
const C1_SCALE: f64 = 1.0e6;
const C2_SCALE: f64 = 1.0e12;

const T0: f64 = 0.0;
const NOUT: sunindextype = 12;
const TWOHR: f64 = 7200.0;
const HALFDAY: f64 = 4.32e4;

const PI: f64 = 3.1415926535898;
// const PI: f64 = 3.1415926535897932384626433;
// use std::f64::consts::PI;

const XMIN: f64 = 0.0;
const XMAX: f64 = 20.0;
const YMIN: f64 = 30.0;
const YMAX: f64 = 50.0;
const XMID: f64 = 10.0;
const YMID: f64 = 40.0;

const MX: i32 = 10;
const MY: i32 = 10;
const NSMX: i32 = 20;
const MM: i32 = MX * MY;

const RTOL: f64 = 1.0e-5;
const FLOOR: f64 = 100.0;
const ATOL: f64 = RTOL * FLOOR;
const NEQ: i32 = NUM_SPECIES * MM;

#[inline]
pub unsafe fn ijkth(vdata: *mut sunrealtype, i: i32, j: i32, k: i32) -> *mut sunrealtype {
    unsafe {
        let idx = i - 1 + j * NUM_SPECIES + k * NSMX;
        vdata.add(idx as usize)
    }
}

#[inline]
pub unsafe fn ijth(a: *mut *mut sunrealtype, i: sunindextype, j: sunindextype) -> *mut sunrealtype {
    unsafe { (*a.add(j as usize - 1)).add(i as usize - 1) }
}

#[repr(C)]
pub struct UserData {
    pub P: [[*mut *mut sunrealtype; MY as usize]; MX as usize],
    pub Jbd: [[*mut *mut sunrealtype; MY as usize]; MX as usize],
    pub pivot: [[*mut sunindextype; MY as usize]; MX as usize],
    pub q4: sunrealtype,
    pub om: sunrealtype,
    pub dx: sunrealtype,
    pub dy: sunrealtype,
    pub hdco: sunrealtype,
    pub haco: sunrealtype,
    pub vdco: sunrealtype,
}

fn allocate_user_data() -> *mut UserData {
    unsafe {
        let mut data = Box::new(UserData {
            P: [[std::ptr::null_mut(); MY as usize]; MX as usize],
            Jbd: [[std::ptr::null_mut(); MY as usize]; MX as usize],
            pivot: [[std::ptr::null_mut(); MY as usize]; MX as usize],
            q4: 0.0,
            om: 0.0,
            dx: 0.0,
            dy: 0.0,
            hdco: 0.0,
            haco: 0.0,
            vdco: 0.0,
        });

        for jx in 0..(MX as usize) {
            for jy in 0..(MY as usize) {
                data.P[jx][jy] = SUNDlsMat_newDenseMat(NUM_SPECIES.into(), NUM_SPECIES.into());
                data.Jbd[jx][jy] = SUNDlsMat_newDenseMat(NUM_SPECIES.into(), NUM_SPECIES.into());
                data.pivot[jx][jy] = SUNDlsMat_newIndexArray(NUM_SPECIES.into());
            }
        }
        Box::into_raw(data)
    }
}

fn init_user_data(data: *mut UserData) {
    unsafe {
        (*data).om = PI / HALFDAY;
        (*data).dx = (XMAX - XMIN) / (MX - 1) as f64;
        (*data).dy = (YMAX - YMIN) / (MY - 1) as f64;
        (*data).hdco = KH / ((*data).dx * (*data).dx);
        (*data).haco = VEL / (2.0 * (*data).dx);
        (*data).vdco = (1.0 / ((*data).dy * (*data).dy)) * KV0;
    }
}

fn set_initial_profiles(u: N_Vector, dx: sunrealtype, dy: sunrealtype) {
    unsafe {
        let udata = N_VGetArrayPointer(u);
        for jy in 0..MY {
            let y = YMIN + (jy as f64) * dy;
            let mut cy = (0.1 * (y - YMID)) * (0.1 * (y - YMID));
            cy = 1.0 - cy + 0.5 * cy * cy;
            for jx in 0..MX {
                let x = XMIN + (jx as f64) * dx;
                let mut cx = (0.1 * (x - XMID)) * (0.1 * (x - XMID));
                cx = 1.0 - cx + 0.5 * cx * cx;
                *ijkth(udata, 1, jx, jy) = C1_SCALE * cx * cy;
                *ijkth(udata, 2, jx, jy) = C2_SCALE * cx * cy;
            }
        }
    }
}

extern "C" fn f(t: sunrealtype, u: N_Vector, udot: N_Vector, data: *mut c_void) -> i32 {
    unsafe {
        println!("=====================");
        let data = data as *mut UserData;
        let udata = N_VGetArrayPointer(u);
        let dudata = N_VGetArrayPointer(udot);
        // println!("val 1 = {:12.6e}", *ijkth(udata, 1, 0, 0));
        // println!("val 2 = {:12.6e}", *ijkth(udata, 2, 0, 0));

        let s = ((*data).om * t).sin();
        let (q3, q4) = if s > 0.0 {
            ((-A3 / s).exp(), (-A4 / s).exp())
        } else {
            (0.0, 0.0)
        };

        (*data).q4 = q4;

        let dely = (*data).dy;
        let verdco = (*data).vdco;
        let hordco = (*data).hdco;
        let horaco = (*data).haco;

        for jy in 0..MY {
            let ydn = YMIN + (jy as f64 - 0.5) * dely;
            let yup = ydn + dely;
            let cydn = verdco * (0.2 * ydn).exp();
            let cyup = verdco * (0.2 * yup).exp();
            let idn = if jy == 0 { 1 } else { -1 };
            let iup = if jy == MY - 1 { -1 } else { 1 };
            for jx in 0..MX {
                let c1 = *ijkth(udata, 1, jx, jy);
                let c2 = *ijkth(udata, 2, jx, jy);
                let qq1 = Q1 * c1 * C3;
                let qq2 = Q2 * c1 * c2;
                let qq3 = q3 * C3;
                let qq4 = q4 * c2;
                let rkin1 = -qq1 - qq2 + 2.0 * qq3 + qq4;
                let rkin2 = qq1 - qq2 - qq4;

                let c1dn = *ijkth(udata, 1, jx, jy + idn);
                let c2dn = *ijkth(udata, 2, jx, jy + idn);
                let c1up = *ijkth(udata, 1, jx, jy + iup);
                let c2up = *ijkth(udata, 2, jx, jy + iup);
                let vertd1 = cyup * (c1up - c1) - cydn * (c1 - c1dn);
                // println!(
                //     "cyup = {:12.10e}, c1up = {:12.10e}, c1 = {:12.10e}, cydn = {:12.10e}, c1dn = {:12.10e}",
                //     cyup, c1up, c1, cydn, c1dn
                // );
                let vertd2 = cyup * (c2up - c2) - cydn * (c2 - c2dn);

                let ileft = if jx == 0 { 1 } else { -1 };
                let iright = if jx == MX - 1 { -1 } else { 1 };
                let c1lt = *ijkth(udata, 1, jx + ileft, jy);
                let c2lt = *ijkth(udata, 2, jx + ileft, jy);
                let c1rt = *ijkth(udata, 1, jx + iright, jy);
                let c2rt = *ijkth(udata, 2, jx + iright, jy);

                let hord1 = hordco * (c1rt - 2.0 * c1 + c1lt);
                let hord2 = hordco * (c2rt - 2.0 * c2 + c2lt);

                let horad1 = horaco * (c1rt - c1lt);
                let horad2 = horaco * (c2rt - c2lt);
                // println!(
                //     "vertd1 = {:12.6e}, hord1 = {:12.6e}, horad1 = {:12.6e}, rkin1 = {:12.6e}",
                //     vertd1, hord1, horad1, rkin1
                // );

                *ijkth(dudata, 1, jx, jy) = vertd1 + hord1 + horad1 + rkin1;
                *ijkth(dudata, 2, jx, jy) = vertd2 + hord2 + horad2 + rkin2;
            }
        }
    }
    0
}

extern "C" fn jtv(
    v: N_Vector,
    Jv: N_Vector,
    t: sunrealtype,
    u: N_Vector,
    _fu: N_Vector,
    user_data: *mut c_void,
    _tmp: N_Vector,
) -> i32 {
    unsafe {
        let data = user_data as *mut UserData;
        let udata = N_VGetArrayPointer(u);
        let vdata = N_VGetArrayPointer(v);
        let Jvdata = N_VGetArrayPointer(Jv);

        let s = ((*data).om * t).sin();
        let q4 = if s > 0.0 { (-A4 / s).exp() } else { 0.0 };

        (*data).q4 = q4;

        let dely = (*data).dy;
        let verdco = (*data).vdco;
        let hordco = (*data).hdco;
        let horaco = (*data).haco;

        for jy in 0..MY {
            let ydn = YMIN + (jy as f64 - 0.5) * dely;
            let yup = ydn + dely;

            let cydn = verdco * (0.2 * ydn).exp();
            let cyup = verdco * (0.2 * yup).exp();
            let idn = if jy == 0 { 1 } else { -1 };
            let iup = if jy == MY - 1 { -1 } else { 1 };
            for jx in 0..MX {
                let mut Jv1 = 0.0;
                let mut Jv2 = 0.0;
                let c1 = *ijkth(udata, 1, jx, jy);
                let c2 = *ijkth(udata, 2, jx, jy);

                let v1 = *ijkth(vdata, 1, jx, jy);
                let v2 = *ijkth(vdata, 2, jx, jy);

                let v1dn = *ijkth(vdata, 1, jx, jy + idn);
                let v2dn = *ijkth(vdata, 2, jx, jy + idn);
                let v1up = *ijkth(vdata, 1, jx, jy + iup);
                let v2up = *ijkth(vdata, 2, jx, jy + iup);

                let ileft = if jx == 0 { 1 } else { -1 };
                let iright = if jx == MX - 1 { -1 } else { 1 };

                let v1lt = *ijkth(vdata, 1, jx + ileft, jy);
                let v2lt = *ijkth(vdata, 2, jx + ileft, jy);
                let v1rt = *ijkth(vdata, 1, jx + iright, jy);
                let v2rt = *ijkth(vdata, 2, jx + iright, jy);

                Jv1 += -(Q1 * C3 + Q2 * c2) * v1 + (q4 - Q2 * c1) * v2;
                Jv2 += (Q1 * C3 - Q2 * c2) * v1 - (q4 + Q2 * c1) * v2;

                Jv1 += -(cyup + cydn) * v1 + cyup * v1up + cydn * v1dn;
                Jv2 += -(cyup + cydn) * v2 + cyup * v2up + cydn * v2dn;

                Jv1 += hordco * (v1rt - 2.0 * v1 + v1lt);
                Jv2 += hordco * (v2rt - 2.0 * v2 + v2lt);

                Jv1 += horaco * (v1rt - v1lt);
                Jv2 += horaco * (v2rt - v2lt);

                *ijkth(Jvdata, 1, jx, jy) = Jv1;
                *ijkth(Jvdata, 2, jx, jy) = Jv2;
            }
        }
    }
    0
}

extern "C" fn precond(
    _tn: sunrealtype,
    u: N_Vector,
    _fu: N_Vector,
    jok: sunbooleantype,
    jcurPtr: *mut sunbooleantype,
    gamma: sunrealtype,
    user_data: *mut c_void,
) -> i32 {
    unsafe {
        let data = user_data as *mut UserData;
        let P = (*data).P;
        let Jbd = (*data).Jbd;
        let pivot = (*data).pivot;
        let udata = N_VGetArrayPointer(u);

        if jok == 1 {
            for jy in 0..(MY as usize) {
                for jx in 0..(MX as usize) {
                    SUNDlsMat_denseCopy(
                        Jbd[jx][jy],
                        P[jx][jy],
                        NUM_SPECIES.into(),
                        NUM_SPECIES.into(),
                    );
                }
            }
            *jcurPtr = 0;
        } else {
            let q4coef = (*data).q4;
            let dely = (*data).dy;
            let verdco = (*data).vdco;
            let hordco = (*data).hdco;

            for jy in 0..MY {
                let ydn = YMIN + (jy as f64 - 0.5) * dely;
                let yup = ydn + dely;
                let cydn = verdco * (0.2 * ydn).exp();
                let cyup = verdco * (0.2 * yup).exp();
                let diag = -(cydn + cyup + 2.0 * hordco);
                for jx in 0..MX {
                    let c1 = *ijkth(udata, 1, jx, jy);
                    let c2 = *ijkth(udata, 2, jx, jy);

                    let j = Jbd[jx as usize][jy as usize];
                    let a = P[jx as usize][jy as usize];
                    *ijth(j, 1, 1) = (-Q1 * C3 - Q2 * c2) + diag;
                    *ijth(j, 1, 2) = -Q2 * c1 + q4coef;
                    *ijth(j, 2, 1) = Q1 * C3 - Q2 * c2;
                    *ijth(j, 2, 2) = (-Q2 * c1 - q4coef) + diag;
                    SUNDlsMat_denseCopy(j, a, NUM_SPECIES.into(), NUM_SPECIES.into());
                }
            }
            *jcurPtr = 1;
        }

        for jy in 0..(MY as usize) {
            for jx in 0..(MX as usize) {
                SUNDlsMat_denseScale(-gamma, P[jx][jy], NUM_SPECIES.into(), NUM_SPECIES.into());
            }
        }

        for jx in 0..(MX as usize) {
            for jy in 0..(MY as usize) {
                SUNDlsMat_denseAddIdentity(P[jx][jy], NUM_SPECIES.into());
                let retval = SUNDlsMat_denseGETRF(
                    P[jx][jy],
                    NUM_SPECIES.into(),
                    NUM_SPECIES.into(),
                    pivot[jx][jy],
                );
                assert!(retval == 0);
            }
        }
    }
    0
}

extern "C" fn psolve(
    _tn: sunrealtype,
    _u: N_Vector,
    _fu: N_Vector,
    r: N_Vector,
    z: N_Vector,
    _gamma: sunrealtype,
    _delta: f64,
    _lr: i32,
    user_data: *mut c_void,
) -> i32 {
    unsafe {
        let data = user_data as *mut UserData;
        let P = (*data).P;
        let pivot = (*data).pivot;
        let zdata = N_VGetArrayPointer(z);

        N_VScale(1.0, r, z);

        for jx in 0..MX {
            for jy in 0..MY {
                let v = &mut *ijkth(zdata, 1, jx, jy);
                SUNDlsMat_denseGETRS(
                    P[jx as usize][jy as usize],
                    NUM_SPECIES.into(),
                    pivot[jx as usize][jy as usize],
                    v,
                );
            }
        }
    }
    0
}

pub unsafe fn print_output(cvode_mem: *mut c_void, u: N_Vector, t: sunrealtype) {
    unsafe {
        let mxh = MX / 2 - 1;
        let myh = MY / 2 - 1;
        let mx1 = MX - 1;
        let my1 = MY - 1;

        let mut nst = 0;
        let mut qu = 0;
        let mut hu = 0.0;

        let udata = N_VGetArrayPointer(u);
        let mut retval = CVodeGetNumSteps(cvode_mem, &mut nst);
        assert!(retval == 0);

        retval = CVodeGetLastOrder(cvode_mem, &mut qu);
        assert!(retval == 0);

        retval = CVodeGetLastStep(cvode_mem, &mut hu);
        assert!(retval == 0);

        println!(
            "t = {:.2e}   no. steps = {}   order = {}   stepsize = {:.2e}",
            t, nst, qu, hu
        );

        println!(
            "c1 (bot.left/middle/top rt.) = {:12.3e}  {:12.3e}  {:12.3e}",
            *ijkth(udata, 1, 0, 0),
            *ijkth(udata, 1, mxh, myh),
            *ijkth(udata, 1, mx1, my1),
        );

        println!(
            "c2 (bot.left/middle/top rt.) = {:12.3e}  {:12.3e}  {:12.3e}\n",
            *ijkth(udata, 2, 0, 0),
            *ijkth(udata, 2, mxh, myh),
            *ijkth(udata, 2, mx1, my1),
        );
    }
}

fn main() {
    unsafe {
        let mut sunctx: SUNContext = ptr::null_mut();
        let mut retval = SUNContext_Create(SUN_COMM_NULL.try_into().unwrap(), &mut sunctx);
        assert!(retval == 0);

        let u = N_VNew_Serial(NEQ.into(), sunctx);
        assert!(!u.is_null());

        let data = allocate_user_data();
        assert!(!data.is_null());

        init_user_data(data);
        set_initial_profiles(u, (*data).dx, (*data).dy);

        let abstol = ATOL;
        let reltol = RTOL;

        let cvode_mem = CVodeCreate(CV_BDF.try_into().unwrap(), sunctx);
        assert!(!cvode_mem.is_null());

        retval = CVodeSetUserData(cvode_mem, data as *mut c_void);
        assert!(retval == 0);

        retval = CVodeInit(cvode_mem, Some(f), T0, u);
        assert!(retval == 0);

        retval = CVodeSStolerances(cvode_mem, reltol, abstol);
        assert!(retval == 0);

        let LS = SUNLinSol_SPGMR(u, SUN_PREC_LEFT.try_into().unwrap(), 0, sunctx);
        assert!(!LS.is_null());

        retval = CVodeSetLinearSolver(cvode_mem, LS, ptr::null_mut());
        assert!(retval == 0);

        retval = CVodeSetJacTimes(cvode_mem, None, Some(jtv));
        assert!(retval == 0);

        retval = CVodeSetPreconditioner(cvode_mem, Some(precond), Some(psolve));
        assert!(retval == 0);

        println!("\n2-species diurnal advection-diffusion problem\n\n");

        let mut t: f64 = T0;
        let mut tout = TWOHR;
        for iout in 1..=NOUT {
            println!("iout = {}", iout);
            let retval = CVode(cvode_mem, tout, u, &mut t, CV_NORMAL.try_into().unwrap());
            print_output(cvode_mem, u, t);
            if retval < 0 {
                break;
            }
            tout += TWOHR;
        }
        N_VDestroy(u);
        CVodeFree(&cvode_mem as *const _ as *mut _);
        SUNLinSolFree(LS);
        SUNContext_Free(&mut sunctx);
    }
}
