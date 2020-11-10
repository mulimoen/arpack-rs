use std::os::raw::c_int;

#[allow(unused)]
#[derive(Copy, Clone, Debug)]
pub enum Which {
    LargestMagnitude,
    SmallestMagnitude,
    LargestRealpart,
    SmallestRealpart,
    LargestImaginarypart,
    SmallestImaginarypart,
}

impl Which {
    fn c_like(self) -> &'static [u8; 2] {
        use Which::*;
        match self {
            LargestMagnitude => b"LM",
            SmallestMagnitude => b"SM",
            LargestRealpart => b"LR",
            SmallestRealpart => b"SR",
            LargestImaginarypart => b"LI",
            SmallestImaginarypart => b"SI",
        }
    }
}

impl Default for Which {
    fn default() -> Self {
        Self::LargestMagnitude
    }
}

// mode: iparam(7) = 1 -> OP = A, B = I
//
// mode: iparam(7) = 3 -> OP = inv[A - sigma*I], B = I
//
// mode: iparam(7) = 2 -> OP = inv[M]*A, B = M
//
// shift mode for Ax = lambda Mx
// -> OP = inv[A - sigma*M]*M, B = M, sigma = real(sigma)
// iparam(7) = 3
//
// shift mode for Ax = lambda Mx
// ->OP = real(inv[A - sigma*m]*b, B = B, sigma = R + I
//
// shift invert Ax = lambda Mx
// -> OP = image(inv[A - sigma*m]*m, B = M, sigma /= real(sigma)
#[derive(Debug, Copy, Clone)]
#[allow(unused)]
enum Mode {
    /// OP = A, B = I
    ///
    /// requires Ax = y, solve for y
    Mode1,
    /// M symmetric positive definite
    /// OP = inv[M]*A, B = M
    Mode2,
    /// M symmetric semi-definite
    /// OP = real(inv[A- sigma*M]*m), B = M
    Mode3,
    /// M symmetric semi-definite
    /// OP = imag(inv[A - sigma*M]*M), B = M
    Mode4,
}

impl Mode {
    fn to_int(self) -> usize {
        match self {
            Mode::Mode1 => 1,
            Mode::Mode2 => 2,
            Mode::Mode3 => 3,
            Mode::Mode4 => 4,
        }
    }
}

#[derive(Default, Debug)]
pub struct InputParameters<'a> {
    /// Tolerance of the ritz values
    pub tolerance: f64,
    /// Initial vector of size `n`, must be non-zero
    pub initial: Option<&'a mut [f64]>,
    /// Size of the problem
    pub n: usize,
    /// Number of eigenvalues wanted
    pub nev: usize,
    /// Buffer for storing the ritz vectors
    pub v: Option<&'a mut [f64]>,
    /// Number of columns wanted/given in v
    pub ncv: usize,
    /// Which egienvectors
    pub which: Which,
    /// Maximum number of iterations
    pub mxiter: usize,
    /// Should be of size (3*n) + (3 * ncv.pow(2) + 6 * ncv)
    pub workbuffer: Option<&'a mut [f64]>,
    /// Whether to return the associated eigenvectors
    pub eigenvectors_wanted: bool,
}

/// Solving Ax = lambda x
/// using lhs :> Ax = y, solving for y
pub fn dnaupd(
    lhs: impl Fn(&[f64], &mut [f64]),
    params: InputParameters<'_>,
) -> ((Vec<f64>, Vec<f64>), Option<Vec<f64>>) {
    let tolerance = params.tolerance;
    let ncv = params.ncv;
    let ldv = params.n;
    let mut _v;
    let v = match params.v {
        Some(v) => {
            assert!(v.len() == params.ncv * params.n);
            v
        }
        None => {
            _v = vec![0.0; ldv * ncv];
            &mut _v
        }
    };

    let mut _residual;
    let mut info = params.initial.is_some() as _;
    let residual = if let Some(res) = params.initial {
        assert_eq!(res.len(), params.n);
        res
    } else {
        _residual = vec![0.0; params.n];
        &mut _residual
    };
    #[repr(C)]
    #[derive(Debug)]
    struct IParam {
        ishift: c_int,
        unused: c_int,
        mxiter: c_int,
        blocksize: c_int,
        nconv: c_int,
        iupd: c_int,
        mode: c_int,
        np: c_int,
        numop: c_int,
        numopb: c_int,
        numreo: c_int,
    }
    assert_eq!(
        std::mem::size_of::<IParam>(),
        11 * std::mem::size_of::<c_int>()
    );
    impl IParam {
        fn new(shift: bool, mxiter: usize, mode: Mode) -> Self {
            Self {
                ishift: shift as _,
                unused: 0,
                mxiter: mxiter as _,
                blocksize: 1,
                nconv: 0,
                iupd: 0,
                mode: mode.to_int() as _,
                np: if matches!(mode, Mode::Mode3) {
                    todo!()
                } else {
                    0
                },
                numop: 0,
                numopb: 0,
                numreo: 0,
            }
        }
    }

    // Only valid for mode = 2,3,4
    // to simplify the right hand side
    let bmat = 'I' as i8;
    let mut iparam = IParam::new(true, params.mxiter, Mode::Mode1);
    let mut ipntr = [0; 14];
    let mut _workd;
    let mut _workl;
    let (workd, workl) = if let Some(wb) = params.workbuffer {
        wb.split_at_mut(3 * params.n)
    } else {
        _workd = vec![0.0; 3 * params.n];
        _workl = vec![0.0; 3 * ncv.pow(2) + 6 * ncv];

        (_workd.as_mut_slice(), _workl.as_mut_slice())
    };

    let mut ido = -2;
    while ido != 99 {
        match ido {
            -2 => {
                // kickstart
                ido = 0
            }
            -1 | 1 => {
                let xstart = ipntr[0] as usize - 1;
                let ystart = ipntr[1] as usize - 1;
                let (x, y) = if ystart > xstart {
                    let (f, s) = workd.split_at_mut(ystart);
                    (&f[xstart..][..params.n], &mut s[..params.n])
                } else {
                    let (f, s) = workd.split_at_mut(xstart);
                    (&s[..params.n], &mut f[ystart..][..params.n])
                };
                lhs(x, y);
            }
            // -1 -> Y = OP * X (initialization phase)
            // 1 => mode 3 & 4 => Y = OP * x in workd(ipntr(3))
            // y = opx
            2 => {
                // Y = B*X
                todo!("ido = 3")
            }
            3 => {
                // IPARAM(8) real and imag parts of the shifts (INPTR(14) ptr
                // in workl for placing shifts)
                todo!("ido = 3")
            }
            _ => unreachable!(),
        }

        unsafe {
            arpack_sys::dnaupd_c(
                &mut ido,
                &bmat,
                params.n as _,
                params.which.c_like().as_ptr().cast(),
                params.nev as _,
                tolerance,
                residual.as_mut_ptr(),
                ncv as _,
                v.as_mut_ptr(),
                ldv as _,
                (&mut iparam as *mut IParam).cast(),
                ipntr.as_mut_ptr(),
                workd.as_mut_ptr(),
                workl.as_mut_ptr(),
                workl.len() as _,
                &mut info,
            );
        }
        assert_eq!(info, 0);
    }
    let mut dr = vec![0.0; params.nev + 1];
    let mut di = vec![0.0; params.nev + 1];
    let mut z = vec![0.0; params.n * (params.nev + 1)];
    let sigmar = 0.0;
    let sigmai = 0.0;
    let mut workev = vec![0.0; 3 * ncv];
    let rvec = params.eigenvectors_wanted;
    let howmny = b"A";

    let select = vec![false as _; ncv];
    unsafe {
        arpack_sys::dneupd_c(
            rvec,
            howmny.as_ptr().cast(),
            select.as_ptr(),
            dr.as_mut_ptr(),
            di.as_mut_ptr(),
            z.as_mut_ptr(),
            params.n as _,
            sigmar,
            sigmai,
            workev.as_mut_ptr(),
            &bmat,
            params.n as _,
            params.which.c_like().as_ptr().cast(),
            params.nev as _,
            tolerance,
            residual.as_mut_ptr(),
            ncv as _,
            v.as_mut_ptr(),
            ldv as _,
            (&mut iparam as *mut IParam).cast(),
            ipntr.as_mut_ptr(),
            workd.as_mut_ptr(),
            workl.as_mut_ptr(),
            workl.len() as _,
            &mut info,
        )
    }

    dr.pop();
    di.pop();
    let lambda = (dr, di);
    if rvec {
        z.truncate(params.n * params.nev);
        (lambda, Some(z))
    } else {
        (lambda, None)
    }
}

#[test]
fn simple_matrices() {
    let mat = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 2.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, -190.0],
        [0.0, 0.0, -190.0, 0.0],
    ];
    let n = mat.len();

    let lhs = |x: &[f64], y: &mut [f64]| {
        for j in 0..n {
            y[j] = 0.0;
        }
        for j in 0..n {
            for i in 0..n {
                y[j] += mat[j][i] * x[i]
            }
        }
    };
    let nev = 1;
    let (lambda, ev) = dnaupd(
        lhs,
        InputParameters {
            which: Which::SmallestRealpart,
            n,
            nev,
            ncv: 2 * nev + 1,
            mxiter: 10,
            ..Default::default()
        },
    );
    assert!(ev.is_none());
    approx::assert_abs_diff_eq!(&lambda.0[..], &[-190.0][..]);
    approx::assert_abs_diff_eq!(&lambda.1[..], &[0.0][..]);
    let (lambda, ev) = dnaupd(
        lhs,
        InputParameters {
            which: Which::SmallestMagnitude,
            n,
            nev,
            mxiter: 100,
            ncv: 2 * nev + 1,
            eigenvectors_wanted: true,
            ..Default::default()
        },
    );
    approx::assert_abs_diff_eq!(&lambda.0[..], &[1.0][..]);
    approx::assert_abs_diff_eq!(&lambda.1[..], &[0.0][..]);
    let ev = ev.unwrap();
    approx::assert_abs_diff_eq!(&ev[..], &[1.0, 0.0, 0.0, 0.0][..]);

    let (lambda, ev) = dnaupd(
        lhs,
        InputParameters {
            which: Which::LargestRealpart,
            n,
            nev,
            mxiter: 100,
            ncv: 2 * nev + 1,
            eigenvectors_wanted: true,
            ..Default::default()
        },
    );
    approx::assert_abs_diff_eq!(&lambda.0[..], &[190.0][..], epsilon = 1e-7);
    approx::assert_abs_diff_eq!(&lambda.1[..], &[0.0][..]);
    let ev = ev.unwrap();
    approx::assert_abs_diff_eq!(
        &ev[..],
        &[
            0.0,
            0.0,
            std::f64::consts::FRAC_1_SQRT_2,
            -std::f64::consts::FRAC_1_SQRT_2
        ][..],
        epsilon = 1e-7
    );
}
