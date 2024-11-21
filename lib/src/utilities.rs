// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use crate::{bessel, quadrature};

//NOTE: using this function to get values for the marcum-q function with
//      very small magnitudes results in very incorrect values
pub fn marcum_q(v: i32, a: f64, b: f64) -> f64 {
    //TODO: use the variant of double-exponential quadrature supporting
    //      improper integration for this and see how it compares
    let (integrated, _) = quadrature::tanh_sinh(
        |x| {
            x.powi(v)
                * ((x.powi(2) + a.powi(2)) / -2.00).exp()
                * bessel::i_n(f64::from(v) - 1.00, x * a)
        },
        //TODO: figure out an appropriate epsilon here
        1e-9,
        (0.00, b),
        6,
    );

    1.00 - a.powi(v - 1).recip() * integrated
}
