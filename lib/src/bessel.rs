// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use rug::{
    float::{Constant, Special},
    Float,
};

use crate::quadrature;

pub fn i_n(n: &Float, z: &Float, precision: u64) -> Float {
    // based on https://dlmf.nist.gov/10.32#E3. probably subpar compared to one
    // not using numerical integration but probably good enough for now. this
    // also has the benefit of working on fractional orders for "free"

    //TODO: same comment as the one in the marcum q implementation
    let epsilon = Float::with_val_64(precision, 1e-28);

    let pi = Float::with_val_64(precision, Constant::Pi);
    let pi_reciprocal = pi.clone().recip();

    let (integrated, _) = quadrature::tanh_sinh(
        |mut theta| {
            let mut exponential = theta.clone();
            exponential.cos_mut();
            exponential *= z;
            exponential.exp_mut();

            theta *= n;
            theta.cos_mut();

            exponential * theta
        },
        (&Float::with_val_64(precision, Special::Zero), &pi),
        &epsilon,
        6,
        precision,
    );

    integrated * pi_reciprocal
}
