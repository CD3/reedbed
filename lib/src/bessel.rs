// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use std::f64;

use crate::quadrature;

pub fn i_n(n: f64, z: f64) -> f64 {
    // based on https://dlmf.nist.gov/10.32#E3. probably subpar compared to one
    // not using numerical integration but probably good enough for now. this
    // also has the benefit of working on fractional orders for "free"

    let (integrated, _) = quadrature::tanh_sinh(
        |theta| (theta.cos() * z).exp() * (theta * n).cos(),
        //TODO: same comment as the one in the marcum q implementation
        1e-9,
        (0.00, f64::consts::PI),
        6,
    );

    integrated * f64::consts::PI.recip()
}
