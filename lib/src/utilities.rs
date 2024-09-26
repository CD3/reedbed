// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use rug::{
    float::Special,
    ops::{CompleteRound, Pow},
    Float,
};

use crate::{bessel, quadrature};

pub fn marcum_q(v: i32, a: &Float, b: &Float, precision: u64) -> Float {
    //TODO: use the variant of double-exponential quadrature supporting
    //      improper integration for this and see how it compares

    //NOTE: using this function to get values for the marcum-q function with
    //      very small magnitudes results in very incorrect values when the
    //      precision is insufficient

    //TODO: figure out an appropriate epsilon here
    //TODO: this should probably be dynamically scaled based on the input
    //      precision
    let epsilon = Float::with_val_64(precision, 1e-18);

    let (integrated, _) = quadrature::tanh_sinh(
        |x| {
            let mut two = x.clone();

            two.square_mut();
            two += a.clone().square();
            two /= -2.0;
            two.exp_mut();

            let mut three = x.clone();

            three *= a;

            three = bessel::i_n(&Float::with_val_64(precision, v - 1), &three, precision);

            x.pow(v) * two * three
        },
        (&Float::with_val_64(precision, Special::Zero), &b),
        &epsilon,
        6,
        precision,
    );

    let mut argument = Float::new_64(precision);
    a.pow(v - 1).complete_into(&mut argument);

    argument.recip_mut();

    1 - argument * integrated
}
