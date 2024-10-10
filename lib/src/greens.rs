// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use rug::{float::Special, Float};
use std::borrow::Cow;

use crate::utilities;

//TODO: we could probably swap the use of [`struct@Float`] for a generic
//      parameter that implements the operation traits in rug::ops in most
//      (if not all) places

/// An abstraction over the various `*AbsorbingLayer` structures
pub trait AbsorbingLayer<'a>: ToOwned + Sized {
    /// Convert this absorbing layer into an owned one
    ///
    /// This is effectively a no-op if all data was already owned
    fn into_owned(self) -> <Self as ToOwned>::Owned;

    /// Get the location at which this layer starts
    ///
    /// This is usually the `z0` parameter on the individual structures
    fn starts_at(&'a self) -> &'a Float;

    /// Evaluate the absorbing layer using the specified parameters
    ///
    /// Not all implementations of [`trait@AbsorbingLayer`] will use all
    /// parameters
    fn evaluate_at(&self, z: &Float, r: &Float, tp: &Float) -> Float;
}

//TODO: document these parameters better
#[derive(Clone, PartialEq)]
pub struct LargeBeamAbsorbingLayer<'a> {
    /// Units: cm^-1
    pub mu_a: Cow<'a, Float>,

    /// Units: g*cm^3
    pub rho: Cow<'a, Float>,

    /// Units: J*g^-1*K^-1
    pub c: Cow<'a, Float>,

    /// Units: W*cm^-1*K^-1
    pub k: Cow<'a, Float>,

    /// Units: cm
    pub d: Cow<'a, Float>,

    /// Units: cm
    pub z0: Cow<'a, Float>,

    /// Units: W*cm^-2
    pub e0: Cow<'a, Float>,

    /// Precision (in bits) of the arbitrary-precision floating point values
    /// used in intermediate calculations (and in the output)
    pub precision: u64,
}

impl<'a> AbsorbingLayer<'a> for LargeBeamAbsorbingLayer<'a> {
    fn into_owned(self) -> <Self as ToOwned>::Owned {
        Self {
            mu_a: Cow::Owned(self.mu_a.into_owned()),
            rho: Cow::Owned(self.rho.into_owned()),
            c: Cow::Owned(self.c.into_owned()),
            k: Cow::Owned(self.k.into_owned()),
            d: Cow::Owned(self.d.into_owned()),
            z0: Cow::Owned(self.z0.into_owned()),
            e0: Cow::Owned(self.e0.into_owned()),
            precision: self.precision,
        }
    }

    fn starts_at(&'a self) -> &'a Float {
        &self.z0
    }

    fn evaluate_at(&self, z: &Float, _r: &Float, tp: &Float) -> Float {
        let LargeBeamAbsorbingLayer {
            mu_a,
            rho,
            c,
            k,
            d,
            z0,
            e0,
            precision,
        } = self;
        let precision = *precision;

        //TODO: make this less naive

        let mut alpha = Float::with_val_64(precision, k.as_ref());
        alpha /= rho.as_ref();
        alpha /= c.as_ref();

        let mut term_1 = Float::with_val_64(precision, mu_a.as_ref());
        term_1 *= e0.as_ref();
        term_1 /= rho.as_ref();
        term_1 /= c.as_ref();
        term_1 /= 2.0;

        let mut term_2 = Float::with_val_64(precision, z);
        term_2 -= z0.as_ref();
        term_2 *= mu_a.as_ref();
        term_2 *= -1;
        term_2.exp_mut();

        if *tp == 0 {
            return term_1 * term_2;
        }

        let mut term_3 = Float::with_val_64(precision, mu_a.as_ref());
        term_3.square_mut();
        term_3 *= tp;
        term_3 *= &alpha;
        term_3.exp_mut();

        let mut reciprocal_sqrt = Float::with_val_64(precision, &alpha);
        reciprocal_sqrt *= tp;
        reciprocal_sqrt *= 4.0;
        reciprocal_sqrt.sqrt_mut();
        reciprocal_sqrt.recip_mut();

        let mut sqrt_mu_a = alpha;
        sqrt_mu_a *= tp;
        sqrt_mu_a.sqrt_mut();
        sqrt_mu_a *= mu_a.as_ref();

        let mut argument_1 = Float::with_val_64(precision, z0.as_ref());
        argument_1 += d.as_ref();
        argument_1 -= z;
        argument_1 *= &reciprocal_sqrt;
        argument_1 += &sqrt_mu_a;
        argument_1.erf_mut();

        let mut argument_2 = Float::with_val_64(precision, z0.as_ref());
        argument_2 -= z;
        argument_2 *= &reciprocal_sqrt;
        argument_2 += &sqrt_mu_a;
        argument_2.erf_mut();

        let mut term_4 = argument_1;
        term_4 -= argument_2;

        term_1 * term_2 * term_3 * term_4
    }
}

//TODO: same todo as above
#[derive(Clone, PartialEq)]
pub struct FlatTopBeamAbsorbingLayer<'a> {
    /// Units: cm^-1
    pub mu_a: Cow<'a, Float>,

    /// Units: g*cm^3
    pub rho: Cow<'a, Float>,

    /// Units: J*g^-1*K^-1
    pub c: Cow<'a, Float>,

    /// Units: W*cm^-1*K^-1
    pub k: Cow<'a, Float>,

    /// Units: cm
    pub d: Cow<'a, Float>,

    /// Units: cm
    pub z0: Cow<'a, Float>,

    /// Units: W*cm^-2
    pub e0: Cow<'a, Float>,

    /// Units: cm
    pub radius: Cow<'a, Float>,

    /// Precision (in bits) of the arbitrary-precision floating point values
    /// used in intermediate calculations (and in the output)
    pub precision: u64,
}

impl<'a> AbsorbingLayer<'a> for FlatTopBeamAbsorbingLayer<'a> {
    fn into_owned(self) -> <Self as ToOwned>::Owned {
        Self {
            mu_a: Cow::Owned(self.mu_a.into_owned()),
            rho: Cow::Owned(self.rho.into_owned()),
            c: Cow::Owned(self.c.into_owned()),
            k: Cow::Owned(self.k.into_owned()),
            d: Cow::Owned(self.d.into_owned()),
            z0: Cow::Owned(self.z0.into_owned()),
            e0: Cow::Owned(self.e0.into_owned()),
            radius: Cow::Owned(self.radius.into_owned()),
            precision: self.precision,
        }
    }

    fn starts_at(&'a self) -> &'a Float {
        &self.z0
    }

    fn evaluate_at(&self, z: &Float, r: &Float, tp: &Float) -> Float {
        let FlatTopBeamAbsorbingLayer {
            mu_a,
            rho,
            c,
            k,
            d,
            z0,
            e0,
            radius,
            precision,
        } = self;
        let precision = *precision;

        if *tp == 0 && r > radius.as_ref() {
            return Float::with_val_64(precision, Special::Zero);
        }

        let z_factor = LargeBeamAbsorbingLayer {
            mu_a: Cow::Borrowed(mu_a.as_ref()),
            rho: Cow::Borrowed(rho.as_ref()),
            c: Cow::Borrowed(c.as_ref()),
            k: Cow::Borrowed(k.as_ref()),
            d: Cow::Borrowed(d.as_ref()),
            z0: Cow::Borrowed(z0.as_ref()),
            e0: Cow::Borrowed(e0.as_ref()),
            precision,
        }
        .evaluate_at(z, &Float::with_val_64(precision, Special::Zero), tp);

        if *tp == 0 {
            return z_factor;
        }

        //TODO: don't duplicate this between large_beam_absorbing_layer and
        //      this function
        let mut alpha = Float::with_val_64(precision, k.as_ref());
        alpha /= rho.as_ref();
        alpha /= c.as_ref();

        z_factor
            * if *r == 0 {
                let mut r_factor = Float::with_val_64(precision, radius.as_ref());
                r_factor.square_mut();
                r_factor /= -4.0;
                r_factor /= alpha;
                r_factor /= tp;
                r_factor.exp_mut();
                r_factor = 1 - r_factor;
                r_factor
            } else {
                //TODO: this is not accurate at all. fix the marcum-q function
                //      implementation

                let mut a = Float::with_val_64(precision, 2.0);
                a *= alpha;
                a *= tp;
                a.recip_mut();

                let mut b = a.clone();
                b *= radius.as_ref();
                a *= r;

                let mut r_factor = utilities::marcum_q(1, &a, &b, precision);
                r_factor = 1 - r_factor;
                r_factor
            }
    }
}

/*pub struct MultiAbsorbingLayer<'a> {
    layers: Vec<Box<dyn AbsorbingLayer<'static>>>,
}*/

#[cfg(test)]
mod tests {
    use super::*;

    #[ctor::ctor]
    static ZERO: Float = Float::with_val_64(64, Special::Zero);

    #[ctor::ctor]
    static ONE: Float = Float::with_val_64(64, 1.0);

    #[ctor::ctor]
    static EPSILON: Float = Float::with_val_64(64, 1e-16);

    #[test]
    fn large_beam_sanity() {
        let layer = LargeBeamAbsorbingLayer {
            mu_a: Cow::Borrowed(&ONE),
            rho: Cow::Borrowed(&ONE),
            c: Cow::Borrowed(&ONE),
            k: Cow::Borrowed(&ONE),
            d: Cow::Borrowed(&ONE),
            z0: Cow::Borrowed(&ZERO),
            e0: Cow::Borrowed(&ONE),
            precision: 64,
        };

        assert_eq!(
            layer.evaluate_at(&ZERO, &ZERO, &ZERO),
            5e-1
        );

        let mut result = layer.evaluate_at(&ONE, &ZERO, &ZERO);
        // reference result: 0.5 * e^-1
        result -= 1.8393972058572116080e-1;
        result.abs_mut();
        assert!(result < *EPSILON);

        let mut result = layer.evaluate_at(&ONE, &ZERO, &ONE);
        // reference result: 0.5 * e^-1 * e^1 * (erf(1) - erf(-1/sqrt(4) + 1))
        result -= 1.6110045756833416583e-1;
        result.abs_mut();
        println!("{}", result);
        assert!(result < *EPSILON);
    }

    #[test]
    fn flat_top_beam_sanity() {
        let layer = FlatTopBeamAbsorbingLayer {
            mu_a: Cow::Borrowed(&ONE),
            rho: Cow::Borrowed(&ONE),
            c: Cow::Borrowed(&ONE),
            k: Cow::Borrowed(&ONE),
            d: Cow::Borrowed(&ONE),
            z0: Cow::Borrowed(&ZERO),
            e0: Cow::Borrowed(&ONE),
            radius: Cow::Borrowed(&ONE),
            precision: 64,
        };

        assert_eq!(
            layer.evaluate_at(&ZERO, &ZERO, &ZERO),
            5e-1
        );

        let mut result = layer.evaluate_at(&ONE, &ZERO, &ZERO);
        // reference result: 0.5 * e^-1 * (1 - 0)
        result -= 1.8393972058572116080e-1;
        result.abs_mut();
        assert!(result < *EPSILON);

        let mut result = layer.evaluate_at(&ONE, &ZERO, &ONE);
        // reference result: 0.5 * e^-1 * e^1 * (erf(1) - erf(-1/sqrt(4) + 1))
        //                       * (1 - e^(-1/4))
        result -= 3.5635295060953884529e-2;
        result.abs_mut();
        println!("{}", result);
        assert!(result < *EPSILON);
    }
}
