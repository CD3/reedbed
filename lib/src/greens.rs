// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use rug::{float::Special, Float};
use std::borrow::Cow;

use crate::utilities;

/// A configuration structure for specific thermal properties
#[derive(Clone, PartialEq)]
pub struct ThermalProperties<'a> {
    /// Units: g*cm^3
    pub rho: Cow<'a, Float>,

    /// Units: J*g^-1*K^-1
    pub c: Cow<'a, Float>,

    /// Units: W*cm^-1*K^-1
    pub k: Cow<'a, Float>,
}

/// A layer of tissue
#[derive(Clone, PartialEq)]
pub struct Layer<'a> {
    /// Units: cm
    pub d: Cow<'a, Float>,

    /// Units: cm
    pub z0: Cow<'a, Float>,

    /// Units: cm^-1
    pub mu_a: Cow<'a, Float>,

    /// Units: W*cm^-2
    pub e0: Cow<'a, Float>,
}

//TODO: we could probably swap the use of [`struct@Float`] for a generic
//      parameter that implements the operation traits in rug::ops in most
//      (if not all) places

/// An abstraction over the various `*Beam` structures
pub trait Beam: ToOwned + Sized {
    /// Convert this [`trait@Beam`] into an owned one
    ///
    /// This is effectively a no-op if all data was already owned
    fn into_owned(self) -> <Self as ToOwned>::Owned;

    /// Run the beam over a given [`struct@Layer`] with the provided
    /// [`struct@ThermalProperties`]
    ///
    /// Not all implementations of [`trait@Beam`] will use all
    /// parameters
    fn evaluate_at<'a>(
        &self,
        precision: u64,
        thermal_properties: &ThermalProperties<'a>,
        layer: &Layer<'a>,
        z: &Float,
        r: &Float,
        tp: &Float,
    ) -> Float;
}

#[derive(Clone, PartialEq)]
pub struct LargeBeam;

impl Beam for LargeBeam {
    fn into_owned(self) -> <Self as ToOwned>::Owned {
        Self
    }

    //TODO: it (might?) be worthwhile to have a specialized method that
    //      doesn't need to take r. however, this could also be addressed with
    //      the genericization of this method at the trait level. see above
    //      for more details
    fn evaluate_at<'a>(
        &self,
        precision: u64,
        thermal_properties: &ThermalProperties<'a>,
        layer: &Layer<'a>,
        z: &Float,
        _r: &Float,
        tp: &Float,
    ) -> Float {
        //TODO: make this less naive

        let mut alpha = Float::with_val_64(precision, thermal_properties.k.as_ref());
        alpha /= thermal_properties.rho.as_ref();
        alpha /= thermal_properties.c.as_ref();

        let mut term_1 = Float::with_val_64(precision, layer.mu_a.as_ref());
        term_1 *= layer.e0.as_ref();
        term_1 /= thermal_properties.rho.as_ref();
        term_1 /= thermal_properties.c.as_ref();
        term_1 /= 2.0;

        let mut term_2 = Float::with_val_64(precision, z);
        term_2 -= layer.z0.as_ref();
        term_2 *= layer.mu_a.as_ref();
        term_2 *= -1;
        term_2.exp_mut();

        if *tp == 0 {
            return term_1 * term_2;
        }

        let mut term_3 = Float::with_val_64(precision, layer.mu_a.as_ref());
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
        sqrt_mu_a *= layer.mu_a.as_ref();

        let mut argument_1 = Float::with_val_64(precision, layer.z0.as_ref());
        argument_1 += layer.d.as_ref();
        argument_1 -= z;
        argument_1 *= &reciprocal_sqrt;
        argument_1 += &sqrt_mu_a;
        argument_1.erf_mut();

        let mut argument_2 = Float::with_val_64(precision, layer.z0.as_ref());
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
pub struct FlatTopBeam<'a> {
    /// Units: cm
    pub radius: Cow<'a, Float>,
}

impl<'a> Beam for FlatTopBeam<'a> {
    fn into_owned(self) -> <Self as ToOwned>::Owned {
        Self {
            radius: Cow::Owned(self.radius.into_owned()),
        }
    }

    fn evaluate_at<'b>(
        &self,
        precision: u64,
        thermal_properties: &ThermalProperties<'b>,
        layer: &Layer<'b>,
        z: &Float,
        r: &Float,
        tp: &Float,
    ) -> Float {
        let radius = self.radius.as_ref();

        if *tp == 0 && r > radius {
            return Float::with_val_64(precision, Special::Zero);
        }

        let z_factor = LargeBeam.evaluate_at(precision, thermal_properties, layer, z, r, tp);

        if *tp == 0 {
            return z_factor;
        }

        //TODO: don't duplicate this between the code in LargeBeam and this
        //      function
        let mut alpha = Float::with_val_64(precision, thermal_properties.k.as_ref());
        alpha /= thermal_properties.rho.as_ref();
        alpha /= thermal_properties.c.as_ref();

        z_factor
            * if *r == 0 {
                let mut r_factor = Float::with_val_64(precision, radius);
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
                b *= radius;
                a *= r;

                let mut r_factor = utilities::marcum_q(1, &a, &b, precision);
                r_factor = 1 - r_factor;
                r_factor
            }
    }
}

/*pub struct MultiAbsorbingLayer {
    layers: (),
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
        let thermal_properties = ThermalProperties {
            rho: Cow::Borrowed(&ONE),
            c: Cow::Borrowed(&ONE),
            k: Cow::Borrowed(&ONE),
        };
        let layer = Layer {
            d: Cow::Borrowed(&ONE),
            z0: Cow::Borrowed(&ZERO),
            mu_a: Cow::Borrowed(&ONE),
            e0: Cow::Borrowed(&ONE),
        };

        assert_eq!(
            LargeBeam.evaluate_at(64, &thermal_properties, &layer, &ZERO, &ZERO, &ZERO),
            5e-1
        );

        let mut result = LargeBeam.evaluate_at(64, &thermal_properties, &layer, &ONE, &ZERO, &ZERO);
        // reference result: 0.5 * e^-1
        result -= 1.8393972058572116080e-1;
        result.abs_mut();
        assert!(result < *EPSILON);

        let mut result = LargeBeam.evaluate_at(64, &thermal_properties, &layer, &ONE, &ZERO, &ONE);
        // reference result: 0.5 * e^-1 * e^1 * (erf(1) - erf(-1/sqrt(4) + 1))
        result -= 1.6110045756833416583e-1;
        result.abs_mut();
        println!("{}", result);
        assert!(result < *EPSILON);
    }

    #[test]
    fn flat_top_beam_sanity() {
        let thermal_properties = ThermalProperties {
            rho: Cow::Borrowed(&ONE),
            c: Cow::Borrowed(&ONE),
            k: Cow::Borrowed(&ONE),
        };
        let layer = Layer {
            d: Cow::Borrowed(&ONE),
            z0: Cow::Borrowed(&ZERO),
            mu_a: Cow::Borrowed(&ONE),
            e0: Cow::Borrowed(&ONE),
        };
        let beam = FlatTopBeam {
            radius: Cow::Borrowed(&ONE),
        };

        assert_eq!(
            beam.evaluate_at(64, &thermal_properties, &layer, &ZERO, &ZERO, &ZERO),
            5e-1
        );

        let mut result = beam.evaluate_at(64, &thermal_properties, &layer, &ONE, &ZERO, &ZERO);
        // reference result: 0.5 * e^-1 * (1 - 0)
        result -= 1.8393972058572116080e-1;
        result.abs_mut();
        assert!(result < *EPSILON);

        let mut result = beam.evaluate_at(64, &thermal_properties, &layer, &ONE, &ZERO, &ONE);
        // reference result: 0.5 * e^-1 * e^1 * (erf(1) - erf(-1/sqrt(4) + 1))
        //                       * (1 - e^(-1/4))
        result -= 3.5635295060953884529e-2;
        result.abs_mut();
        println!("{}", result);
        assert!(result < *EPSILON);
    }
}
