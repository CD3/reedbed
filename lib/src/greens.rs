// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use serde::{Deserialize, Serialize};
use statrs::function::erf::erf;

use crate::{errors::Greens as Error, quadrature::Quadrature, utilities};

/// A configuration structure for specific thermal properties
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct ThermalProperties {
    /// Units: g*cm^3
    pub rho: f64,

    /// Units: J*g^-1*K^-1
    pub c: f64,

    /// Units: W*cm^-1*K^-1
    pub k: f64,
}

/// A layer of tissue
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct Layer {
    /// Units: cm
    pub d: f64,

    /// Units: cm
    pub z0: f64,

    /// Units: cm^-1
    pub mu_a: f64,

    /// Irradiance. Units: W*cm^-2
    pub e0: f64,
}

/// Multiple layers of tissue
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct Layers {
    /// The layers this [`struct@Layers`] is composed of
    layers: Vec<Layer>,
}

impl Layers {
    /// Creates a new [`struct@Layers`] from multiple [`struct@Layer`]s
    ///
    /// If the input layers are not sorted in order of incidence, they are
    /// sorted. Irradiance is taken from the topmost layer and propagated
    /// downward according to Beer's Law
    ///
    /// If the input layers overlap in any way, an error is returned
    pub fn new(
        input_layers: impl IntoIterator<Item = Layer>,
    ) -> Result<Self, Error> {
        let mut layers = input_layers.into_iter().collect::<Vec<_>>();
        layers.sort_by(|a, b| a.z0.total_cmp(&b.z0));

        if let Some(layer) = layers.first() {
            let mut e0 = layer.e0;
            let mut z0 = layer.z0 + layer.d;
            let mut b = (layer.d * layer.mu_a * -1.00).exp();
            e0 *= b;

            for layer in layers.iter_mut().skip(1) {
                if layer.z0 < z0 {
                    return Err(Error::LayerOverlap);
                }

                layer.e0 = e0;
                z0 = layer.z0 + layer.d;
                b = (layer.d * layer.mu_a * -1.00).exp();
                e0 *= b;
            }
        }

        Ok(Self { layers })
    }

    /// Updates the irradiance value, in W*cm^-2
    ///
    /// The input irradiance is set as-is on the topmost layer and is
    /// propagated downward according to Beer's Law
    pub fn set_e0(&mut self, mut e0: f64) {
        if let Some(layer) = self.layers.first_mut() {
            layer.e0 = e0;
            let mut b = (layer.d * layer.mu_a * -1.00).exp();
            e0 *= b;

            for layer in self.layers.iter_mut().skip(1) {
                layer.e0 = e0;

                b = (layer.d * layer.mu_a * -1.00).exp();
                e0 *= b;
            }
        }
    }

    /// Runs the given [`trait@Beam`] over the contained [`struct@Layer`]s
    /// with the provided [`struct@ThermalProperties`]
    ///
    /// Not all implementations of [`trait@Beam`] will use all parameters
    pub fn evaluate_with(
        &self,
        beam: &impl Beam,
        thermal_properties: &ThermalProperties,
        z: f64,
        r: f64,
        tp: f64,
    ) -> f64 {
        let mut sum = 0.00;

        for layer in &self.layers {
            sum += beam.evaluate_with(thermal_properties, layer, z, r, tp);
        }

        sum
    }

    /// Calculates the temperature rise over the interval a..b
    ///
    /// Similar to [`fn@temperature_rise`], this is really just a convenience
    /// wrapper over `Quadrature::integrate`
    #[allow(clippy::too_many_arguments)]
    pub fn temperature_rise(
        &self,
        quadrature: &impl Quadrature<f64>,
        beam: &impl Beam,
        thermal_properties: &ThermalProperties,
        z: f64,
        r: f64,
        epsilon: f64,
        bounds: (f64, f64),
    ) -> (f64, f64) {
        quadrature.integrate(
            |t| self.evaluate_with(beam, thermal_properties, z, r, t),
            epsilon,
            bounds,
        )
    }
}

//TODO: we could probably swap the use of [`struct@Float`] for a generic
//      parameter that implements the operation traits in rug::ops in most
//      (if not all) places

/// An abstraction over the various `*Beam` structures
pub trait Beam {
    /// Run the beam over a given [`struct@Layer`] with the provided
    /// [`struct@ThermalProperties`]
    ///
    /// Not all implementations of [`trait@Beam`] will use all parameters
    fn evaluate_with(
        &self,
        thermal_properties: &ThermalProperties,
        layer: &Layer,
        z: f64,
        r: f64,
        tp: f64,
    ) -> f64;
}

//TODO: add documentation describing what this is
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct LargeBeam;

impl Beam for LargeBeam {
    //TODO: it (might?) be worthwhile to have a specialized method that
    //      doesn't need to take r. however, this could also be addressed with
    //      the genericization of this method at the trait level. see above
    //      for more details
    fn evaluate_with(
        &self,
        thermal_properties: &ThermalProperties,
        layer: &Layer,
        z: f64,
        _r: f64,
        tp: f64,
    ) -> f64 {
        //TODO: make this less naive

        let mut alpha = thermal_properties.k
            / thermal_properties.rho
            / thermal_properties.c;
        let mut term_1 = (layer.mu_a * layer.e0)
            / thermal_properties.rho
            / thermal_properties.c
            / 2.0;
        let mut term_2 = ((z - layer.z0) * layer.mu_a * -1.00).exp();

        if tp < 1e-9 {
            return term_1 * term_2;
        }

        let mut term_3 = (layer.mu_a.powi(2) * tp * alpha).exp();
        let mut reciprocal_sqrt = (alpha * tp * 4.00).sqrt().recip();
        let mut sqrt_mu_a = (alpha * tp).sqrt() * layer.mu_a;
        let mut argument_1 =
            erf(((layer.z0 + layer.d - z) * reciprocal_sqrt) + sqrt_mu_a);
        let mut argument_2 =
            erf(((layer.z0 - z) * reciprocal_sqrt) + sqrt_mu_a);
        let mut term_4 = argument_1 - argument_2;

        term_1 * term_2 * term_3 * term_4
    }
}

#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct FlatTopBeam {
    /// Units: cm
    pub radius: f64,
}

impl Beam for FlatTopBeam {
    //TODO: same todo as above
    fn evaluate_with(
        &self,
        thermal_properties: &ThermalProperties,
        layer: &Layer,
        z: f64,
        r: f64,
        tp: f64,
    ) -> f64 {
        if tp < 1e-9 && r > self.radius {
            return 0.00;
        }

        let z_factor =
            LargeBeam.evaluate_with(thermal_properties, layer, z, r, tp);

        if tp < 1e-9 {
            return z_factor;
        }

        let mut alpha = thermal_properties.k
            / thermal_properties.rho
            / thermal_properties.c;

        z_factor
            * if r < 1e-9 {
                1.00 - (self.radius.powi(2) / -4.00 / alpha / tp).exp()
            } else {
                //TODO: this is not accurate at all. fix the marcum-q function
                //      implementation
                let a = (2.0 * alpha * tp).recip();
                1.00 - utilities::marcum_q(1, a, a * self.radius * r)
            }
    }
}

/// Calculates the temperature rise over the interval a..b
///
/// This is really just a convenience wrapper around `Quadrature::integrate`
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn temperature_rise(
    quadrature: &impl Quadrature<f64>,
    beam: &impl Beam,
    thermal_properties: &ThermalProperties,
    layer: &Layer,
    z: f64,
    r: f64,
    epsilon: f64,
    bounds: (f64, f64),
) -> (f64, f64) {
    quadrature.integrate(
        |t| beam.evaluate_with(thermal_properties, layer, z, r, t),
        epsilon,
        bounds,
    )
}

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
            LargeBeam.evaluate_with(
                64,
                &thermal_properties,
                &layer,
                &ZERO,
                &ZERO,
                &ZERO
            ),
            5e-1
        );

        let mut result = LargeBeam.evaluate_with(
            64,
            &thermal_properties,
            &layer,
            &ONE,
            &ZERO,
            &ZERO,
        );
        // reference result: 0.5 * e^-1
        result -= 1.8393972058572116080e-1;
        result.abs_mut();
        assert!(result < *EPSILON);

        let mut result = LargeBeam.evaluate_with(
            64,
            &thermal_properties,
            &layer,
            &ONE,
            &ZERO,
            &ONE,
        );
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
            beam.evaluate_with(
                64,
                &thermal_properties,
                &layer,
                &ZERO,
                &ZERO,
                &ZERO
            ),
            5e-1
        );

        let mut result = beam.evaluate_with(
            64,
            &thermal_properties,
            &layer,
            &ONE,
            &ZERO,
            &ZERO,
        );
        // reference result: 0.5 * e^-1 * (1 - 0)
        result -= 1.8393972058572116080e-1;
        result.abs_mut();
        assert!(result < *EPSILON);

        let mut result = beam.evaluate_with(
            64,
            &thermal_properties,
            &layer,
            &ONE,
            &ZERO,
            &ONE,
        );
        // reference result: 0.5 * e^-1 * e^1 * (erf(1) - erf(-1/sqrt(4) + 1))
        //                       * (1 - e^(-1/4))
        result -= 3.5635295060953884529e-2;
        result.abs_mut();
        println!("{}", result);
        assert!(result < *EPSILON);
    }

    #[test]
    fn layers_sanity() {
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
        let layers = Layers::new([layer.clone()])
            .expect("unable to construct a Layers structure");

        let mut result = layers.evaluate_with(
            64,
            &LargeBeam,
            &thermal_properties,
            &ONE,
            &ZERO,
            &ONE,
        );
        result -= LargeBeam.evaluate_with(
            64,
            &thermal_properties,
            &layer,
            &ONE,
            &ZERO,
            &ONE,
        );
        assert!(result < *EPSILON);

        let layers = Layers::new([
            Layer {
                d: Cow::Borrowed(&ONE),
                z0: Cow::Borrowed(&ZERO),
                mu_a: Cow::Borrowed(&ONE),
                e0: Cow::Borrowed(&ONE),
            },
            Layer {
                d: Cow::Borrowed(&ONE),
                z0: Cow::Borrowed(&ONE),
                mu_a: Cow::Borrowed(&ONE),
                e0: Cow::Borrowed(&ZERO),
            },
        ])
        .expect("unable to construct a Layers structure");

        let layer = Layer {
            d: Cow::Owned(Float::with_val_64(64, 2)),
            z0: Cow::Borrowed(&ZERO),
            mu_a: Cow::Borrowed(&ONE),
            e0: Cow::Borrowed(&ONE),
        };

        let beam = FlatTopBeam {
            radius: Cow::Borrowed(&ONE),
        };

        let small = Float::with_val_64(64, 1e-6);

        let mut result = layers.evaluate_with(
            64,
            &beam,
            &thermal_properties,
            &ZERO,
            &ZERO,
            &small,
        );
        result -= beam.evaluate_with(
            64,
            &thermal_properties,
            &layer,
            &ZERO,
            &ZERO,
            &small,
        );
        assert!(result < *EPSILON);
    }
}
