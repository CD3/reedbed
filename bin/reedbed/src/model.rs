// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use rug::Float;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

use crate::utilities::FloatParseDisplay;
use reedbed_lib::{
    errors,
    greens::{self, LargeBeam},
};

//NOTE: this module contains copies of all of the structures in
//      lib/src/greens.rs so that we can manually control serialization and
//      deserialization behavior
//
//      of particular note are the changes done to rug::Float handling. we
//      override the way in which these are handled because rug's serde
//      implementations for Floats (and presumably the other types) serialize
//      them like this (in JSON):
//
//      {
//          "prec": /* precision of the floating point value in bits */,
//          "radix": /* base in which the floating point value is stored */,
//          "value": /* the value of the float*/,
//      }
//
//      the override merely changes this behavior to storing the float as a
//      string instead, and assumes when reading that it is 64-bit precision.
//      if this behavior needs to be changed it can be revisited in the future

/// Configuration of the entire simulation
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct Configuration {
    /// Thermal properties
    pub thermal: ThermalProperties,

    /// Retinal tissue layers
    pub layers: Layers,

    /// Beam parameters
    pub laser: Beam,

    /// Simulation parameters
    pub simulation: Simulation,
}

/// A configuration structure for specific thermal properties
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct ThermalProperties {
    /// Units: g*cm^3
    #[serde(with = "FloatParseDisplay")]
    pub rho: Float,

    /// Units: J*g^-1*K^-1
    #[serde(with = "FloatParseDisplay")]
    pub c: Float,

    /// Units: W*cm^-1*K^-1
    #[serde(with = "FloatParseDisplay")]
    pub k: Float,
}

impl ThermalProperties {
    /// Convert this structure into its equivalent in `reedbed_lib`
    pub fn into_lib(self) -> greens::ThermalProperties<'static> {
        greens::ThermalProperties {
            rho: Cow::Owned(self.rho),
            c: Cow::Owned(self.c),
            k: Cow::Owned(self.k),
        }
    }
}

/// A layer of tissue
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct Layer {
    /// Units: cm
    #[serde(with = "FloatParseDisplay")]
    pub d: Float,

    /// Units: cm
    #[serde(with = "FloatParseDisplay")]
    pub z0: Float,

    /// Units: cm^-1
    #[serde(with = "FloatParseDisplay")]
    pub mu_a: Float,

    /// Irradiance. Units: W*cm^-2
    #[serde(with = "FloatParseDisplay")]
    pub e0: Float,
}

impl Layer {
    /// Convert this structure into its equivalent in `reedbed_lib`
    pub fn into_lib(self) -> greens::Layer<'static> {
        greens::Layer {
            d: Cow::Owned(self.d),
            z0: Cow::Owned(self.z0),
            mu_a: Cow::Owned(self.mu_a),
            e0: Cow::Owned(self.e0),
        }
    }
}

/// Multiple layers of tissue
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(transparent)]
pub struct Layers {
    /// The layers this [`struct@Layers`] is composed of
    pub layers: Vec<Layer>,
}

impl Layers {
    /// Convert this structure into its equivalent in `reedbed_lib`
    pub fn into_lib(self) -> Result<greens::Layers, errors::Greens> {
        greens::Layers::new_static(self.layers.into_iter().map(|layer| layer.into_lib()))
    }
}

/// Container for different kinds of [`trait@Beam`]s
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(tag = "kind")]
pub enum Beam {
    LargeBeam(LargeBeam),
    FlatTopBeam(FlatTopBeam),
}

//TODO: get a better description both here and in lib/src/greens.rs
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct FlatTopBeam {
    /// Units: cm
    #[serde(with = "FloatParseDisplay")]
    pub radius: Float,
}

impl FlatTopBeam {
    /// Convert this structure into its equivalent in `reedbed_lib`
    pub fn into_lib(self) -> greens::FlatTopBeam<'static> {
        greens::FlatTopBeam {
            radius: Cow::Owned(self.radius),
        }
    }
}

/// Simulation parameters
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct Simulation {
    /// Sensor parameters
    pub sensor: Sensor,

    /// Specification of time intervals to sample values at
    pub time: Time,
    //TODO: put quadrature and epsilon configuration for numerical quadrature
    //      here in a nested struct
}

/// Sensor parameters
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct Sensor {
    /// Units: cm
    #[serde(with = "FloatParseDisplay")]
    pub z: Float,

    /// Units: cm
    #[serde(with = "FloatParseDisplay")]
    pub r: Float,
}

/// Specification of time intervals
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(tag = "kind")]
pub enum Time {
    /// Intervals created by partitioning a region with a width
    Uniform {
        /// Units: s
        #[serde(with = "FloatParseDisplay")]
        t0: Float,

        /// Units: s
        #[serde(with = "FloatParseDisplay")]
        tmax: Float,

        /// Units: s
        #[serde(with = "FloatParseDisplay")]
        dt: Float,
    },

    /// Intervals specified by their bounds
    Intervals(Vec<Interval>),
}

impl IntoIterator for Time {
    type Item = (Float, Float);
    type IntoIter = TimeIterator;

    fn into_iter(self) -> Self::IntoIter {
        match self {
            Self::Uniform { t0, tmax, dt } => TimeIterator::Uniform {
                t: t0.clone() + &dt,
                t0,
                tmax,
                dt,
            },
            Self::Intervals(intervals) => TimeIterator::Intervals(intervals),
        }
    }
}

pub enum TimeIterator {
    /// Intervals created by partitioning a region with a width
    Uniform {
        /// Units: s
        t0: Float,

        /// Units: s
        tmax: Float,

        /// Units: s
        dt: Float,

        /// The current time. Units: s
        t: Float,
    },

    /// Intervals specified by their bounds
    Intervals(Vec<Interval>),
}

impl Iterator for TimeIterator {
    type Item = (Float, Float);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Uniform { t0, tmax, dt, t } => {
                if t > tmax {
                    return None;
                }
                let told = t.clone();
                *t += &*dt;

                Some((t0.clone(), told))
            }
            Self::Intervals(intervals) => {
                intervals.pop().map(|Interval { from, to }| (from, to))
            }
        }
    }
}

/// An interval [from, to) specified by its bounds
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct Interval {
    /// Units: s
    #[serde(with = "FloatParseDisplay")]
    pub from: Float,

    /// Units: s
    #[serde(with = "FloatParseDisplay")]
    pub to: Float,
}
