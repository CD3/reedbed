// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use serde::{Deserialize, Serialize};

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
    pub z: f64,

    /// Units: cm
    pub r: f64,
}

/// Specification of time intervals
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(tag = "kind")]
pub enum Time {
    /// Intervals created by partitioning a region with a width
    Uniform {
        /// Units: s
        t0: f64,

        /// Units: s
        tmax: f64,

        /// Units: s
        dt: f64,
    },

    /// Intervals specified by their bounds
    Intervals(Vec<Interval>),
}

impl IntoIterator for Time {
    type Item = (f64, f64);
    type IntoIter = TimeIterator;

    fn into_iter(self) -> Self::IntoIter {
        match self {
            Self::Uniform { t0, tmax, dt } => TimeIterator::Uniform {
                t: t0 + dt,
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
        t0: f64,

        /// Units: s
        tmax: f64,

        /// Units: s
        dt: f64,

        /// The current time. Units: s
        t: f64,
    },

    /// Intervals specified by their bounds
    Intervals(Vec<Interval>),
}

impl Iterator for TimeIterator {
    type Item = (f64, f64);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Uniform { t0, tmax, dt, t } => {
                if t > tmax {
                    return None;
                }
                let told = *t;
                *t += *dt;

                Some((*t0, told))
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
    pub from: f64,

    /// Units: s
    pub to: f64,
}
