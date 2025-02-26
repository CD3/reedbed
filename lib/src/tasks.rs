// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use serde::{Deserialize, Serialize};

use crate::{
    greens::{FlatTopBeam, LargeBeam, Layers, ThermalProperties},
    multiple_pulse::Pulse,
    simulation::Simulation,
};

/// Description of a simulation task
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(tag = "kind")]
pub enum Operation {
    /// Computes the temperature rise resulting from laser exposure in retinal
    /// tissue at given points in time
    TemperatureRise(TemperatureRise),

    /// Computes the temperature rise resulting from multiple pulses
    MultiplePulse(MultiplePulse),
}

/// Structure containing the information necessary to perform a temperature
/// rise computation
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct TemperatureRise {
    /// Thermal properties
    pub thermal: ThermalProperties,

    /// Retinal tissue layers
    pub layers: Layers,

    /// Beam parameters
    pub laser: Beam,

    /// Simulation parameters
    pub simulation: Simulation,
}

/// Structure containing the information necessary to perform a multiple pulse
/// computation
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct MultiplePulse {
    /// Thermal properties
    pub thermal: ThermalProperties,

    /// Retinal tissue layers
    pub layers: Layers,

    /// Beam parameters
    pub laser: Beam,

    /// Simulation parameters
    pub simulation: Simulation,

    /// Pulse information
    pub pulses: Vec<Pulse>,
}

/// Container for different kinds of [`trait@Beam`]s
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
#[serde(tag = "kind")]
pub enum Beam {
    LargeBeam(LargeBeam),
    FlatTopBeam(FlatTopBeam),
}

pub fn init() {
    // do something here with crossbeam and some thread spawning
}
