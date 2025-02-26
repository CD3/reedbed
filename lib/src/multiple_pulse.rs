// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use serde::{Deserialize, Serialize};

/// Description of a pulse done by a beam
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct Pulse {
    /// Units: s
    pub arrival_time: f64,

    /// Units: s
    pub duration: f64,

    /// Unitless scaling factor
    pub scale: f64,
}

/// Extract temperature rise computation tasks given a slice of pulses and an
/// iterator over time intervals of the form [a, b)
///
/// Tasks are returned as a `Vec<Vec<Pulse>` in which the first index
/// represents the interval which the vector of pulses is associated with
pub fn extract_temperature_rise_computations(
    pulses: &[Pulse],
    intervals: impl IntoIterator<Item = (f64, f64)>,
) -> Vec<Vec<Pulse>> {
    let mut output = vec![];

    for (_, b) in intervals {
        let mut contained_pulses = vec![];

        for Pulse {
            arrival_time,
            duration,
            scale,
        } in pulses
        {
            match (arrival_time, arrival_time + duration) {
                (arrival_time, _) if *arrival_time >= b => continue,
                (arrival_time, end_time) if end_time <= b => {
                    // the pulse ends within or before the interval
                    contained_pulses.push(Pulse {
                        arrival_time: *arrival_time,
                        duration: b - arrival_time,
                        scale: *scale,
                    });
                    contained_pulses.push(Pulse {
                        arrival_time: end_time,
                        duration: b - end_time,
                        scale: -scale,
                    });
                }
                (arrival_time, _) => {
                    // the pulse ends after the interval ends
                    contained_pulses.push(Pulse {
                        arrival_time: *arrival_time,
                        duration: b - arrival_time,
                        scale: *scale,
                    });
                }
            }
        }

        output.push(contained_pulses);
    }

    output
}
