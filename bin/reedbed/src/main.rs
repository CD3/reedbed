// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

#![allow(clippy::cognitive_complexity)]
#![warn(clippy::cargo_common_metadata)]
#![warn(clippy::dbg_macro)]
#![warn(clippy::explicit_deref_methods)]
#![warn(clippy::filetype_is_file)]
#![warn(clippy::imprecise_flops)]
#![warn(clippy::large_stack_arrays)]
#![warn(clippy::todo)]
#![warn(clippy::unimplemented)]
#![deny(clippy::await_holding_lock)]
#![deny(clippy::cast_lossless)]
#![deny(clippy::clone_on_ref_ptr)]
#![deny(clippy::doc_markdown)]
#![deny(clippy::empty_enum)]
#![deny(clippy::enum_glob_use)]
#![deny(clippy::exit)]
#![deny(clippy::explicit_into_iter_loop)]
#![deny(clippy::explicit_iter_loop)]
#![deny(clippy::fallible_impl_from)]
#![deny(clippy::inefficient_to_string)]
#![deny(clippy::large_digit_groups)]
#![deny(clippy::wildcard_dependencies)]
#![deny(clippy::wildcard_imports)]
#![deny(clippy::unused_self)]
#![deny(clippy::single_match_else)]
#![deny(clippy::option_option)]
#![deny(clippy::mut_mut)]

use anyhow::Context;
use clap::{Parser, Subcommand, ValueEnum};
use std::{
    fs::File,
    io::{self, BufReader, Read, Write},
    str::FromStr,
};
use strum_macros::{Display, EnumString};

use reedbed_lib::{
    multiple_pulse::{self, Pulse},
    quadrature,
    tasks::{Beam, MultiplePulse, Operation, TemperatureRise},
};

#[global_allocator]
static GLOBAL_ALLOCATOR: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Command line interface for a Green's function based model for calculating
/// temperature rise resulting from laser exposure in retinal tissue
#[derive(Parser, Debug)]
#[clap(
    author = "superwhiskers <whiskerdev@protonmail.com>",
    version = "0.0.0"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// The size of the thread pool to run jobs in
    ///
    /// If no value is provided, either the value in `RAYON_NUM_THREADS` or
    /// the number of logical cores will be used
    #[arg(short = 'j', long = "jobs", value_name = "JOBS")]
    threads: Option<usize>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Start executing configurations input either through a file or stdin
    Start {
        /// The input file describing the computation. If none is provided,
        /// then stdin is used
        #[arg(value_name = "IN")]
        input: Option<String>,

        /// The format to read computation descriptions in
        ///
        /// If none is provided, it will be automatically detected from the
        /// file extension (if there is one). If a format is not able to be
        /// determined, no computations will be performed
        #[arg(value_enum, short = 'f', long = "from", value_name = "FROM")]
        from: Option<InputFormat>,

        /// The output file to write computation results to. If none is
        /// provided, stdout is used
        #[arg(short = 'o', long = "output", value_name = "OUT")]
        output: Option<String>,

        /// The format to write computation results in
        ///
        /// If none is provided, it will be automatically detected from the
        /// file extension (if there is one). If a format is not able to be
        /// determined, no computations will be performed
        #[arg(value_enum, short = 't', long = "to", value_name = "TO")]
        to: Option<OutputFormat>,
    },
}

#[derive(ValueEnum, Debug, Clone, Display, EnumString)]
#[strum(serialize_all = "kebab-case")]
pub enum InputFormat {
    /// JSON
    Json,
}

#[derive(ValueEnum, Debug, Clone, Display, EnumString)]
#[strum(serialize_all = "kebab-case")]
pub enum OutputFormat {
    /// JSON
    Json,
}

#[derive(thiserror::Error, Debug)]
enum ReedbedError {
    #[error("The input format is ambiguous")]
    AmbiguousInput,

    #[error("The output format is ambiguous")]
    AmbiguousOutput,
}

fn main() -> anyhow::Result<()> {
    let app = Cli::parse();

    match app.command {
        Commands::Start {
            ref input,
            from,
            ref output,
            to,
        } => {
            let input_stream = BufReader::new(if let Some(name) = input {
                Box::new(File::open(name).context(
                    "Unable to open the file containing the computation",
                )?) as Box<dyn Read>
            } else {
                Box::new(io::stdin()) as Box<dyn Read>
            });

            let input = match from
                .or_else(|| {
                    input.as_ref().and_then(|n| n.rsplit_once('.')).and_then(
                        |s| <InputFormat as FromStr>::from_str(s.1).ok(),
                    )
                })
                .ok_or(ReedbedError::AmbiguousInput)?
            {
                InputFormat::Json => {
                    serde_json::Deserializer::from_reader(input_stream)
                        .into_iter::<Operation>()
                }
            };

            let output_stream = if let Some(name) = output {
                Box::new(File::create(name).context(
                    "Unable to open the file intended to hold the results of the computation",
                )?) as Box<dyn Write>
            } else {
                Box::new(io::stdout()) as Box<dyn Write>
            };

            let output = match to
                .or_else(|| {
                    output.as_ref().and_then(|n| n.rsplit_once('.')).and_then(
                        |s| <OutputFormat as FromStr>::from_str(s.1).ok(),
                    )
                })
                .ok_or(ReedbedError::AmbiguousOutput)?
            {
                OutputFormat::Json => {
                    serde_json::Serializer::new(output_stream)
                }
            };

            // initialize task system here

            for configuration in input {
                let configuration = configuration.context(
                    "Unable to deserialize an input configuration",
                )?;
                let quadrature = quadrature::TanhSinh { iteration_limit: 6 };
                let epsilon = 1e-5;

                match configuration {
                    Operation::TemperatureRise(TemperatureRise {
                        thermal,
                        layers,
                        laser,
                        simulation,
                    }) => match laser {
                        Beam::LargeBeam(_) => {}
                        Beam::FlatTopBeam(beam) => {
                            for (a, b) in simulation.time {
                                let (rise, _err) = layers.temperature_rise(
                                    &quadrature,
                                    &beam,
                                    &thermal,
                                    simulation.sensor.z,
                                    simulation.sensor.r,
                                    epsilon,
                                    (a, b),
                                );

                                println!("{b} {rise}");
                            }
                        }
                    },
                    Operation::MultiplePulse(MultiplePulse {
                        thermal,
                        layers,
                        laser,
                        simulation,
                        pulses,
                    }) => {
                        let regions = multiple_pulse::extract_temperature_rise_computations(&pulses, simulation.time.clone());

                        match laser {
                            Beam::LargeBeam(_) => {}
                            Beam::FlatTopBeam(beam) => {
                                for ((_, b), pulses) in
                                    simulation.time.into_iter().zip(regions)
                                {
                                    let mut rise = 0.0;

                                    for Pulse {
                                        arrival_time,
                                        duration,
                                        scale,
                                    } in pulses
                                    {
                                        let (individual_rise, _err) = layers
                                            .temperature_rise(
                                                &quadrature,
                                                &beam,
                                                &thermal,
                                                simulation.sensor.z,
                                                simulation.sensor.r,
                                                epsilon,
                                                (0.0, duration),
                                            );

                                        rise += scale * individual_rise;
                                    }

                                    println!("{b} {rise}");
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    Ok(())
}
