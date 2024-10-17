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

use clap::{Parser, Subcommand};
use rug::Float;

use reedbed_lib::utilities;

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
}

#[derive(Subcommand, Debug)]
enum Commands {}

fn main() -> anyhow::Result<()> {
    let _app = Cli::parse();

    let a = Float::with_val_64(64, 20.0);
    let b = Float::with_val_64(64, 20.0);

    println!("{:?}", utilities::marcum_q(1, &a, &b, 64));

    Ok(())
}
