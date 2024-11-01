// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

#[derive(thiserror::Error, Debug)]
pub enum Greens {
    #[error("One or more layers overlap")]
    LayerOverlap,
}
