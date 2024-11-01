// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use rug::{float::ParseFloatError, Float};

serde_with::serde_conv!(
    pub FloatParseDisplay,
    Float,
    |float: &Float| float.to_string(),
    |value: String| -> Result<Float, ParseFloatError> {
        Float::parse(value).map(|f| Float::with_val(64, f))
    }
);
