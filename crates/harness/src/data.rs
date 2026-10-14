//! Reading reference datasets: JSON Lines, schema 1 (see `docs/data-format.md`).

use std::collections::BTreeMap;
use std::fmt;
use std::path::Path;

use serde::Deserialize;

/// The only schema version this harness understands.
const SCHEMA: u32 = 1;

/// A malformed or unreadable dataset.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DataError(pub String);

impl fmt::Display for DataError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for DataError {}

/// A rule-assigned label on a point (assigned by the generator, never by hand).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Tag {
    /// The function has a pole or is undefined here; mpmath is not evaluated.
    Pole,
    /// The true value is finite but larger in magnitude than `f64::MAX`.
    Overflow,
    /// mpmath gave different doubles at its two precisions; point not judged.
    UnstableReference,
}

/// Accuracy demanded of every implementation of a function.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct Tolerance {
    /// Maximum relative error (see [`crate::rel_error`]).
    pub rel: f64,
    /// Why this tolerance is right for this function.
    pub reason: String,
}

/// Version of a reference tool.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct ToolVersion {
    /// Version string, e.g. `1.4.1`.
    pub version: String,
}

/// The first line of a dataset: what it contains and how it was made.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct Header {
    /// Format version.
    pub schema: u32,
    /// Dataset name, e.g. `erfcx`.
    pub function: String,
    /// What the function is.
    pub description: String,
    /// Input names, in order.
    pub inputs: Vec<String>,
    /// Output names, in order.
    pub outputs: Vec<String>,
    /// Number of point records that follow.
    pub points: usize,
    /// Accuracy every implementation must meet.
    pub tolerance: Tolerance,
    /// How the inputs were chosen.
    pub grid: String,
    /// mpmath version used as the arbiter.
    pub mpmath: ToolVersion,
    /// SciPy version used as the practical reference.
    pub scipy: ToolVersion,
}

/// A value in the data: a JSON number, or `"nan"`, `"inf"`, `"-inf"`.
#[derive(Deserialize)]
#[serde(untagged)]
enum Num {
    Finite(f64),
    Special(String),
}

impl Num {
    fn value(&self) -> Result<f64, DataError> {
        match self {
            Self::Finite(v) => Ok(*v),
            Self::Special(s) => match s.as_str() {
                "nan" => Ok(f64::NAN),
                "inf" => Ok(f64::INFINITY),
                "-inf" => Ok(f64::NEG_INFINITY),
                other => Err(DataError(format!("bad number {other:?}"))),
            },
        }
    }
}

type RawValues = BTreeMap<String, Vec<Num>>;

#[derive(Deserialize)]
struct RawPoint {
    #[serde(rename = "in")]
    inputs: RawValues,
    mpmath: RawValues,
    scipy: RawValues,
    #[serde(default)]
    tags: Vec<Tag>,
}

/// Named values: each name maps to a list (length 1 for scalars).
pub type Values = BTreeMap<String, Vec<f64>>;

/// One reference point.
#[derive(Debug, Clone, PartialEq)]
pub struct Point {
    /// Input values by name.
    pub inputs: Values,
    /// mpmath (arbiter) outputs by name, correctly rounded to `f64`.
    pub mpmath: Values,
    /// SciPy outputs by name.
    pub scipy: Values,
    /// Rule-assigned tags.
    pub tags: Vec<Tag>,
}

impl Point {
    /// The scalar input `name`.
    ///
    /// # Panics
    ///
    /// If `name` is not an input of length 1 (a bug in the check crate).
    #[must_use]
    pub fn scalar(&self, name: &str) -> f64 {
        match self.inputs.get(name).map(Vec::as_slice) {
            Some([v]) => *v,
            _ => panic!("no scalar input {name:?}"),
        }
    }

    /// All inputs flattened in the header's input order.
    #[must_use]
    pub fn flat_inputs(&self, header: &Header) -> Vec<f64> {
        header
            .inputs
            .iter()
            .flat_map(|n| self.inputs[n].iter().copied())
            .collect()
    }
}

/// A parsed, validated dataset.
#[derive(Debug, Clone, PartialEq)]
pub struct Dataset {
    /// Metadata line.
    pub header: Header,
    /// Point records, in file order.
    pub points: Vec<Point>,
}

fn convert(raw: &RawValues, names: &[String], what: &str) -> Result<Values, DataError> {
    if raw.len() != names.len() || names.iter().any(|n| !raw.contains_key(n)) {
        return Err(DataError(format!("{what}: expected names {names:?}")));
    }
    raw.iter()
        .map(|(k, v)| {
            Ok((
                k.clone(),
                v.iter().map(Num::value).collect::<Result<_, _>>()?,
            ))
        })
        .collect()
}

impl Dataset {
    /// Parses and validates dataset text.
    ///
    /// # Errors
    ///
    /// On malformed JSON, an unknown schema, a wrong point count, missing
    /// names, or a non-finite mpmath value at an untagged point.
    pub fn parse(text: &str) -> Result<Self, DataError> {
        let mut lines = text.lines();
        let first = lines
            .next()
            .ok_or_else(|| DataError("empty dataset".into()))?;
        let header: Header =
            serde_json::from_str(first).map_err(|e| DataError(format!("header: {e}")))?;
        if header.schema != SCHEMA {
            return Err(DataError(format!("unsupported schema {}", header.schema)));
        }
        let mut points = Vec::new();
        for (i, line) in lines.enumerate() {
            let at = |m: String| DataError(format!("{} point {i}: {m}", header.function));
            let raw: RawPoint = serde_json::from_str(line).map_err(|e| at(e.to_string()))?;
            let point = Point {
                inputs: convert(&raw.inputs, &header.inputs, "inputs").map_err(|e| at(e.0))?,
                mpmath: convert(&raw.mpmath, &header.outputs, "mpmath").map_err(|e| at(e.0))?,
                scipy: convert(&raw.scipy, &header.outputs, "scipy").map_err(|e| at(e.0))?,
                tags: raw.tags,
            };
            let finite = point.mpmath.values().flatten().all(|v| v.is_finite());
            if !finite && point.tags.is_empty() {
                return Err(at("non-finite mpmath value at an untagged point".into()));
            }
            points.push(point);
        }
        if points.len() != header.points {
            return Err(DataError(format!(
                "header says {} points, found {}",
                header.points,
                points.len()
            )));
        }
        Ok(Self { header, points })
    }

    /// Loads `reference/data/<name>.jsonl` from `root`.
    ///
    /// # Errors
    ///
    /// If the file cannot be read or fails [`Dataset::parse`].
    pub fn load(root: &Path, name: &str) -> Result<Self, DataError> {
        let path = root.join("reference/data").join(format!("{name}.jsonl"));
        let text = std::fs::read_to_string(&path)
            .map_err(|e| DataError(format!("{}: {e}", path.display())))?;
        Self::parse(&text)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub(crate) const HEADER: &str = r#"{"schema":1,"function":"f","description":"d","inputs":["x"],"outputs":["y"],"points":2,"tolerance":{"rel":1e-13,"reason":"r"},"grid":"g","mpmath":{"version":"1.4.1","dps":[50,100]},"scipy":{"version":"1.18.1","function":"s"},"numpy":"2","python":"3.13"}"#;

    pub(crate) fn sample() -> String {
        format!(
            "{HEADER}\n{}\n{}\n",
            r#"{"in":{"x":[0.5]},"mpmath":{"y":[0.1]},"scipy":{"y":[0.1]}}"#,
            r#"{"in":{"x":[-30]},"mpmath":{"y":["inf"]},"scipy":{"y":["-inf"]},"tags":["overflow"]}"#
        )
    }

    #[test]
    fn parses_numbers_specials_and_tags() {
        let d = Dataset::parse(&sample()).unwrap();
        assert_eq!(d.header.function, "f");
        assert_eq!(d.points.len(), 2);
        assert_eq!(d.points[0].scalar("x").to_bits(), 0.5f64.to_bits());
        assert_eq!(d.points[1].mpmath["y"], vec![f64::INFINITY]);
        assert_eq!(d.points[1].scipy["y"], vec![f64::NEG_INFINITY]);
        assert_eq!(d.points[1].tags, vec![Tag::Overflow]);
        assert_eq!(d.points[1].flat_inputs(&d.header), vec![-30.0]);
    }

    #[test]
    fn parses_decimals_exactly() {
        // 0.1 must become the nearest double; float_roundtrip guarantees it.
        let text = sample().replace("0.1]", "2.2250738585072011e-308]");
        let d = Dataset::parse(&text).unwrap();
        assert_eq!(
            d.points[0].mpmath["y"][0].to_bits(),
            2.225_073_858_507_201e-308_f64.to_bits()
        );
    }

    #[test]
    fn nan_string_parses() {
        let text = sample().replace(r#""inf"]"#, r#""nan"]"#);
        assert!(Dataset::parse(&text).unwrap().points[1].mpmath["y"][0].is_nan());
    }

    fn err(text: &str) -> String {
        Dataset::parse(text).unwrap_err().to_string()
    }

    #[test]
    fn rejects_malformed_data() {
        assert_eq!(err(""), "empty dataset");
        assert!(err("{").starts_with("header:"));
        assert!(err(&sample().replace(r#""schema":1"#, r#""schema":2"#)).contains("schema 2"));
        assert!(err(&sample().replace(r#""points":2"#, r#""points":3"#)).contains("found 2"));
        assert!(err(&sample().replace(r#""inf"]"#, r#""big"]"#)).contains("bad number"));
        assert!(err(&sample().replace(r#""x":[0.5]"#, r#""z":[0.5]"#)).contains("inputs"));
        assert!(
            err(&sample().replace(r#"{"y":[0.1]},"scipy""#, r#"{},"scipy""#)).contains("mpmath")
        );
        assert!(
            err(&sample().replace(r#""scipy":{"y":[0.1]}"#, r#""scipy":{}"#)).contains("scipy")
        );
        assert!(err(&sample().replace(r#","tags":["overflow"]"#, "")).contains("untagged"));
        assert!(err(&format!("{HEADER}\nnot json\n")).contains("point 0"));
    }

    #[test]
    #[should_panic(expected = "no scalar input \"q\"")]
    fn scalar_panics_on_unknown_input() {
        let d = Dataset::parse(&sample()).unwrap();
        let _ = d.points[0].scalar("q");
    }

    #[test]
    fn load_reads_from_reference_dir_and_reports_missing_files() {
        let dir = std::env::temp_dir().join(format!("rsv-load-{}", std::process::id()));
        std::fs::create_dir_all(dir.join("reference/data")).unwrap();
        std::fs::write(dir.join("reference/data/f.jsonl"), sample()).unwrap();
        assert_eq!(Dataset::load(&dir, "f").unwrap().points.len(), 2);
        assert!(
            Dataset::load(&dir, "missing")
                .unwrap_err()
                .0
                .contains("missing.jsonl")
        );
        std::fs::remove_dir_all(dir).unwrap();
    }
}
