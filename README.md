# rust-sci-verify

Verification of Rust scientific-computing crates against independent references.

Every function we check is evaluated at thousands of points and compared with

- **mpmath** at 50 and 100 significant digits: the arbiter of the true value;
- **SciPy**: the practical reference most users compare against.

Each point is classified as *pass*, *crate wrong*, *SciPy wrong*, *both wrong*,
or *convention difference* (see [docs/conventions.md](docs/conventions.md)).
Confirmed bugs become reproducible reports in [`findings/`](findings/), and the
per-crate results are in [`reports/`](reports/).

## Status

See [docs/status.md](docs/status.md) for what is done and what is next, and
[docs/roadmap.md](docs/roadmap.md) for the plan.

| Crate | Report |
|---|---|
| puruspe 0.4.4 | [reports/puruspe-0.4.4.md](reports/puruspe-0.4.4.md) |
| special 0.14.2 | [reports/special-0.14.2.md](reports/special-0.14.2.md) |
| statrs 0.19.1 | [reports/statrs-0.19.1.md](reports/statrs-0.19.1.md) |

## Running it

```sh
cargo test                                   # uses committed reference data; no Python needed
cd reference && uv sync && uv run python -m rsv_reference   # regenerate the data
```

The data format is described in [docs/data-format.md](docs/data-format.md).

## AI disclosure

The code, data generators and documentation in this repository are written by
Claude (an AI model by Anthropic) working under the direction of the project
owner, who reviews and merges every change. Findings are reported to upstream
projects only with the owner's explicit approval for each report, and those
reports disclose AI assistance.

## License

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or
[MIT license](LICENSE-MIT) at your option. Unless you explicitly state
otherwise, any contribution intentionally submitted for inclusion in this work
shall be dual licensed as above, without any additional terms or conditions.
