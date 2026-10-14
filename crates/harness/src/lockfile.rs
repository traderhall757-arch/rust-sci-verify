//! Finding a dependency's exact version in `Cargo.lock`.

/// The single locked version of package `name`, or `None` if it is absent
/// or locked at more than one version (ambiguous).
#[must_use]
pub fn locked_version(lock: &str, name: &str) -> Option<String> {
    let wanted = format!("name = \"{name}\"");
    let mut found = None;
    let mut lines = lock.lines();
    while let Some(line) = lines.next() {
        if line.trim() != wanted {
            continue;
        }
        let version = lines
            .next()?
            .trim()
            .strip_prefix("version = \"")?
            .strip_suffix('"')?
            .to_string();
        if found.replace(version).is_some() {
            return None;
        }
    }
    found
}

#[cfg(test)]
mod tests {
    use super::locked_version;

    const LOCK: &str = "[[package]]\nname = \"puruspe\"\nversion = \"0.4.4\"\n\n[[package]]\nname = \"dup\"\nversion = \"1.0.0\"\n\n[[package]]\nname = \"dup\"\nversion = \"2.0.0\"\n\n[[package]]\nname = \"broken\"\nsource = \"x\"\n";

    #[test]
    fn finds_single_versions_only() {
        assert_eq!(locked_version(LOCK, "puruspe").as_deref(), Some("0.4.4"));
        assert_eq!(locked_version(LOCK, "dup"), None);
        assert_eq!(locked_version(LOCK, "absent"), None);
        assert_eq!(locked_version(LOCK, "broken"), None);
        assert_eq!(locked_version("name = \"x\"", "x"), None);
    }
}
