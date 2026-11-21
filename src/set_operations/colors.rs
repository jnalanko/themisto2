//! Color mappings of set operations: how the colors of the two input indexes become the colors of
//! the result, matched by name where needed.

use std::collections::{HashMap, HashSet};

pub(crate) fn check_unique_names(names: &[String], which: &str) -> Result<(), String> {
    let mut seen = HashSet::<&str>::with_capacity(names.len());
    for name in names {
        if !seen.insert(name.as_str()) {
            return Err(format!("Color name {:?} appears more than once in the {} index. Color names must be unique to merge shared colors.", name, which));
        }
    }
    Ok(())
}

/// Returns the merged color id of each color of the second coloring, and the color names of
/// the merged coloring. The colors of the first coloring keep their ids. If merge_shared is
/// true, a color of the second coloring with the same name as a color of the first coloring is
/// mapped to that color. All other colors of the second coloring get new ids after the colors
/// of the first coloring, in their original order. Returns an error if merge_shared is true and
/// a name appears twice within one coloring.
pub(crate) fn merged_color_mapping(names1: &[String], names2: &[String], merge_shared: bool) -> Result<(Vec<usize>, Vec<String>), String> {
    let n1 = names1.len();
    let mut merged_names = names1.to_vec();

    if !merge_shared {
        merged_names.extend(names2.iter().cloned());
        return Ok(((n1..n1 + names2.len()).collect(), merged_names));
    }

    check_unique_names(names1, "first")?;
    check_unique_names(names2, "second")?;

    let name_to_id1: HashMap<&str, usize> = names1.iter().enumerate().map(|(i, name)| (name.as_str(), i)).collect();
    let color2_to_merged: Vec<usize> = names2.iter().map(|name| {
        name_to_id1.get(name.as_str()).copied().unwrap_or_else(|| {
            merged_names.push(name.clone());
            merged_names.len() - 1
        })
    }).collect();

    log::info!("{} colors are shared between the two indexes", n1 + names2.len() - merged_names.len());
    Ok((color2_to_merged, merged_names))
}

/// Returns, for the intersection of the color sets, the result color id of each color of the
/// first and of the second coloring (None if the color is not in the result), and the color names
/// of the result. The result has the colors whose names are in both colorings, in the order of the
/// first coloring. Returns an error if a name appears twice within one coloring.
#[allow(clippy::type_complexity)]
pub(crate) fn intersected_color_mapping(names1: &[String], names2: &[String]) -> Result<(Vec<Option<usize>>, Vec<Option<usize>>, Vec<String>), String> {
    check_unique_names(names1, "first")?;
    check_unique_names(names2, "second")?;

    let in_names2: std::collections::HashSet<&str> = names2.iter().map(|name| name.as_str()).collect();
    let mut result_names = Vec::<String>::new();
    let mut name_to_result = HashMap::<&str, usize>::new();
    let color1_to_result: Vec<Option<usize>> = names1.iter().map(|name| {
        in_names2.contains(name.as_str()).then(|| {
            name_to_result.insert(name.as_str(), result_names.len());
            result_names.push(name.clone());
            result_names.len() - 1
        })
    }).collect();
    let color2_to_result: Vec<Option<usize>> = names2.iter().map(|name| name_to_result.get(name.as_str()).copied()).collect();

    log::info!("{} colors are shared between the two indexes", result_names.len());
    Ok((color1_to_result, color2_to_result, result_names))
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_merge_shared_colors_rejects_duplicate_names() {
        let a = vec!["A".to_string(), "A".to_string()];
        let b = vec!["B".to_string()];
        assert!(super::merged_color_mapping(&a, &b, true).is_err());
        assert!(super::merged_color_mapping(&b, &a, true).is_err());
        assert!(super::merged_color_mapping(&a, &b, false).is_ok()); // Names are not matched without the flag
    }

    #[test]
    fn test_merged_color_mapping() {
        let to_strings = |v: &[&str]| -> Vec<String> { v.iter().map(|s| s.to_string()).collect() };
        let names1 = to_strings(&["A", "B", "C"]);
        let names2 = to_strings(&["D", "C", "E", "A"]);

        let (map, names) = super::merged_color_mapping(&names1, &names2, true).unwrap();
        assert_eq!(map, vec![3, 2, 4, 0]);
        assert_eq!(names, to_strings(&["A", "B", "C", "D", "E"]));

        let (map, names) = super::merged_color_mapping(&names1, &names2, false).unwrap();
        assert_eq!(map, vec![3, 4, 5, 6]);
        assert_eq!(names, to_strings(&["A", "B", "C", "D", "C", "E", "A"]));
    }

    #[test]
    fn test_intersected_color_mapping() {
        let to_strings = |v: &[&str]| -> Vec<String> { v.iter().map(|s| s.to_string()).collect() };
        let names1 = to_strings(&["A", "B", "C", "D"]);
        let names2 = to_strings(&["D", "E", "A"]);
        let (map1, map2, names) = super::intersected_color_mapping(&names1, &names2).unwrap();
        assert_eq!(map1, vec![Some(0), None, None, Some(1)]);
        assert_eq!(map2, vec![Some(1), None, Some(0)]);
        assert_eq!(names, to_strings(&["A", "D"]));

        let duplicates = to_strings(&["A", "A"]);
        assert!(super::intersected_color_mapping(&duplicates, &names2).is_err());
        assert!(super::intersected_color_mapping(&names2, &duplicates).is_err());
    }
}
