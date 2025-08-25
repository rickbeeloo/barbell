/// Get edit distance based on formula presented in paper
pub fn get_edit_cut_off(l: usize) -> usize {
    let a = l as f64;
    let value = 0.5100_f64 * a - 1.7312_f64 * a.sqrt();
    let ceil_value = value.ceil();
    if ceil_value > 0.0 {
        ceil_value as usize
    } else {
        0
    }
}
