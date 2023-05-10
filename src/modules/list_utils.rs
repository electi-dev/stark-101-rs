use super::field::Zero;
use std::ops::Mul;

/// Removes trailing occurrences of a specified element from a given list.
///
/// # Arguments
///
/// * `list_of_elements` - A reference to a vector of elements of type `F`.
/// * `element_to_remove` - A reference to the element of type `F` to be removed from the end of the list.
///
/// # Returns
///
/// A new vector with the trailing occurrences of the specified element removed.
pub fn remove_trailing_elements<F: PartialEq + Clone>(list_of_elements: &Vec<F>, element_to_remove: &F) -> Vec<F> {
    if let Some(last_element) = list_of_elements.last() {
        if last_element == element_to_remove {
            let elements_without_trailing = list_of_elements
                .iter()
                .rev()
                .skip_while(|x| x == &element_to_remove)
                .cloned()
                .collect::<Vec<F>>()
                .into_iter()
                .rev()
                .collect();
            return elements_without_trailing;
        }
    }
    list_of_elements.clone()
}

/// Applies a specified function to corresponding elements of two lists, padding the shorter list with zeros.
///
/// # Arguments
///
/// * `f` - A reference to the first vector of elements of type `T`.
/// * `g` - A reference to the second vector of elements of type `T`.
/// * `func` - A function of type `F` that takes two references to elements of type `T` and returns a value of type `T`.
///
/// # Returns
///
/// A new vector with the results of applying the specified function to corresponding elements of the input vectors.
pub fn two_lists_tuple_operation<F, T: Clone + Zero>(f: &Vec<T>, g: &Vec<T>, func: F) -> Vec<T> 
    where
        F: Fn(&T, &T) -> T
{

    let mut padded_f = f.clone();
    let mut padded_g = g.clone();

    if f.len() < g.len() {
        padded_f.extend(vec![T::zero(); g.len() - f.len()]);
    } else if g.len() < f.len() {
        padded_g.extend(vec![T::zero(); f.len() - g.len()]);
    }

    padded_f
        .iter()
        .zip(padded_g.iter())
        .map(|(elem1, elem2)| func(elem1, elem2))
        .collect()
}

/// Multiplies each element of a given list by a specified scalar.
///
/// # Arguments
///
/// * `list_of_elements` - A reference to a vector of elements of type `T`.
/// * `scalar` - A reference to the scalar of type `T` to multiply each element by.
///
/// # Returns
///
/// A new vector with the results of multiplying each element of the input vector by the specified scalar.
pub fn scalar_multiplication<T: Mul<Output = T> + Clone>(list_of_elements: &Vec<T>, scalar: &T) -> Vec<T> {
    list_of_elements
        .iter()
        .map(|elem| elem.clone() * scalar.clone())
        .collect()
}

/// Calculates the product of all elements in a given slice.
///
/// # Arguments
///
/// * `values` - A slice of elements of type `T`.
///
/// # Returns
///
/// The product of all elements in the slice. If the slice is empty, the function returns the default value of type `T`.
pub fn prod<T>(values: &[T]) -> T
where
    T: Mul<Output = T> + Clone + Default,
{
    let len_values = values.len();
    if len_values == 0 {
        return T::default();
    }
    if len_values == 1 {
        return values[0].clone();
    }
    let mid = len_values / 2;
    prod(&values[..mid]) * prod(&values[mid..])
}

#[cfg(test)]
mod list_utils_test {
    use crate::FieldElement;

    #[test]
    fn test_remove_trailing_elements(){
        let a = vec![
            FieldElement::new(0),
            FieldElement::new(0),
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(1),
            FieldElement::new(0),
            FieldElement::new(0),
            FieldElement::new(0),];

        super::remove_trailing_elements(&a, &FieldElement::new(0));

        assert_eq!(a.last(), Some(&FieldElement::zero()));
    }

    #[test]
    fn test_two_lists_tuple_operation() {
        let vec1 = vec![FieldElement::new(1), FieldElement::new(2), FieldElement::new(3)];
        let vec2 = vec![FieldElement::new(3), FieldElement::new(4)];
        let res = super::two_lists_tuple_operation(&vec1, &vec2, |a, b| a.clone() + b.clone());
        assert_eq!(vec![FieldElement::new(4), FieldElement::new(6), FieldElement::new(3)], res);
    }

    #[test]
    fn test_scalar_multiplication() {
        let vec1 = vec![FieldElement::new(1), FieldElement::new(2), FieldElement::new(3)];
        let scalar = FieldElement::new(6);
        let res = super::scalar_multiplication(&vec1, &scalar);
        assert_eq!(vec![FieldElement::new(6), FieldElement::new(12), FieldElement::new(18)], res);
    }
}
