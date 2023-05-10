use std::fmt;
use rand::Rng;
use primitive_types::U256;

#[derive(Clone, Copy)]
pub struct FieldElement {
    pub val: i128
}

/// A trait to represent zero of a numeric type.
pub trait Zero {
    fn zero() -> Self;
}

impl FieldElement {
    // Constants and constructors

    const K_MODULUS: i128 = 3 * 2i128.pow(30u32) + 1;
    const GENERATOR_VAL: i128 = 5i128;

    /// Creates a new `FieldElement` with the specified value.
    ///
    /// # Arguments
    ///
    /// * `val` - An `i128` value to initialize the field element.
    ///
    /// # Returns
    ///
    /// A new `FieldElement` instance with the normalized value.
    pub fn new(val: i128) -> Self {
        // Normalize the value
        let normalized_val = ((val % Self::K_MODULUS) + Self::K_MODULUS) % Self::K_MODULUS;
        FieldElement { val: normalized_val }
    }

    // Utility methods

    /// Returns a `FieldElement` representing zero.
    pub fn zero() -> FieldElement {
        FieldElement::new(0)
    }

    /// Returns a `FieldElement` representing one.
    pub fn one() -> FieldElement {
        FieldElement::new(1)
    }

    /// Returns a `FieldElement` representing the generator value.
    pub fn generator() -> FieldElement {
        FieldElement::new(Self::GENERATOR_VAL)
    }

    /// Returns the modulus of the field as an `i128`.
    pub fn modulus() -> i128 {
        Self::K_MODULUS
    }

    /// Returns a random field element
    pub fn random_element() -> FieldElement {
        let mut rng = rand::thread_rng();
        let random_integer = rng.gen_range(1..Self::K_MODULUS);
        FieldElement::new(random_integer)
    }

    /// Returns the multiplicative inverse of the field element.
    pub fn inverse(self) -> FieldElement {
        let (mut t, mut new_t) = (0i128, 1i128);
        let (mut r, mut new_r) = (Self::K_MODULUS, self.val.clone());

        while new_r != 0i128 {
            let quotient = r.div_euclid(new_r);
            (t, new_t) = (new_t, (t - (quotient * new_t)));
            (r, new_r) = (new_r, (r - (quotient * new_r)))
        }
        return FieldElement::new(t)

    }

    /// Raises the field element to the power of the given exponent.
    ///
    /// # Arguments
    ///
    /// * `exponent` - An `i128` exponent value.
    ///
    /// # Returns
    ///
    /// A new `FieldElement` instance with the result of the exponentiation.
    pub fn pow(&self, mut exponent: i128) -> FieldElement {
        let mut cur_pow = self.clone();
        let mut res = FieldElement::new(1);

        while exponent > 0i128 {
            if exponent % 2 != 0 {
                res = res * cur_pow;
            }
            exponent = exponent.div_euclid(2);
            cur_pow = cur_pow * cur_pow;
        }
        res
    }

    /// Checks if the field element has the specified order.
    ///
    /// # Arguments
    ///
    /// * `n` - An `i128` value representing the order to check.
    ///
    /// # Returns
    ///
    /// A boolean value indicating if the field element has the specified order.
    pub fn is_order(&self, n: i128) -> bool {
        let mut h = FieldElement::one();

        for _ in 1..n {
            h = h * *self;
            if h == FieldElement::one() {
                return false;
            }
        }
        return h * *self == FieldElement::one();
    }
}

// Implementations of standard traits

impl fmt::Display for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.val)
    }
}

impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        let normalized_self_val = ((self.val % Self::K_MODULUS) + Self::K_MODULUS) % Self::K_MODULUS;
        let normalized_other_val = ((other.val % Self::K_MODULUS) + Self::K_MODULUS) % Self::K_MODULUS;
        normalized_self_val == normalized_other_val
    }
}

impl std::ops::Neg for FieldElement {

    type Output = FieldElement;

    fn neg(self) -> FieldElement {
       FieldElement{
        val: FieldElement::zero().val - self.val,
       }
    }
}

impl std::ops::Add for FieldElement {
    type Output = FieldElement;

    fn add(self, other: FieldElement) -> FieldElement {
        FieldElement { val: (self.val + other.val) % Self::K_MODULUS }
    }
}

impl std::ops::Sub for FieldElement {
    type Output = FieldElement;

    fn sub(self, rhs: Self) -> Self::Output {
        FieldElement {
            val: (self.val - rhs.val) % Self::K_MODULUS
        }
    }
}

impl std::ops::Mul for FieldElement {
    type Output = FieldElement;

    fn mul(self, rhs: Self) -> Self::Output {
        FieldElement::new(self.val * rhs.val)
    }
}

impl std::ops::Div for FieldElement {
    type Output = FieldElement;

    fn div(self, rhs: Self) -> Self::Output {
        FieldElement::new(rhs.inverse().val * self.val)
    }
}

impl fmt::Debug for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", (self.val + Self::K_MODULUS.div_euclid(2)) % Self::K_MODULUS - Self::K_MODULUS.div_euclid(2))
    }
}

impl Zero for FieldElement {
    fn zero() -> FieldElement {
        Self::zero()
    }
}

impl From<i128> for FieldElement {
    fn from(value: i128) -> Self {
        // Assuming FieldElement has a constructor that takes i32
        FieldElement::new(value)
    }
}

impl From<U256> for FieldElement {
    fn from(value: U256) -> Self {
        // Assuming FieldElement has a constructor that takes i32
        let new_value = value % U256::from(Self::K_MODULUS);
        FieldElement::new(new_value.as_u128() as i128)
    }
}

impl Default for FieldElement {
    fn default() -> Self {
        FieldElement { val: (0) }
    }
}

#[cfg(test)]
mod field {
    use super::FieldElement;

    #[test]
    fn test_zero(){
        assert_eq!(FieldElement::zero(), FieldElement::new(0), "Should have been 0");
    }

    #[test]
    fn test_one(){
        assert_eq!(FieldElement::one(), FieldElement::new(1), "Should have been 1");
    }

    #[test]
    fn test_generator(){
        assert_eq!(FieldElement::generator().val, FieldElement::GENERATOR_VAL, "Wrong generator value");
    }

    #[test]
    fn test_modulus(){
        assert_eq!(FieldElement::modulus(), FieldElement::K_MODULUS, "Wrong modulus value");
        assert_eq!(FieldElement::K_MODULUS % FieldElement::modulus(), FieldElement::zero().val, "Wrong modulus value");
    }

    #[test]
    fn test_random_element(){
        assert!(FieldElement::random_element().val < FieldElement::K_MODULUS, "Wrong random element range");
    }

    #[test]
    fn test_inverse(){
        let r = FieldElement::random_element();
        assert_eq!(((r * r.inverse()).val + FieldElement::modulus()) % FieldElement::K_MODULUS, 1, "Multiplication of an element with its inverse should produce 1");
    }

    #[test]
    fn test_pow(){
        assert_eq!(FieldElement::generator().pow(FieldElement::modulus() - 1), FieldElement::one(), "It should have been 1");
    }

    #[test]
    fn test_is_order(){
        let g = FieldElement::generator().pow(3 * 2i128.pow(20));
        assert_eq!(g.is_order(1024), true, "Invalid order");
    }
}