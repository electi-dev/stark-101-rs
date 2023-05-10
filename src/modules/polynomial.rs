// use std::vec;
use crate::{
    list_utils,
    FieldElement};

/// A struct representing a polynomial.
#[derive(Clone, Debug)]
pub struct Polynomial {
    pub poly: Vec<FieldElement>,
    pub var: char,
}

impl Polynomial {

    /// Create a new polynomial with the given coefficients and variable.
    /// 
    /// # Arguments
    ///
    /// * `coeffs` - A vector of coefficients for the polynomial
    /// * `var` - A character representing the variable of the polynomial
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance
    pub fn new(coeffs: &Vec<FieldElement>, var: char) -> Self {
        let modified_coeffs = list_utils::remove_trailing_elements(&coeffs, &FieldElement::zero());
        Polynomial { poly: modified_coeffs, var: var }
    }

    /// Return a LaTeX representation of the polynomial.
    ///
    /// # Returns
    ///
    /// A string containing the LaTeX representation of the polynomial
    pub fn repr_latex(&self) -> String {
        if self.poly.is_empty() {
            return String::from("$0$");
        }
    
        let mut res = String::from("$");
        let mut first = true;
    
        for (exponent, coef) in self.poly.iter().enumerate() {
            if coef == &FieldElement::zero() {
                continue;
            }
    
            let mut monomial = Polynomial::latex_monomial(exponent, coef, self.var);
            if first {
                first = false;
                res.push_str(&monomial);
                continue;
            }
    
            let oper = if monomial.starts_with('-') {
                monomial.remove(0); // Remove the '-' sign
                '-'
            } else {
                '+'
            };
    
            let formatted_monomial = format!(" {}{}", oper, monomial);
            res.push_str(&formatted_monomial);
        }
    
        res.push('$');
        res
    }

    /// Create a LaTeX representation of a monomial given an exponent, coefficient, and variable.
    ///
    /// # Arguments
    ///
    /// * `exponent` - The exponent of the monomial
    /// * `coef` - The coefficient of the monomial
    /// * `var` - The variable of the monomial
    ///
    /// # Returns
    ///
    /// A string containing the LaTeX representation of the monomial
    pub fn latex_monomial(exponent: usize, coef: &FieldElement, var: char) -> String {
        if exponent == 0 {
            return coef.to_string();
        }
        let coef_str = if coef == &FieldElement::one() {
            String::new()
        } else if coef == &FieldElement::new(-1) {
            String::from("-")
        } else {
            coef.to_string()
        };
        if exponent == 1 {
            format!("{}{}", coef_str, var)
        } else {
            format!("{}{}^{{{}}}", coef_str, var, exponent)
        }
    }

    /// Return the degree of the polynomial.
    ///
    /// # Returns
    ///
    /// An usize representing the degree of the polynomial
    pub fn degree(&self) -> usize {
        if self.poly.len() == 0 {
            return 0;
        } else {
            return list_utils::remove_trailing_elements(&self.poly, &FieldElement::zero()).len() - 1;
        }
    }      

    /// Compose two polynomials.
    ///
    /// # Arguments
    ///
    /// * `other` - A reference to another `Polynomial` instance
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the composition of the two polynomials
    pub fn compose(&self, other: &Polynomial) -> Polynomial {
        let mut res = Polynomial::new(&vec![], self.var);
    
        for coef in self.poly.iter().rev() {
            let coef_poly = Polynomial::new(&vec![coef.clone()], self.var);
            res = (res * other.clone()) + coef_poly;
        }
    
        res
    }
    
    /// Divide two polynomials and return the quotient and remainder as a tuple.
    ///
    /// # Arguments
    ///
    /// * `other` - A reference to another `Polynomial` instance representing the divisor
    ///
    /// # Returns
    ///
    /// A tuple containing the quotient and remainder as `Polynomial` instances
    pub fn qdiv(&self, other: &Polynomial) -> (Polynomial, Polynomial) {
        let pol2 = other.poly.clone();
        assert!(!pol2.is_empty(), "Dividing by zero polynomial.");

        let pol1 = self.poly.clone();
        if pol1.is_empty() {
            return (Polynomial::new(&vec![], self.var), Polynomial::new(&vec![], self.var));
        }

        let mut rem: Vec<FieldElement> = pol1.clone();
        let mut deg_dif: i128 = (rem.len() - pol2.len()) as i128;
        let mut quotient = vec![FieldElement::zero() ; (deg_dif + 1) as usize];

        let g_msc_inv = pol2.iter().last().unwrap().inverse();

        while deg_dif >= 0 {
            let tmp = *rem.last().unwrap() * g_msc_inv;
            quotient[deg_dif as usize] = quotient[deg_dif as usize] + tmp;
    
            let mut last_non_zero = deg_dif - 1;
            for (i, coef) in pol2.iter().enumerate().map(|(index, coef)| (index + deg_dif as usize, coef)) {//.skip(deg_dif as usize) {
                rem[i] = rem[i] - (tmp * *coef);
                if rem[i] != FieldElement::zero() {
                    last_non_zero = i as i128;
                }
            }
    
            // Eliminate trailing zeroes (i.e. make r end with its last non-zero coefficient).
            rem.truncate((last_non_zero + 1) as usize);
            deg_dif = rem.len() as i128 - pol2.len() as i128;
        }
    
        (
            Polynomial::new(&quotient, self.var),
            Polynomial::new(&rem, self.var),
        )

    }

    /// Get the coefficient of the nth degree term of the polynomial.
    ///
    /// # Arguments
    ///
    /// * `n` - The degree of the term whose coefficient is to be returned
    ///
    /// # Returns
    ///
    /// A `FieldElement` representing the coefficient of the nth degree term
    pub fn get_nth_degree_coefficient(&self, n: usize) -> FieldElement {
        match n > self.degree() {
            true => FieldElement::zero(),
            false => self.poly[n]
        }
    }

    /// Multiply the polynomial by a scalar.
    ///
    /// # Arguments
    ///
    /// * `scalar` - A reference to a `FieldElement` representing the scalar
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the result of the scalar multiplication
    pub fn scalar_mul(&self, scalar: &FieldElement) -> Polynomial {
        Polynomial { poly: list_utils::scalar_multiplication::<FieldElement>(&self.poly, scalar), var: 'x' }
    }

    /// Evaluate the polynomial at the given point.
    ///
    /// # Arguments
    ///
    /// * `point` - A reference to a `FieldElement` representing the point at which the polynomial is to be evaluated
    ///
    /// # Returns
    ///
    /// A `FieldElement` representing the result of the evaluation
    pub fn eval(&self, point: &FieldElement) -> FieldElement {
        self.poly
            .iter()
            .rev()
            .fold(FieldElement::zero(), |acc, coeff| (acc * *point + *coeff))
    }

    /// Divide two polynomials and return the result as a new `Polynomial` instance. 
    /// The remainder should be zero, otherwise, it will panic.
    ///
    /// # Arguments
    ///
    /// * `other` - A reference to another `Polynomial` instance representing the divisor
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the result of the division
    pub fn truediv(&self, other: &Polynomial) -> Polynomial {
        let (div, modulo): (Polynomial, Polynomial) = self.qdiv(&other);
        assert!(modulo == Polynomial::new(&vec![FieldElement::zero()], self.var));
        div
    }

    /// Calculate the remainder of the division of two polynomials.
    ///
    /// # Arguments
    ///
    /// * `other` - A reference to another `Polynomial` instance representing the divisor
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the remainder of the division
    pub fn modulo(&self, other: &Polynomial) -> Polynomial {
        self.qdiv(&other).1
    }

    /// Create a monomial with the given degree and coefficient.
    ///
    /// # Arguments
    ///
    /// * `degree` - The degree of the monomial
    /// * `coeff` - A reference to a `FieldElement` representing the coefficient
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the monomial
    pub fn monomial(degree: i128, coeff: &FieldElement) -> Polynomial {
        Polynomial { poly: vec![FieldElement::zero(); degree as usize].into_iter().chain(vec![*coeff]).collect(), var: 'x' }
    }

    /// Generate a linear term of the form x - point.
    ///
    /// # Arguments
    ///
    /// * `point` - A reference to a `FieldElement` representing the point
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the linear term
    pub fn gen_linear_term(point: &FieldElement) -> Polynomial {
        Polynomial { poly: vec![FieldElement::zero() - point.clone(), FieldElement::one()], var: 'x' }
    }

    /// Raise the polynomial to a given power.
    ///
    /// # Arguments
    ///
    /// * `other` - The exponent as an i128 value
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the result of the exponentiation
    pub fn pow(&self, other: i128) -> Polynomial {
        let mut res = Polynomial::new(&vec![FieldElement::one()], 'x');
        let mut cur = self.clone();
        let mut exponent = other.clone();
        loop {
            if exponent % 2 != 0 {
                res = res * cur.clone();
            }
            exponent = exponent >> 1;
            if exponent == 0 {
                break;
            }
            cur = cur.clone() * cur;
        }
        res
    }
    
    /// Calculate the Lagrange polynomials for the given x_values.
    ///
    /// # Arguments
    ///
    /// * `x_values` - A reference to a vector of `FieldElement` instances representing the x_values
    ///
    /// # Returns
    ///
    /// A vector of `Polynomial` instances representing the Lagrange polynomials
    pub fn calculate_lagrange_polynomials(x_values: &Vec<FieldElement>) -> Vec<Polynomial> {
        let monomials: Vec<Polynomial> = x_values
            .iter()
            .map(|x| {
                Polynomial::monomial(1, &FieldElement::one()) - Polynomial::monomial(0, x)
            })
            .collect();
    
        let numerator = list_utils::prod(&monomials);
    
        let mut lagrange_polynomials = Vec::with_capacity(x_values.len());
    
        for j in 0..x_values.len() {
            let denominator = list_utils::prod(
                &x_values
                .iter()
                .enumerate()
                .filter(|&(i, _)| i != j)
                .map(|(_, x)| x_values[j].clone() - x.clone())
                .collect::<Vec<FieldElement>>());
    
            let (cur_poly, _) = numerator.qdiv(
                &monomials[j]
                    .scalar_mul(&denominator),
            );
    
            lagrange_polynomials.push(cur_poly);
        }
    
        lagrange_polynomials
    }
    
    /// Interpolate a polynomial using the given y_values and Lagrange polynomials.
    ///
    /// # Arguments
    ///
    /// * `y_values` - A reference to a vector of `FieldElement` instances representing the y_values
    /// * `lagrange_polynomials` - A reference to a vector of `Polynomial` instances representing the Lagrange polynomials
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the interpolated polynomial
    pub fn interpolate_poly_lagrange(y_values: &Vec<FieldElement>, lagrange_polynomials: &Vec<Polynomial>) -> Polynomial {
        let mut poly = Polynomial::new(&Vec::new(), 'x');
    
        for (j, y_value) in y_values.iter().enumerate() {
            poly = poly + lagrange_polynomials[j].scalar_mul(y_value);
        }
    
        poly
    }

    /// Interpolate a polynomial using the given x_values and y_values.
    ///
    /// # Arguments
    ///
    /// * `x_values` - A reference to a vector of `FieldElement` instances representing the x_values
    /// * `y_values` - A reference to a vector of `FieldElement` instances representing the y_values
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the interpolated polynomial
    pub fn interpolate_poly(x_values: &Vec<FieldElement>, y_values: &Vec<FieldElement>) -> Polynomial {
        assert_eq!(x_values.len(), y_values.len());
        let lagrange_poly = Self::calculate_lagrange_polynomials(&x_values);
        let poly = Self::interpolate_poly_lagrange(y_values, &lagrange_poly);
        let ret = list_utils::remove_trailing_elements(&poly.poly, &FieldElement::zero());
        
        Polynomial { poly: ret, var: 'x' }
    }
    
}

impl PartialEq for Polynomial {
    ///
    /// # Arguments
    ///
    /// * `other` - A reference to another `Polynomial` instance
    ///
    /// # Returns
    ///
    /// `true` if the two polynomials are equal, otherwise `false`
    fn eq(&self, other: &Self) -> bool {
        self.poly == other.poly
    }
}


impl std::ops::Add for Polynomial {
    type Output = Polynomial;

    /// Add two polynomials.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The right-hand side polynomial
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the sum of the two polynomials
    fn add(self, rhs: Self) -> Self::Output {
        Polynomial {
            poly: list_utils::two_lists_tuple_operation(&self.poly, &rhs.poly, |a, b| a.clone() + b.clone()), var: self.var
        }
    }
}

impl std::ops::Sub for Polynomial {
    type Output = Polynomial;

    /// Subtract one polynomial from another.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The right-hand side polynomial
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the difference of the two polynomials
    fn sub(self, rhs: Self) -> Self::Output {
        Polynomial {
            poly: list_utils::two_lists_tuple_operation(&self.poly, &rhs.poly, |a, b| a.clone() - b.clone()), var: self.var
        }
    }
}

impl std::ops::Neg for Polynomial {
    type Output = Polynomial;

    /// Negate a polynomial.
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the negation of the polynomial
    fn neg(self) -> Self::Output {
        Polynomial::new(&vec![FieldElement::new(0)], self.var) - self
    }
}

impl std::ops::Mul for Polynomial {
    type Output = Polynomial;

    /// Multiply two polynomials.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The right-hand side polynomial
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the product of the two polynomials
    fn mul(self, rhs: Self) -> Self::Output {
        let res_len = self.poly.len() + rhs.poly.len() - 1;
        let mut res = vec![FieldElement::zero(); res_len];

        for (i, c1) in self.poly.iter().enumerate() {
            for (j, c2) in rhs.poly.iter().enumerate() {
                res[i + j] = res[i + j] + (*c1 * *c2);
            }
        }
        let res: Vec<FieldElement> = res.into_iter().map(FieldElement::from).collect();

        Polynomial { poly: res, var: self.var }
    }
}

impl std::ops::Div for Polynomial {
    type Output = Polynomial;

    /// Divide two polynomials.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The right-hand side polynomial
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the quotient of the two polynomials
    fn div(self, rhs: Self) -> Self::Output {
        let (quotient, remainder) = self.qdiv(&rhs);
        assert!(remainder != Polynomial::from(0));
        quotient
    }
}

impl From<i128> for Polynomial {
    /// Create a `Polynomial` instance from an `i128` value.
    ///
    /// # Arguments
    ///
    /// * `value` - An `i128` value
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the given value
    fn from(value: i128) -> Self {
        Polynomial { poly: vec![FieldElement::from(value)], var: 'x' }
    }
}

impl From<FieldElement> for Polynomial {
    /// Create a `Polynomial` instance from a `FieldElement` value.
    ///
    /// # Arguments
    ///
    /// * `value` - A `FieldElement` value
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the given value
    fn from(value: FieldElement) -> Self {
        Polynomial { poly: vec![value], var: 'x' }
    }
}

impl Default for Polynomial {
    /// Create a default `Polynomial` instance with a zero polynomial.
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the zero polynomial
    fn default() -> Self {
        Self { poly: vec![FieldElement::zero()], var: 'x' }
    }
}

impl std::ops::Rem for Polynomial {
    type Output = Polynomial;

    /// Calculate the remainder of the division of two polynomials.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The right-hand side polynomial
    ///
    /// # Returns
    ///
    /// A new `Polynomial` instance representing the remainder of the division
    fn rem(self, rhs: Self) -> Self::Output {
        return self.qdiv(&rhs).1.clone();
    }
}

#[cfg(test)]
mod polynomial {
    use super::*;

    #[test]
    fn test_new() {
        let coeffs = vec![
            FieldElement::new(0),
            FieldElement::new(0),
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(1),
            FieldElement::new(0),
            FieldElement::new(0),
            FieldElement::new(0),];
        
        let polynomial = Polynomial::new(&coeffs, 'x');
        assert_eq!(polynomial.poly, 
            vec![FieldElement::new(0), 
            FieldElement::new(0),
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(1)])
    }

    #[test]
    fn test_latex_repr() {
        let coeffs = vec![
            FieldElement::new(0),
            FieldElement::new(0),
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(1)];
        
        let polynomial = Polynomial::new(&coeffs, 'x');
        assert_eq!(polynomial.repr_latex(), "$x^{2} +2x^{3} +3x^{4} +x^{5}$");
    }

    #[test]
    fn test_add() {
        let coeffs = vec![
            FieldElement::new(0),
            FieldElement::new(0),
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(1)];
        
        let polynomial_1 = Polynomial::new(&coeffs, 'x');
        let polynomial_2 = Polynomial::new(&coeffs, 'x');
        assert_eq!((polynomial_1.clone() + polynomial_2).repr_latex(), "$2x^{2} +4x^{3} +6x^{4} +2x^{5}$");
    }

    #[test]
    fn test_sub() {
        let coeffs = vec![
            FieldElement::new(0),
            FieldElement::new(0),
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(1)];
        
        let polynomial_1 = Polynomial::new(&coeffs, 'x');
        let polynomial_2 = Polynomial::new(&coeffs, 'x');
        assert_eq!((polynomial_1.clone() - polynomial_2).repr_latex(), "$$");

        let coeffs2 = vec![
            FieldElement::one()
        ];

        let polynomial_3 = Polynomial::new(&coeffs2, 'x'); 
        assert_eq!((polynomial_1 - polynomial_3).repr_latex(), "$-1 +x^{2} +2x^{3} +3x^{4} +x^{5}$")
    }

    #[test]
    fn test_neg() {
        let coeffs = vec![
            FieldElement::new(0),
            FieldElement::new(0),
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(1)];
        
        let polynomial = Polynomial::new(&coeffs, 'y');
        assert_eq!((-polynomial).repr_latex(), "$-y^{2} -2y^{3} -3y^{4} -y^{5}$");
    }

    #[test]
    fn test_mul() {
        let poly1 = Polynomial::new(&vec![FieldElement::from(1), FieldElement::from(2), FieldElement::from(1)], 'x');
        let poly2 = Polynomial::new(&vec![FieldElement::from(1), FieldElement::from(1)], 'x');

        let result = poly1 * poly2;

        let expected = Polynomial::new(&vec![FieldElement::from(1), FieldElement::from(3), FieldElement::from(3), FieldElement::from(1)], 'x');
        assert_eq!(result, expected);
    }

    #[test]
    fn test_compose() {
        let f = Polynomial::new(&vec![FieldElement::new(1), FieldElement::new(1), FieldElement::new(1)], 'x');
        let g = Polynomial::new(&vec![FieldElement::new(1), FieldElement::new(1)], 'x');
        assert_eq!(f.compose(&g).repr_latex(), "$3 +3x +x^{2}$");
    }

    #[test]
    fn test_qdiv() {
        let poly1 = Polynomial::new(&vec![
            FieldElement::new(6),
            FieldElement::new(5),
            FieldElement::new(4),
            FieldElement::new(3),
            FieldElement::new(0),
            FieldElement::new(0),
        ], 'x');

        let poly2 = Polynomial::new(&vec![
            FieldElement::new(1),
            FieldElement::new(2),
        ], 'x');

        let expected_quotient = Polynomial::new(&vec![
            FieldElement::new(402653186),
            FieldElement::new(-805306367),
            FieldElement::new(-1610612735),
        ], 'x');

        let expected_remainder = Polynomial::new(&vec![
            FieldElement::new(-402653180),
        ], 'x');

        let (quotient, remainder) = poly1.qdiv(&poly2);

        // Compare the calculated quotient and remainder with the expected values
        assert_eq!(quotient, expected_quotient);
        assert_eq!(remainder, expected_remainder);
    }

    #[test]
    fn test_monomial() {
        let monomial = Polynomial::monomial(10, &FieldElement::new(3));
        assert_eq!(monomial.repr_latex(), "$3x^{10}$")
    }

    #[test]
    fn test_gen_linear_term() {
        let point = FieldElement::new(3);
        let linear_term = Polynomial::gen_linear_term(&point);

        let expected_coeffs = vec![FieldElement::new(-3), FieldElement::one()];
        let expected_poly = Polynomial::new(&expected_coeffs, 'x');

        assert_eq!(linear_term, expected_poly);
    }

    #[test]
    fn test_get_nth_degree_coefficient() {
        let poly = Polynomial::new(&vec![
            FieldElement::from(3),
            FieldElement::from(2),
            FieldElement::from(1),
        ], 'x');
        
        assert_eq!(poly.get_nth_degree_coefficient(0), FieldElement::from(3));
        assert_eq!(poly.get_nth_degree_coefficient(1), FieldElement::from(2));
        assert_eq!(poly.get_nth_degree_coefficient(2), FieldElement::from(1));
        assert_eq!(poly.get_nth_degree_coefficient(3), FieldElement::zero());
    }

    #[test]
    fn test_polynomial_scalar_mul() {
        let poly = Polynomial::new(&vec![FieldElement::from(1), FieldElement::from(2), FieldElement::from(1)], 'x');
        let scalar = FieldElement::from(2);

        let result = poly.scalar_mul(&scalar);

        let expected = Polynomial::new(&vec![FieldElement::from(2), FieldElement::from(4), FieldElement::from(2)], 'x');
        assert_eq!(result, expected);
    }

    #[test]
    fn test_eval() {
        let poly = Polynomial {
            poly: vec![
                FieldElement::from(5), 
                FieldElement::from(3), 
                FieldElement::from(2)],
            var: 'x',
        }; //5 + 3x + 2x^2
        
        let point = FieldElement::from(3);
        let result = poly.eval(&point);

        assert_eq!(result, FieldElement::from(32));
    }

    #[test]
    fn test_calculate_lagrange_polynomials() {
        let x_values = vec![
            FieldElement::from(4),
            FieldElement::from(3),
            FieldElement::from(2),
            FieldElement::from(1),
        ];

        let lagrange_polynomials = Polynomial::calculate_lagrange_polynomials(&x_values);

        // Expected lagrange polynomials
        let expected = vec![
            Polynomial::new(&vec![FieldElement::from(-1), FieldElement::from(536870914), FieldElement::from(-1), FieldElement::from(-536870912),], 'x'),
            Polynomial::new(&vec![FieldElement::from(4), FieldElement::from(-7), FieldElement::from(-1610612733), FieldElement::from(1610612736),], 'x'),
            Polynomial::new(&vec![FieldElement::from(-6), FieldElement::from(-1610612727), FieldElement::from(-4), FieldElement::from(-1610612736),], 'x'),
            Polynomial::new(&vec![FieldElement::from(4), FieldElement::from(1073741820), FieldElement::from(-1610612735), FieldElement::from(536870912),], 'x'),
            ];

        assert_eq!(lagrange_polynomials, expected);
    }

    #[test]
    fn test_interpolate_poly_lagrange() {
        let x_values = vec![
            FieldElement::from(1),
            FieldElement::from(2),
            FieldElement::from(3),
            FieldElement::from(4),
        ];
        let y_values = vec![
            FieldElement::from(2),
            FieldElement::from(5),
            FieldElement::from(10),
            FieldElement::from(20),
        ];

        let lagrange_polynomials = Polynomial::calculate_lagrange_polynomials(&x_values);
        let interpolated_poly = Polynomial::interpolate_poly_lagrange(&y_values, &lagrange_polynomials);

        let expected_poly = Polynomial::new(&vec![
            FieldElement::from(-2),
            FieldElement::from(-1610612731),
            FieldElement::from(-2),
            FieldElement::from(-1610612736),
        ], 'x');

        assert_eq!(interpolated_poly, expected_poly);

    }
    

    #[test]
    fn test_interpolate_poly() {
        let x_values = vec![
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(4),
        ];

        let original_poly = Polynomial::new(&vec![
            FieldElement::new(5),
            FieldElement::new(3),
            FieldElement::new(2),
        ], 'x');

        // Calculate the corresponding y values for the given x values
        let y_values: Vec<FieldElement> = x_values
            .iter()
            .map(|x| original_poly.eval(x))
            .collect();

        let interpolated_poly = Polynomial::interpolate_poly(&x_values, &y_values);

        // Compare the interpolated polynomial with the original polynomial
        assert_eq!(interpolated_poly, original_poly);
    }

    #[test]
    fn test_interpolate_poly2() {
        let x_values = vec![
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(4),
        ];

        let y_values = vec![
            FieldElement::new(5),
            FieldElement::new(6),
            FieldElement::new(7),
            FieldElement::new(21),
        ];

        let interpolated_poly = Polynomial::interpolate_poly(&x_values, &y_values);

        // Compare the interpolated polynomial with the original polynomial
        assert_eq!(interpolated_poly.repr_latex(), "$3221225464 +536870937x +3221225460x^{2} +2684354563x^{3}$");
    }
}
