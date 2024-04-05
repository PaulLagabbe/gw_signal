/* --------------------------------------------------------------------------------------------- *
 * Frequency series definition
 * --------------------------------------------------------------------------------------------- */

use num::complex::{
	Complex,
	ComplexFloat
};
use std::ops::{Add, Neg, Sub, Mul, Div, Index, IndexMut};
use more_asserts::assert_gt;

/* --------------------------------------------------------------------------------------------- *
 * Define structures
 * --------------------------------------------------------------------------------------------- */
/// Frequency series object
/// The data vector is a complex value vector indexed by the frequency
/// 
/// The frequency vector start from 0 Hz.
#[derive(Debug, Clone)]
pub struct Specrogram {
    f_max: f64,
	dt: f64,
    data: Vec<Vec<Complex<f64>>>,
}



/* --------------------------------------------------------------------------------------------- *
 * Constructors
 * --------------------------------------------------------------------------------------------- */

impl FrequencySeries {
	/// Constructor
	pub fn from_vector(f_max: f64, input_data: Vec<Vec<Complex<f64>>>) -> Self {

        FrequencySeries {
            f_max,
            data: input_data
        }
    }
	
	/// get data size
	pub fn get_size(&self) -> (usize, usize) {
		(self.data.len(), self.data[0].len())
	}
	/// get maximum frequency
	pub fn get_f_max(&self) -> f64 {
		self.f_max
	}
	/// get time step
	pub fn get_f_max(&self) -> f64 {
		self.dt
	}
	/// get data vector
	pub fn get_data(&self) -> Vec<Vec<Complex<f64>>> {
		self.data.clone()
	}
}
/*

	/// compute real part of the data vector
	pub fn real(mut self) -> FrequencySeries {
		for i in 0..self.data.len() {
			self.data[i] = Complex{re: self.data[i].re, im: 0.0};
		}
		self
	}
	/// compute imaginary part of the data vector
	pub fn imag(mut self) -> FrequencySeries {
		for i in 0..self.data.len() {
			self.data[i] = Complex{re: self.data[i].im, im: 0.0};
		}
		self
	}
	/// compute absolute value of the data vector
	pub fn abs(mut self) -> FrequencySeries {
		for i in 0..self.data.len() {
			self.data[i] = Complex{re: self.data[i].norm(), im: 0.0};
		}
		self
	}
	/// compute square of absolute value of the data vector
	pub fn abs2(mut self) -> FrequencySeries {
		for i in 0..self.data.len() {
			self.data[i] = self.data[i] * self.data[i].conj();
		}
		self
	}
	/// compute phase of the data vector
	pub fn arg(mut self) -> FrequencySeries {
		for i in 0..self.data.len() {
			self.data[i] = Complex{re: self.data[i].arg(), im: 0.0};
		}
		self
	}
	/// compute complex conjugate of the data vector
	pub fn conj(mut self) -> FrequencySeries {
		
		for i in 0..self.data.len() {
			self.data[i] = self.data[i].conj();
		}
		self
	}
	/// compute square root of the data vector
	pub fn sqrt(mut self) -> FrequencySeries {
		
		for i in 0..self.data.len() {
			self.data[i] = self.data[i].sqrt();
		}
		self
	}
	/// compute inverse of the data vector
	pub fn inv(&mut self) -> &mut FrequencySeries {
		
		for i in 0..self.data.len() {
			self.data[i] = self.data[i].inv();
		}
		self
	}
	pub fn max_abs(&self) -> f64 {
		
		let mut maximum: f64 = 0f64;
		for i in 0..self.get_size() {
			if maximum < self[i].abs() {
				maximum = self[i].abs();
			}
		}
		maximum
	}
	pub fn min_abs(&self) -> f64 {
		
		let mut minimum: f64 = f64::INFINITY;
		for i in 0..self.get_size() {
			if minimum > self[i].abs() {
				minimum = self[i].abs();
			}
		}
		minimum
	}
}

/* --------------------------------------------------------------------------------------------- *
 * Operator overloading:
 * notes:
 *      The operators DOES NOT create a new time series, they modify one of the time series
 *      parameter.
 * --------------------------------------------------------------------------------------------- */


/* operator+ ----------------------------------------------------------------------------------- */
impl<'a> Add<&FrequencySeries> for &'a mut FrequencySeries {

    type Output = &'a mut FrequencySeries;

    fn add(self, other: &FrequencySeries) -> &'a mut FrequencySeries {
        // verify if the length of the time series are equal
        assert_eq!(self.data.len(), other.data.len());
        // Modify data vector of the left time series
        for i in 0..self.data.len() {
            self.data[i] += other.data[i];
        }
        // return modified left timeseries
        self
    }

}
impl<'a> Add<Complex<f64>> for &'a mut FrequencySeries {

    type Output = &'a mut FrequencySeries;

    fn add(self, other: Complex<f64>) -> &'a mut FrequencySeries {
        // Modify data vector of the left time series
        for i in 0..self.data.len() {
            self.data[i] += other;
        }
        // return modified left timeseries
        self
    }

}
impl<'a> Add<&'a mut FrequencySeries> for Complex<f64> {

    type Output = &'a mut FrequencySeries;

    fn add(self, other: &'a mut FrequencySeries) -> &'a mut FrequencySeries {
        // uses previous operation
        other.add(self)
    }

}

/* operator- ----------------------------------------------------------------------------------- */
impl<'a> Neg for &'a mut FrequencySeries {

	type Output = &'a mut FrequencySeries;

	fn neg(self) -> &'a mut FrequencySeries {
		// return modified
		self.mul(Complex{re: -1., im: 0.})
	}
}
impl<'a> Sub<&FrequencySeries> for &'a mut FrequencySeries {

	type Output = &'a mut FrequencySeries;

	fn sub(self, other: &FrequencySeries) -> &'a mut FrequencySeries {
        // verify if the length of the time series are equal
        assert_eq!(self.data.len(), other.data.len());
        // Modify data vector of the left time series
        for i in 0..self.data.len() {
            self.data[i] -= other.data[i];
        }
        // return modified left timeseries
        self
	}
}
impl<'a> Sub<Complex<f64>> for &'a mut FrequencySeries {

	type Output = &'a mut FrequencySeries;

	fn sub(self, other: Complex<f64>) -> &'a mut FrequencySeries {
        // return modified left timeseries
        self.add(-other)
	}
}
impl<'a> Sub<&'a mut FrequencySeries> for Complex<f64> {

	type Output = &'a mut FrequencySeries;

	fn sub(self, other: &'a mut FrequencySeries) -> &'a mut FrequencySeries {
        // return modified left timeseries
        other.neg().add(self)
	}
}

/* operator* ----------------------------------------------------------------------------------- */
impl<'a> Mul<&FrequencySeries> for &'a mut FrequencySeries {

    type Output = &'a mut FrequencySeries;

    fn mul(self, other: &FrequencySeries) -> &'a mut FrequencySeries {
        // verify if the length of the time series are equal
        assert_eq!(self.data.len(), other.data.len());
        // Modify data vector of the left time series
        for i in 0..self.data.len() {
            self.data[i] *= other.data[i];
        }
        // return modified left timeseries
        self
    }
}
impl<'a> Mul<Complex<f64>> for &'a mut FrequencySeries {

    type Output = &'a mut FrequencySeries;

    fn mul(self, other: Complex<f64>) -> &'a mut FrequencySeries {
        // Modify data vector of the left time series
        for i in 0..self.data.len() {
            self.data[i] *= other;
        }
        // return modified left timeseries
        self
    }
}
impl<'a> Mul<&'a mut FrequencySeries> for Complex<f64> {

    type Output = &'a mut FrequencySeries;

    fn mul(self, other: &'a mut FrequencySeries) -> &'a mut FrequencySeries {
        // uses previous operation
        other.mul(self)
    }

}
impl<'a> Mul<f64> for &'a mut FrequencySeries {

    type Output = &'a mut FrequencySeries;

    fn mul(self, other: f64) -> &'a mut FrequencySeries {
        // use complex multimplication
		self.mul(Complex{re: other, im: 0.})
    }
}
impl<'a> Mul<&'a mut FrequencySeries> for f64 {

    type Output = &'a mut FrequencySeries;

    fn mul(self, other: &'a mut FrequencySeries) -> &'a mut FrequencySeries {
        // uses previous operation
        other.mul(self)
    }
}

/* operator/ ----------------------------------------------------------------------------------- */
impl<'a> Div<&FrequencySeries> for &'a mut FrequencySeries {

	type Output = &'a mut FrequencySeries;

	fn div(self, other: &FrequencySeries) -> &'a mut FrequencySeries {
        // verify if the length of the time series are equal
        assert_eq!(self.data.len(), other.data.len());
        // Modify data vector of the left time series
        for i in 0..self.data.len() {
            self.data[i] /= other.data[i];
        }
        // return modified left timeseries
        self
	}
}
impl<'a> Div<Complex<f64>> for &'a mut FrequencySeries {

	type Output = &'a mut FrequencySeries;

	fn div(self, other: Complex<f64>) -> &'a mut FrequencySeries {
        // return modified left timeseries
        self.mul(1. / other)
	}
}
impl<'a> Div<&'a mut FrequencySeries> for Complex<f64> {

	type Output = &'a mut FrequencySeries;

	fn div(self, other: &'a mut FrequencySeries) -> &'a mut FrequencySeries {
        // return modified left timeseries
        other.inv().mul(self)
	}
}
impl<'a> Div<f64> for &'a mut FrequencySeries {

	type Output = &'a mut FrequencySeries;

	fn div(self, other: f64) -> &'a mut FrequencySeries {
        // return modified left timeseries
        self.mul(1. / other)
	}
}
impl<'a> Div<&'a mut FrequencySeries> for f64 {

	type Output = &'a mut FrequencySeries;

	fn div(self, other: &'a mut FrequencySeries) -> &'a mut FrequencySeries {
        // return modified left timeseries
        other.inv().mul(self)
	}
}

/* operator[] ---------------------------------------------------------------------------------- */
impl Index<usize> for FrequencySeries {
	
	type Output = Complex<f64>;
	fn index(&self, i: usize) -> &Complex<f64> {
		assert_gt!(self.data.len(), i);
		&self.data[i]
	}
}
impl IndexMut<usize> for FrequencySeries {
		
	fn index_mut(&mut self, i: usize) -> &mut Complex<f64> {
		assert_gt!(self.data.len(), i);
		&mut self.data[i]
	}
}

*/
