/* --------------------------------------------------------------------------------------------- *
 * Frequency series definition
 * --------------------------------------------------------------------------------------------- */

use num::complex::{
	Complex,
	ComplexFloat
};
use std::ops::{Add, Neg, Sub, Mul, Div};

/* --------------------------------------------------------------------------------------------- *
 * Define structures
 * --------------------------------------------------------------------------------------------- */
/// Frequency series object
/// The data vector is a complex value vector indexed by the frequency
/// 
/// The frequency vector start from 0 Hz.
#[derive(Debug, Clone)]
pub struct Spectrogram {
    f_max: f64,
	t0: f64,
	dt: f64,
    data: Vec<Vec<Complex<f64>>>,
}



/* --------------------------------------------------------------------------------------------- *
 * Constructors
 * --------------------------------------------------------------------------------------------- */

/// Getter functions
impl Spectrogram {
	/// Constructor
	pub fn from_vector(f_max: f64, t0: f64, dt: f64, input_data: Vec<Vec<Complex<f64>>>) -> Self {

        Spectrogram {
            f_max,
			t0,
			dt,
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
	/// get starting time
	pub fn get_t0(&self) -> f64 {
		self.t0
	}

	/// get time step
	pub fn get_dt(&self) -> f64 {
		self.dt
	}
	/// get data vector
	pub fn get_data(&self) -> Vec<Vec<Complex<f64>>> {
		self.data.clone()
	}
}

/// Mathematical function
impl Spectrogram {
	/// compute real part of the data vector
	pub fn real(mut self) -> Spectrogram {
		// iterate over data
        for it_vec in self.data.iter_mut() {
        	// Modify data vector of the left time series
			for it_val in it_vec.iter_mut() {
				*it_val = Complex{re: it_val.re, im: 0f64};
			}
        }
        // return modified left timeseries
		self
	}
	/// compute imaginary part of the data vector
	pub fn imag(mut self) -> Spectrogram {
		// iterate over data
        for it_vec in self.data.iter_mut() {
        	// Modify data vector of the left time series
			for it_val in it_vec.iter_mut() {
				*it_val = Complex{re: it_val.im, im: 0f64};
			}
        }
        // return modified left timeseries
		self
	}
	/// compute absolute value of the data vector
	pub fn abs(mut self) -> Spectrogram {
		// iterate over data
        for it_vec in self.data.iter_mut() {
        	// Modify data vector of the left time series
			for it_val in it_vec.iter_mut() {
				*it_val = Complex{re: it_val.norm(), im: 0f64};
			}
        }
        // return modified left timeseries
		self
	}
	/// compute square of absolute value of the data vector
	pub fn abs2(mut self) -> Spectrogram {
		// iterate over data
        for it_vec in self.data.iter_mut() {
        	// Modify data vector of the left time series
			for it_val in it_vec.iter_mut() {
				*it_val = *it_val * it_val.conj();
			}
        }
        // return modified left timeseries
		self
	}
	/// compute phase of the data vector
	pub fn arg(mut self) -> Spectrogram {
		// iterate over data
        for it_vec in self.data.iter_mut() {
        	// Modify data vector of the left time series
			for it_val in it_vec.iter_mut() {
				*it_val = Complex{re: it_val.arg(), im: 0f64};
			}
        }
        // return modified left timeseries
		self
	}
	/// compute complex conjugate of the data vector
	pub fn conj(mut self) -> Spectrogram {
		// iterate over data
        for it_vec in self.data.iter_mut() {
        	// Modify data vector of the left time series
			for it_val in it_vec.iter_mut() {
				*it_val = it_val.conj();
			}
        }
        // return modified left timeseries
		self
	}
	/// compute square root of the data vector
	pub fn sqrt(mut self) -> Spectrogram {
		// iterate over data
        for it_vec in self.data.iter_mut() {
        	// Modify data vector of the left time series
			for it_val in it_vec.iter_mut() {
				*it_val = it_val.sqrt();
			}
        }
        // return modified left timeseries
		self
	}
	/// compute inverse of the data vector
	pub fn inv(&mut self) -> &mut Spectrogram {
		// iterate over data
        for it_vec in self.data.iter_mut() {
        	// Modify data vector of the left time series
			for it_val in it_vec.iter_mut() {
				*it_val = it_val.inv();
			}
        }
        // return modified left timeseries
		self
	}
}


/* --------------------------------------------------------------------------------------------- *
 * Operator overloading:
 * notes:
 *      The operators DOES NOT create a new time series, they modify one of the time series
 *      parameter.
 * --------------------------------------------------------------------------------------------- */


/* operator+ ----------------------------------------------------------------------------------- */
impl<'a> Add<&Spectrogram> for &'a mut Spectrogram {

    type Output = &'a mut Spectrogram;

    fn add(self, other: &Spectrogram) -> &'a mut Spectrogram {
		// iterate over data
        for it_vec in other.data.iter().zip(self.data.iter_mut()) {
			let (vec1, vec2) = it_vec;
        	// verify if the length of the time series are equal
        	assert_eq!(vec1.len(), vec2.len());
        	// Modify data vector of the left time series
			for it_val in vec1.iter().zip(vec2.iter_mut()) {
            	let (val1, val2) = it_val;
				*val2 += *val1;
			}
        }
        // return modified left timeseries
        self
    }

}
impl<'a> Add<Complex<f64>> for &'a mut Spectrogram {

    type Output = &'a mut Spectrogram;

    fn add(self, other: Complex<f64>) -> &'a mut Spectrogram {
		// iterate over data
        for it_vec in self.data.iter_mut() {
        	// Modify data vector of the left time series
			for it_val in it_vec.iter_mut() {
				*it_val += other;
			}
        }
        // return modified left timeseries
        self
    }

}
impl<'a> Add<&'a mut Spectrogram> for Complex<f64> {

    type Output = &'a mut Spectrogram;

    fn add(self, other: &'a mut Spectrogram) -> &'a mut Spectrogram {
        // uses previous operation
        other.add(self)
    }

}

/* operator- ----------------------------------------------------------------------------------- */
impl<'a> Neg for &'a mut Spectrogram {

	type Output = &'a mut Spectrogram;

	fn neg(self) -> &'a mut Spectrogram {
		// return modified
		self.mul(Complex{re: -1., im: 0.})
	}
}
impl<'a> Sub<&Spectrogram> for &'a mut Spectrogram {

	type Output = &'a mut Spectrogram;

	fn sub(self, other: &Spectrogram) -> &'a mut Spectrogram {
  		// iterate over data
        for it_vec in other.data.iter().zip(self.data.iter_mut()) {
			let (vec1, vec2) = it_vec;
        	// verify if the length of the time series are equal
        	assert_eq!(vec1.len(), vec2.len());
        	// Modify data vector of the left time series
			for it_val in vec1.iter().zip(vec2.iter_mut()) {
            	let (val1, val2) = it_val;
				*val2 -= *val1;
			}
        }
        // return modified left timeseries
        self
	}
}
impl<'a> Sub<Complex<f64>> for &'a mut Spectrogram {

	type Output = &'a mut Spectrogram;

	fn sub(self, other: Complex<f64>) -> &'a mut Spectrogram {
        // return modified left timeseries
        self.add(-other)
	}
}
impl<'a> Sub<&'a mut Spectrogram> for Complex<f64> {

	type Output = &'a mut Spectrogram;

	fn sub(self, other: &'a mut Spectrogram) -> &'a mut Spectrogram {
        // return modified left timeseries
        other.neg().add(self)
	}
}

/* operator* ----------------------------------------------------------------------------------- */
impl<'a> Mul<&Spectrogram> for &'a mut Spectrogram {

    type Output = &'a mut Spectrogram;

    fn mul(self, other: &Spectrogram) -> &'a mut Spectrogram {
   		// iterate over data
        for it_vec in other.data.iter().zip(self.data.iter_mut()) {
			let (vec1, vec2) = it_vec;
        	// verify if the length of the time series are equal
        	assert_eq!(vec1.len(), vec2.len());
        	// Modify data vector of the left time series
			for it_val in vec1.iter().zip(vec2.iter_mut()) {
            	let (val1, val2) = it_val;
				*val2 *= *val1;
			}
        }
        // return modified left timeseries
        self
    }
}
impl<'a> Mul<Complex<f64>> for &'a mut Spectrogram {

    type Output = &'a mut Spectrogram;

    fn mul(self, other: Complex<f64>) -> &'a mut Spectrogram {
  		// iterate over data
        for it_vec in self.data.iter_mut() {
        	// Modify data vector of the left time series
			for it_val in it_vec.iter_mut() {
				*it_val *= other;
			}
        }
        // return modified left timeseries
        self
    }
}
impl<'a> Mul<&'a mut Spectrogram> for Complex<f64> {

    type Output = &'a mut Spectrogram;

    fn mul(self, other: &'a mut Spectrogram) -> &'a mut Spectrogram {
        // uses previous operation
        other.mul(self)
    }

}
impl<'a> Mul<f64> for &'a mut Spectrogram {

    type Output = &'a mut Spectrogram;

    fn mul(self, other: f64) -> &'a mut Spectrogram {
        // use complex multimplication
		self.mul(Complex{re: other, im: 0.})
    }
}
impl<'a> Mul<&'a mut Spectrogram> for f64 {

    type Output = &'a mut Spectrogram;

    fn mul(self, other: &'a mut Spectrogram) -> &'a mut Spectrogram {
        // uses previous operation
        other.mul(self)
    }
}

/* operator/ ----------------------------------------------------------------------------------- */
impl<'a> Div<&Spectrogram> for &'a mut Spectrogram {

	type Output = &'a mut Spectrogram;

	fn div(self, other: &Spectrogram) -> &'a mut Spectrogram {
  		// iterate over data
        for it_vec in other.data.iter().zip(self.data.iter_mut()) {
			let (vec1, vec2) = it_vec;
        	// verify if the length of the time series are equal
        	assert_eq!(vec1.len(), vec2.len());
        	// Modify data vector of the left time series
			for it_val in vec1.iter().zip(vec2.iter_mut()) {
            	let (val1, val2) = it_val;
				*val2 /= *val1;
			}
        }
        // return modified left timeseries
        self
	}
}
impl<'a> Div<Complex<f64>> for &'a mut Spectrogram {

	type Output = &'a mut Spectrogram;

	fn div(self, other: Complex<f64>) -> &'a mut Spectrogram {
        // return modified left timeseries
        self.mul(1. / other)
	}
}
impl<'a> Div<&'a mut Spectrogram> for Complex<f64> {

	type Output = &'a mut Spectrogram;

	fn div(self, other: &'a mut Spectrogram) -> &'a mut Spectrogram {
        // return modified left timeseries
        other.inv().mul(self)
	}
}
impl<'a> Div<f64> for &'a mut Spectrogram {

	type Output = &'a mut Spectrogram;

	fn div(self, other: f64) -> &'a mut Spectrogram {
        // return modified left timeseries
        self.mul(1. / other)
	}
}
impl<'a> Div<&'a mut Spectrogram> for f64 {

	type Output = &'a mut Spectrogram;

	fn div(self, other: &'a mut Spectrogram) -> &'a mut Spectrogram {
        // return modified left timeseries
        other.inv().mul(self)
	}
}

/* operator[] ---------------------------------------------------------------------------------- */
/*
impl Index<usize> for Spectrogram {
	
	type Output = Complex<f64>;
	fn index(&self, i: usize) -> &Complex<f64> {
		assert_gt!(self.data.len(), i);
		&self.data[i]
	}
}
impl IndexMut<usize> for Spectrogram {
		
	fn index_mut(&mut self, i: usize) -> &mut Complex<f64> {
		assert_gt!(self.data.len(), i);
		&mut self.data[i]
	}
}
*/
