/* --------------------------------------------------------------------------------------------- *
 * Time series and frequency series definition
 * --------------------------------------------------------------------------------------------- */

use rand::thread_rng;
use rand_distr::{Normal, Distribution};
use std::f64::consts::PI;
use std::ops::{Add, Neg, Sub, Mul, Div, Index, IndexMut};
use more_asserts::assert_gt;

/* --------------------------------------------------------------------------------------------- *
 * Define structures
 * --------------------------------------------------------------------------------------------- */

/// Time series object:
/// Consists in a vector of data indexed by time.
#[derive(Debug, Clone)]
pub struct TimeSeries {
    fs: f64,
    t0: f64,
    data: Vec<f64>,
}



impl TimeSeries {
    /// Signal generators:
    /// 
    /// Generates a white noise signal with a given size, sampling frequency and the noise amplitude
    /// 
    /// # Examples
    ///
    /// ```
    /// use ::timeseries::timeseries as ts;
    /// 
    /// // creates a white noise signal with 20000 points, sampled at 1 kHz,
    /// // the variance of the noise as 0.1 V
    /// let mut signal: ts::TimeSeries = ts::TimeSeries::white_noise(20000, 1e3, 1e-1);
    ///
    /// ```
    pub fn white_noise(size: usize, fs: f64, sigma: f64) -> Self {

        let normal = Normal::new(0., sigma).unwrap();
        let mut data_vec: Vec<f64> = Vec::new();

        // fill data vector
        for _i in 0..size {
            data_vec.push(normal.sample( &mut thread_rng() ));
        }
        // initialize TimeSeries
        TimeSeries::from_vector(fs, 0., data_vec)
    }

    /// Generates a sinusoidal signal.
    /// 
    /// # Examples
    ///
    /// ```
    /// use ::timeseries::timeseries as ts;
    ///
    /// // creates a sinusoidal signal with 20000 points, sampled at 1 kHz
    /// // The frequency, amplitudes and phase at the origin are respectively:
    /// //  5 Hz, 10 V and 0 rad.
    /// let mut signal: ts::TimeSeries = ts::TimeSeries::wave(20000, 1e3, 5., 10., 0.);
    ///
    /// ```
    pub fn wave(size: usize, fs: f64, freq: f64, ampl: f64, phase: f64) -> Self {

        let mut phi: f64;
        let mut data_vec: Vec<f64> = Vec::new();

        // fill data vector
        for i in 0..size {
            phi = 2. * PI * freq * (i as f64 / fs) + phase;
            data_vec.push(ampl * phi.cos());
        }
        // initialize TimeSeries
        TimeSeries::from_vector(fs, 0., data_vec)
    }

    /// Generates a constant signal.
    /// 
    /// # Examples
    ///
    /// ```
    /// use ::timeseries::timeseries as ts;
    ///
    /// // creates a constant signal with 20000 points, sampled at 1 kHz
    /// let mut signal: ts::TimeSeries = ts::TimeSeries::constant(20000, 1e3, 1.);
    ///
    /// ```
    pub fn constant(size: usize, fs: f64, value: f64) -> Self {

        let data_vec: Vec<f64> = vec![value; size];

        // initialize TimeSeries
        TimeSeries::from_vector(fs, 0., data_vec)
    }
	/// Base onstructor, called by the other functions
    pub fn from_vector(fs: f64,
                       t0: f64,
                       input_data: Vec<f64>) -> Self {

        // define time series
        TimeSeries {
            fs,
            t0,
            data: input_data
        }

    }


}

/* --------------------------------------------------------------------------------------------- *
 * traits
 * --------------------------------------------------------------------------------------------- */

/* Getter trait -------------------------------------------------------------------------------- */
impl TimeSeries {
	
	pub fn get_size(&self) -> usize {
		self.data.len()
	}

	pub fn get_fs(&self) -> f64 {
		self.fs
	}

	pub fn get_t0(&self) -> f64 {
		self.t0
	}
	/// get data vector
	pub fn get_data(&self) -> Vec<f64> {
		self.data.clone()
	}
	/// compute the inverse value of the data
	pub fn inv(&mut self) -> &mut TimeSeries {
		
		for i in 0..self.data.len() {
			self.data[i] = 1. / self.data[i];
		}
		self
	}
	/// compute square root of the data
	pub fn sqrt(mut self) -> TimeSeries {
		
		for i in 0..self.data.len() {
			self.data[i] = self.data[i].sqrt();
		}
		self
	}
}


/// Operator overloading:
/// 
/// notes:
/// 
/// 	The operators DOES NOT create a new time series, they modify one of the time series parameter.


/* operator+ ----------------------------------------------------------------------------------- */
impl<'a> Add<&TimeSeries> for &'a mut TimeSeries {

    type Output = &'a mut TimeSeries;

    fn add(self, other: &TimeSeries) -> &'a mut TimeSeries {
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
impl<'a> Add<f64> for &'a mut TimeSeries {

    type Output = &'a mut TimeSeries;

    fn add(self, other: f64) -> &'a mut TimeSeries {
        // Modify data vector of the left time series
        for i in 0..self.data.len() {
            self.data[i] += other;
        }
        // return modified left timeseries
        self
    }

}
impl<'a> Add<&'a mut TimeSeries> for f64 {

    type Output = &'a mut TimeSeries;

    fn add(self, other: &'a mut TimeSeries) -> &'a mut TimeSeries {
        // uses previous operation
        other.add(self)
    }

}


/* operator- ----------------------------------------------------------------------------------- */
impl<'a> Neg for &'a mut TimeSeries {

	type Output = &'a mut TimeSeries;

	fn neg(self) -> &'a mut TimeSeries {
		// return modified
		self.mul(-1.)
	}
}
impl<'a> Sub<&TimeSeries> for &'a mut TimeSeries {

	type Output = &'a mut TimeSeries;

	fn sub(self, other: &TimeSeries) -> &'a mut TimeSeries {
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
impl<'a> Sub<f64> for &'a mut TimeSeries {

	type Output = &'a mut TimeSeries;

	fn sub(self, other: f64) -> &'a mut TimeSeries {
        // return modified left timeseries
        self.add(-other)
	}
}
impl<'a> Sub<&'a mut TimeSeries> for f64 {

	type Output = &'a mut TimeSeries;

	fn sub(self, other: &'a mut TimeSeries) -> &'a mut TimeSeries {
        // return modified left timeseries
        other.neg().add(self)
	}
}


/* operator* ----------------------------------------------------------------------------------- */
impl<'a> Mul<&TimeSeries> for &'a mut TimeSeries {

    type Output = &'a mut TimeSeries;

    fn mul(self, other: &TimeSeries) -> &'a mut TimeSeries {
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
impl<'a> Mul<f64> for &'a mut TimeSeries {

    type Output = &'a mut TimeSeries;

    fn mul(self, other: f64) -> &'a mut TimeSeries {
        // Modify data vector of the left time series
        for i in 0..self.data.len() {
            self.data[i] *= other;
        }
        // return modified left timeseries
        self
    }
}
impl<'a> Mul<&'a mut TimeSeries> for f64 {

    type Output = &'a mut TimeSeries;

    fn mul(self, other: &'a mut TimeSeries) -> &'a mut TimeSeries {
        // uses previous operation
        other.mul(self)
    }
}


/* operator/ ----------------------------------------------------------------------------------- */
impl<'a> Div<&TimeSeries> for &'a mut TimeSeries {

	type Output = &'a mut TimeSeries;

	fn div(self, other: &TimeSeries) -> &'a mut TimeSeries {
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
impl<'a> Div<f64> for &'a mut TimeSeries {

	type Output = &'a mut TimeSeries;

	fn div(self, other: f64) -> &'a mut TimeSeries {
        // return modified left timeseries
        self.mul(1. / other)
	}
}
impl<'a> Div<&'a mut TimeSeries> for f64 {

	type Output = &'a mut TimeSeries;

	fn div(self, other: &'a mut TimeSeries) -> &'a mut TimeSeries {
        // return modified left timeseries
        other.inv().mul(self)
	}
}


/* operator[] ---------------------------------------------------------------------------------- */
impl Index<usize> for TimeSeries {
	
	type Output = f64;
	fn index(&self, i: usize) -> &f64 {
		assert_gt!(self.data.len(), i);
		&self.data[i]
	}
}

impl IndexMut<usize> for TimeSeries {
		
	fn index_mut(&mut self, i: usize) -> &mut f64 {
		assert_gt!(self.data.len(), i);
		&mut self.data[i]
	}
}




