/* --------------------------------------------------------------------------------------------- *
 * Time series and frequency series definition
 * --------------------------------------------------------------------------------------------- */
pub mod data;

use rand::thread_rng;
use rand_distr::{Normal, StandardNormal, Distribution};
use std::f64::consts::PI;
use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, MulAssign, Div, DivAssign, Index, IndexMut};
use more_asserts::assert_gt;
use num::{Float, Complex, complex::ComplexFloat, NumCast};
use data::Data;


/* --------------------------------------------------------------------------------------------- *
 * Define structures
 * --------------------------------------------------------------------------------------------- */

/// Time series object:
/// Consists in a vector of data indexed by time.
#[derive(Debug, Clone)]
pub struct TimeSeries<D> {
    fs: f64,
    t0: f64,
    data: Vec<D>,
}



impl<D: ComplexFloat + Float + Data> TimeSeries<D> {
    /// Real signal generators:
    /// 
    /// Generates a white noise signal with a given size, sampling frequency and the noise amplitude
    /// 
    /// # Examples
    ///
    /// ```
    /// use gw_signal::timeseries::*;
    /// 
    /// // creates a white noise signal with 20000 points, sampled at 1 kHz,
    /// // the variance of the noise as 0.1 V
    /// let mut signal: TimeSeries = TimeSeries::white_noise(20000, 1e3, 0f64, 1f64);
    ///
    /// ```
    pub fn white_noise(size: usize, fs: f64, mu: D, sigma: D) -> Self 
		where StandardNormal: Distribution<D>
		{
		
		let rng = thread_rng();
        let normal = Normal::new(mu, sigma).unwrap();

        // fill data vector
        let data_vec: Vec<D> = normal.sample_iter(rng).take(size).collect();
		
        // initialize TimeSeries
        TimeSeries::from_vector(fs, 0., data_vec)
    }
}
impl<D: ComplexFloat + Data> TimeSeries<D> {
    /// Generates a sinusoidal signal.
    /// 
    /// # Examples
    ///
    /// ```
    /// use gw_signal::timeseries::*;
    ///
    /// // creates a sinusoidal signal with 20000 points, sampled at 1 kHz, using 8 bytes for each sample
    /// // The frequency, amplitudes and phase at the origin are respectively:
    /// //  5 Hz, 10 V and 0 rad.
    /// let mut signal: TimeSeries = TimeSeries::wave(20000, 1e3, 5f64, 10f64, 0f64);
    ///
    /// ```
    pub fn wave(size: usize, fs: f64, freq: D, ampl: D, phase: D) -> Self {

        let mut phi: D;
        let mut data_vec: Vec<D> = Vec::new();

        // fill data vector
        for i in 0..size {
            phi = freq * NumCast::from(2. * PI * i as f64 / fs).unwrap() + phase;
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
    /// use gw_signal::timeseries::*;
    ///
    /// // creates a constant signal with 20000 points, sampled at 1 kHz
    /// let mut signal: TimeSeries = TimeSeries::constant(20000, 1e3, 1f64);
    ///
    /// ```
    pub fn constant(size: usize, fs: f64, value: D) -> Self {

        let data_vec: Vec<D> = vec![value; size];

        // initialize TimeSeries
        TimeSeries::from_vector(fs, 0., data_vec)
    }

	/// Base onstructor, called by the other functions
    pub fn from_vector(fs: f64,
                       t0: f64,
                       input_data: Vec<D>) -> Self {

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
/// Getter functions
impl<D: ComplexFloat<Real = F>, F> TimeSeries<D> {
	
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
	pub fn get_data(&self) -> Vec<D> {
		self.data.clone()
	}
}

/// Math functions
impl<D: ComplexFloat<Real = F>, F> TimeSeries<D> {
	/// compute the inverse value of the data
	pub fn inv(&mut self) -> &mut TimeSeries<D> {
		
		for i in 0..self.data.len() {
			self.data[i] = self.data[i].recip();
		}
		self
	}
	/// compute square root of the data
	pub fn sqrt(&mut self) -> &mut TimeSeries<D> {
		
		for i in 0..self.data.len() {
			self.data[i] = self.data[i].sqrt();
		}
		self
	}
	pub fn real(&self) -> TimeSeries<F> {
		
		let mut output: Vec<F> = Vec::new();
		for i in 0..self.data.len() {
			output.push(self.data[i].re());
		}
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: output
		}
	}
	pub fn imag(&self) -> TimeSeries<F> {
		
		let mut output: Vec<F> = Vec::new();
		for i in 0..self.data.len() {
			output.push(self.data[i].im());
		}
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: output
		}
	}
	pub fn abs(&self) -> TimeSeries<F> {
		
		let mut output: Vec<F> = Vec::new();
		for i in 0..self.data.len() {
			output.push(self.data[i].abs());
		}
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: output
		}
	}
	pub fn arg(&self) -> TimeSeries<F> {
		
		let mut output: Vec<F> = Vec::new();
		for i in 0..self.data.len() {
			output.push(self.data[i].arg());
		}
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: output
		}
	}
}

/// Math functions for real time series
impl<D: Float> TimeSeries<D> {
	/// Get the maximum value of the time series
	pub fn max(&self) -> D {
		let mut output: D = D::neg_infinity();

		for sample in self.data.iter() {
			if *sample > output {
				output = *sample;
			}
		}
		output
	}
	/// Get the maximum value of the time series
	pub fn min(&self) -> D {
		let mut output: D = D::infinity();

		for sample in self.data.iter() {
			if *sample < output {
				output = *sample;
			}
		}
		output
	}
}



/* --------------------------------------------------------------------------------------------- */
/// Transform the data type into another
pub trait ToType {
	/// Transform time series into a time series of f32. Copy the data vector
	/// If the initial time series is complex, takes the real part of the data
	///
	/// # Example
	/// 
	/// ```
    /// use gw_signal::timeseries::*;
	/// 	
    /// 
    /// // creates a white noise signal
    /// let signal_f64: TimeSeries<f64> = TimeSeries::white_noise(20000, 1e3, 0f64, 1f64);
	/// let signal_f32: TimeSeries<f32> = signal_f64.to_f32();
	/// ```
	fn to_f32(&self) -> TimeSeries<f32>;
	/// Transform time series into a time series of f64. Copy the data vector
	/// If the initial time series is complex, takes the real part of the data
	///
	/// # Example
	/// 
	/// ```
    /// use gw_signal::timeseries::*;
	/// 	
    /// 
    /// // creates a white noise signal
    /// let signal_f32: TimeSeries<f32> = TimeSeries::white_noise(20000, 1e3, 0f32, 1f32);
	/// let signal_f64: TimeSeries<f64> = signal_f32.to_f64();
	/// ```
	fn to_f64(&self) -> TimeSeries<f64>;
	/// Transform time series into a time series of Complex<f32>. Copy the data vector
	///
	/// # Example
	/// 
	/// ```
    /// use gw_signal::timeseries::*;
	/// 	
    /// 
    /// // creates a white noise signal
    /// let signal_f64: TimeSeries<f64> = TimeSeries::white_noise(20000, 1e3, 0f64, 1f64);
	/// let signal_c32: TimeSeries<Complex<f32>> = signal_f64.to_c32();
	/// ```
	fn to_c32(&self) -> TimeSeries<Complex<f32>>;
	/// Transform time series into a time series of Complex<f32>. Copy the data vector
	///
	/// # Example
	/// 
	/// ```
    /// use gw_signal::timeseries::*;
	/// 	
    /// 
    /// // creates a white noise signal
    /// let signal_f64: TimeSeries<f64> = TimeSeries::white_noise(20000, 1e3, 0f64, 1f64);
	/// let signal_c64: TimeSeries<Complex<f64>> = signal_f64.to_c64();
	/// ```
	fn to_c64(&self) -> TimeSeries<Complex<f64>>;
}
impl ToType for TimeSeries<f32> {
	fn to_f32(&self) -> TimeSeries<f32> {
		self.clone()
	}
	fn to_f64(&self) -> TimeSeries<f64> {
		// initialize data vector
		let mut data: Vec<f64> = Vec::new();
		for i in 0..self.get_size() {
			data.push(self[i] as f64);
		}
		// return option with vector in it
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: data
		}
	}
	fn to_c32(&self) -> TimeSeries<Complex<f32>> {
		// initialize data vector
		let mut data: Vec<Complex<f32>> = Vec::new();
		for i in 0..self.get_size() {
			data.push(Complex{re:self[i], im: 0f32});
		}
		// return option with vector in it
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: data
		}
	}
	fn to_c64(&self) -> TimeSeries<Complex<f64>> {
		// initialize data vector
		let mut data: Vec<Complex<f64>> = Vec::new();
		for i in 0..self.get_size() {
			data.push(Complex{re:self[i] as f64, im: 0f64});
		}
		// return option with vector in it
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: data
		}
	}
}
impl ToType for TimeSeries<f64> {
	fn to_f32(&self) -> TimeSeries<f32> {
		// initialize data vector
		let mut data: Vec<f32> = Vec::new();
		for i in 0..self.get_size() {
			if self[i].abs() < f32::MAX as f64 {
				data.push(self[i] as f32);
			} else {
				data.push(f32::MAX * self[i].signum() as f32);
			}
		}
		// return option with vector in it
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: data
		}
	}
	fn to_f64(&self) -> TimeSeries<f64> {
		self.clone()
	}
	fn to_c32(&self) -> TimeSeries<Complex<f32>> {
		// initialize data vector
		let mut data: Vec<Complex<f32>> = Vec::new();
		for i in 0..self.get_size() {
			if self[i].abs() < f32::MAX as f64 {
				data.push(Complex{re: self[i] as f32, im: 0f32});
			} else {
				data.push(Complex{re: f32::MAX * self[i].signum() as f32, im: 0f32});
			}

		}
		// return option with vector in it
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: data
		}
	}
	fn to_c64(&self) -> TimeSeries<Complex<f64>> {
		// initialize data vector
		let mut data: Vec<Complex<f64>> = Vec::new();
		for i in 0..self.get_size() {
			data.push(Complex{re:self[i], im: 0f64});
		}
		// return option with vector in it
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: data
		}
	}
}
impl ToType for TimeSeries<Complex<f32>> {
	fn to_f32(&self) -> TimeSeries<f32> {
		// initialize data vector
		let mut data: Vec<f32> = Vec::new();
		for i in 0..self.get_size() {
			data.push(self[i].re);
		}
		// return option with vector in it
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: data
		}
	}
	fn to_f64(&self) -> TimeSeries<f64> {
		// initialize data vector
		let mut data: Vec<f64> = Vec::new();
		for i in 0..self.get_size() {
			data.push(self[i].re as f64);
		}
		// return option with vector in it
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: data
		}
	}
	fn to_c32(&self) -> TimeSeries<Complex<f32>> {
		self.clone()
	}
	fn to_c64(&self) -> TimeSeries<Complex<f64>> {
		// initialize data vector
		let mut data: Vec<Complex<f64>> = Vec::new();
		for i in 0..self.get_size() {
			data.push(Complex{re: self[i].re as f64, im: self[i].im as f64});
		}
		// return option with vector in it
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: data
		}
	}
}
impl ToType for TimeSeries<Complex<f64>> {
	fn to_f32(&self) -> TimeSeries<f32> {
		// initialize data vector
		let mut data: Vec<f32> = Vec::new();
		for i in 0..self.get_size() {
			if self[i].abs() < f32::MAX as f64 {
				data.push(self[i].re as f32);
			} else {
				data.push(f32::MAX * self[i].re.signum() as f32);
			}
		}
		// return option with vector in it
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: data
		}
	}
	fn to_f64(&self) -> TimeSeries<f64> {
		// initialize data vector
		let mut data: Vec<f64> = Vec::new();
		for i in 0..self.get_size() {
			data.push(self[i].re);
		}
		// return option with vector in it
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: data
		}

	}
	fn to_c32(&self) -> TimeSeries<Complex<f32>> {
		// initialize data vector
		let mut data: Vec<Complex<f32>> = Vec::new();
		for i in 0..self.get_size() {
			let (real, imag): (f32, f32);
			if self[i].re.abs() < f32::MAX as f64 {
				real = self[i].re as f32;
			} else {
				real = f32::MAX * self[i].re.signum() as f32;
			}
			if self[i].im.abs() < f32::MAX as f64 {
				imag = self[i].im as f32;
			} else {
				imag = f32::MAX * self[i].im.signum() as f32;
			}
			data.push(Complex{re: real, im: imag});
		}
		// return option with vector in it
		TimeSeries{
			fs: self.fs,
			t0: self.t0,
			data: data
		}
	}
	fn to_c64(&self) -> TimeSeries<Complex<f64>> {
		self.clone()
	}
}



/* --------------------------------------------------------------------------------------------- */
/// Operator overloading:
/// 
/// notes:
/// 
/// 	The operators DOES NOT create a new time series, they modify one of the time series parameter.


/* operator+ ----------------------------------------------------------------------------------- */
impl<'a, D> Add<&TimeSeries<D>> for &'a mut TimeSeries<D> 
	where D: AddAssign + ComplexFloat,
{
    type Output = &'a mut TimeSeries<D>;

    fn add(self, other: &TimeSeries<D>) -> &'a mut TimeSeries<D> {
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
impl<'a, D> Add<D> for &'a mut TimeSeries<D>
	where D: AddAssign + ComplexFloat,
{
    type Output = &'a mut TimeSeries<D>;

    fn add(self, other: D) -> &'a mut TimeSeries<D> {
        // Modify data vector of the left time series
        for i in 0..self.data.len() {
            self.data[i] += other;
        }
        // return modified left timeseries
        self
    }

}
impl<'a> Add<&'a mut TimeSeries<f32>> for f32 {
    type Output = &'a mut TimeSeries<f32>;
    fn add(self, other: &'a mut TimeSeries<f32>) -> &'a mut TimeSeries<f32> {
        other.add(self)
    }
}
impl<'a> Add<&'a mut TimeSeries<f64>> for f64 {
    type Output = &'a mut TimeSeries<f64>;
    fn add(self, other: &'a mut TimeSeries<f64>) -> &'a mut TimeSeries<f64> {
        other.add(self)
    }
}
impl<'a> Add<&'a mut TimeSeries<Complex<f32>>> for Complex<f32> {
    type Output = &'a mut TimeSeries<Complex<f32>>;
    fn add(self, other: &'a mut TimeSeries<Complex<f32>>) -> &'a mut TimeSeries<Complex<f32>> {
        other.add(self)
    }
}
impl<'a> Add<&'a mut TimeSeries<Complex<f64>>> for Complex<f64> {
    type Output = &'a mut TimeSeries<Complex<f64>>;
    fn add(self, other: &'a mut TimeSeries<Complex<f64>>) -> &'a mut TimeSeries<Complex<f64>> {
        other.add(self)
    }
}
/* operator- ----------------------------------------------------------------------------------- */
impl<'a, D> Neg for &'a mut TimeSeries<D>
	where D: MulAssign + ComplexFloat,
{
	type Output = &'a mut TimeSeries<D>;

    fn neg(self) -> &'a mut TimeSeries<D> {
        // Modify data vector of the left time series
        for i in 0..self.data.len() {
            self.data[i] *= NumCast::from(-1).unwrap();
        }
        // return modified left timeseries
        self
    }

}
impl<'a, D> Sub<&TimeSeries<D>> for &'a mut TimeSeries<D>
	where D: SubAssign + ComplexFloat,
{
	type Output = &'a mut TimeSeries<D>;

	fn sub(self, other: &TimeSeries<D>) -> &'a mut TimeSeries<D> {
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
impl<'a, D> Sub<D> for &'a mut TimeSeries<D>
	where D: SubAssign + ComplexFloat,
{
	type Output = &'a mut TimeSeries<D>;

    fn sub(self, other: D) -> &'a mut TimeSeries<D> {
        // Modify data vector of the left time series
        for i in 0..self.data.len() {
            self.data[i] -= other;
        }
        // return modified left timeseries
        self
    }

}
impl<'a> Sub<&'a mut TimeSeries<f32>> for f32 {
	type Output = &'a mut TimeSeries<f32>;
	fn sub(self, other: &'a mut TimeSeries<f32>) -> &'a mut TimeSeries<f32> {
        // return modified left timeseries
        other.neg().add(self)
	}
}
impl<'a> Sub<&'a mut TimeSeries<f64>> for f64 {
	type Output = &'a mut TimeSeries<f64>;
	fn sub(self, other: &'a mut TimeSeries<f64>) -> &'a mut TimeSeries<f64> {
        // return modified left timeseries
        other.neg().add(self)
	}
}
impl<'a> Sub<&'a mut TimeSeries<Complex<f32>>> for Complex<f32> {
	type Output = &'a mut TimeSeries<Complex<f32>>;
	fn sub(self, other: &'a mut TimeSeries<Complex<f32>>) -> &'a mut TimeSeries<Complex<f32>> {
        // return modified left timeseries
        other.neg().add(self)
	}
}
impl<'a> Sub<&'a mut TimeSeries<Complex<f64>>> for Complex<f64> {
	type Output = &'a mut TimeSeries<Complex<f64>>;
	fn sub(self, other: &'a mut TimeSeries<Complex<f64>>) -> &'a mut TimeSeries<Complex<f64>> {
        // return modified left timeseries
        other.neg().add(self)
	}
}

/* operator* ----------------------------------------------------------------------------------- */
impl<'a, D> Mul<&TimeSeries<D>> for &'a mut TimeSeries<D>
	where D: MulAssign + ComplexFloat,
{
    type Output = &'a mut TimeSeries<D>;

    fn mul(self, other: &TimeSeries<D>) -> &'a mut TimeSeries<D> {
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
impl<'a, D> Mul<D> for &'a mut TimeSeries<D>
	where D: MulAssign + ComplexFloat,
{
    type Output = &'a mut TimeSeries<D>;

    fn mul(self, other: D) -> &'a mut TimeSeries<D> {
        // Modify data vector of the left time series
        for i in 0..self.data.len() {
            self.data[i] *= other;
        }
        // return modified left timeseries
        self
    }
}
impl<'a> Mul<&'a mut TimeSeries<f32>> for f32 {
    type Output = &'a mut TimeSeries<f32>;
    fn mul(self, other: &'a mut TimeSeries<f32>) -> &'a mut TimeSeries<f32> {
        other.mul(self)
    }
}
impl<'a> Mul<&'a mut TimeSeries<f64>> for f64 {
    type Output = &'a mut TimeSeries<f64>;
    fn mul(self, other: &'a mut TimeSeries<f64>) -> &'a mut TimeSeries<f64> {
        other.mul(self)
    }
}
impl<'a> Mul<&'a mut TimeSeries<Complex<f32>>> for Complex<f32> {
    type Output = &'a mut TimeSeries<Complex<f32>>;
    fn mul(self, other: &'a mut TimeSeries<Complex<f32>>) -> &'a mut TimeSeries<Complex<f32>> {
        other.mul(self)
    }
}
impl<'a> Mul<&'a mut TimeSeries<Complex<f64>>> for Complex<f64> {
    type Output = &'a mut TimeSeries<Complex<f64>>;
    fn mul(self, other: &'a mut TimeSeries<Complex<f64>>) -> &'a mut TimeSeries<Complex<f64>> {
        other.mul(self)
    }
}


/* operator/ ----------------------------------------------------------------------------------- */
impl<'a, D> Div<&TimeSeries<D>> for &'a mut TimeSeries<D> 
	where D: DivAssign + ComplexFloat,
{
	type Output = &'a mut TimeSeries<D>;

	fn div(self, other: &TimeSeries<D>) -> &'a mut TimeSeries<D> {

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
impl<'a, D> Div<D> for &'a mut TimeSeries<D> 
	where D: DivAssign + ComplexFloat,
{
	type Output = &'a mut TimeSeries<D>;

	fn div(self, other: D) -> &'a mut TimeSeries<D> {
		// Modify data vector
		for i in 0..self.data.len() {
			self.data[i] /= other;
		}
		// return modified left timeseries
		self
	}
}
impl<'a> Div<&'a mut TimeSeries<f32>> for f32 {
	type Output = &'a mut TimeSeries<f32>;
	fn div(self, other: &'a mut TimeSeries<f32>) -> &'a mut TimeSeries<f32> {
        // return modified left timeseries
        other.inv().mul(self)
	}
}
impl<'a> Div<&'a mut TimeSeries<f64>> for f64 {
	type Output = &'a mut TimeSeries<f64>;
	fn div(self, other: &'a mut TimeSeries<f64>) -> &'a mut TimeSeries<f64> {
        // return modified left timeseries
        other.inv().mul(self)
	}
}
impl<'a> Div<&'a mut TimeSeries<Complex<f32>>> for Complex<f32> {
	type Output = &'a mut TimeSeries<Complex<f32>>;
	fn div(self, other: &'a mut TimeSeries<Complex<f32>>) -> &'a mut TimeSeries<Complex<f32>> {
        // return modified left timeseries
        other.inv().mul(self)
	}
}
impl<'a> Div<&'a mut TimeSeries<Complex<f64>>> for Complex<f64> {
	type Output = &'a mut TimeSeries<Complex<f64>>;
	fn div(self, other: &'a mut TimeSeries<Complex<f64>>) -> &'a mut TimeSeries<Complex<f64>> {
        // return modified left timeseries
        other.inv().mul(self)
	}
}


/* operator[] ---------------------------------------------------------------------------------- */
impl<D: ComplexFloat + ComplexFloat> Index<usize> for TimeSeries<D> {
	
	type Output = D;
	fn index(&self, i: usize) -> &D {
		assert_gt!(self.data.len(), i);
		&self.data[i]
	}
}

impl<D: ComplexFloat + ComplexFloat> IndexMut<usize> for TimeSeries<D> {
		
	fn index_mut(&mut self, i: usize) -> &mut D {
		assert_gt!(self.data.len(), i);
		&mut self.data[i]
	}
}




