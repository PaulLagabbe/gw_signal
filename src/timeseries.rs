/* --------------------------------------------------------------------------------------------- *
 * Time series and frequency series definition
 * --------------------------------------------------------------------------------------------- */

use rand::thread_rng;
use rand_distr::{Normal, StandardNormal, Distribution};
use std::f64::consts::PI;
use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, MulAssign, Div, DivAssign, Index, IndexMut};
use more_asserts::assert_gt;
use num::{Float, Complex, complex::ComplexFloat, NumCast};



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



impl<D: Float> TimeSeries<D> {

    /// Real signal generators:
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
    /// use ::timeseries::timeseries as ts;
    ///
    /// // creates a constant signal with 20000 points, sampled at 1 kHz
    /// let mut signal: ts::TimeSeries = ts::TimeSeries::constant(20000, 1e3, 1.);
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
impl<D: ComplexFloat> TimeSeries<D> {
	
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
	/// compute the inverse value of the data
	pub fn inv(&mut self) -> &mut TimeSeries<D> {
		
		for i in 0..self.data.len() {
			self.data[i] = self.data[i].recip();
		}
		self
	}
	/// compute square root of the data
	pub fn sqrt(mut self) -> TimeSeries<D> {
		
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




