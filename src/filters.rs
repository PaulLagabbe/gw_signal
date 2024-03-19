/* --------------------------------------------------------------------------------------------- *
 * Libraries
 * --------------------------------------------------------------------------------------------- */

use std::f64::consts::PI;
use rustfft::num_complex::Complex;
use more_asserts as ma;


/// IIR Filter object. This object modelizes both the IIR and the FIR filters.
/// A filter can be defined in two ways:
///
/// - By calling a standard filter constructors such as Butterworth or Chebyshev
/// 
/// - By initializing an allpass filter and adding poles and zeros at the wanted frequencies
///
/// A filter is defined by its Laplace transform. When the filter is applied to a signal,
/// the bilinear transform function as called in order to compute its z-transform.
#[derive(Debug, Clone)]
pub struct Filter {
	gain: f64,
	poles: Vec<Complex<f64>>,
	zeros: Vec<Complex<f64>>,
	fs: f64
}
/// Bandwidth type, to define the filter bandpass and cutoff frequency(ies)
/// low pass, high pass and band pass
pub enum BType {
	LowPass(f64),
	HighPass(f64),
	BandPass(f64, f64),
}


impl Filter {

	/// Bilinear transform function
	/// Compute the z-tranform of the filter from its Laplace transform, using the Tustin's method.
	/// This function is called when the filter is applied to a time series using the [apply_filter]
	/// method.
	/// [apply_filter]: <https://docs.rs/gw_signal/0.1.2/gw_signal/timeseries/struct.TimeSeries.html>
	pub fn bilinear_transform(&mut self) {

		// there should be more pole than zeros
		ma::assert_ge!(self.poles.len(), self.zeros.len());

		// compute the bilinear transform of the pole and zeros
		let mut gain_factor: Complex<f64> = Complex{re: 1., im: 0.};
		for i in 0..self.poles.len() {
			gain_factor /= 2. * self.fs - self.poles[i];
			self.poles[i] = (2. * self.fs + self.poles[i]) / (2. * self.fs - self.poles[i]);
		}
		for i in 0..self.zeros.len() {
			gain_factor *= 2. * self.fs - self.zeros[i];
			self.zeros[i] = (2. * self.fs + self.zeros[i]) / (2. * self.fs - self.zeros[i]);
		}

		self.gain *= gain_factor.re;
		
		// complete the list of zeros with -1., so the two lists have the same size
		self.zeros.append(&mut vec![Complex{re:-1., im: 0.}; self.poles.len()-self.zeros.len()]);
	}

	/* ----------------------------------------------------------------------------------------- */
	fn adapt_frequencies_lp(&mut self, factor: f64) {
	
		/* ------------------------------------------------------------------------------------- *
		 * Scale the pole and zeros frequencies, and the gain of the filter by a factor given in
		 * parameters
		 * ------------------------------------------------------------------------------------- */

        let dim: usize = self.poles.len() - self.zeros.len();

        // adapt the poles, zeros and gain values to the cutoff frequency
        for i in 0..self.zeros.len() {
            self.zeros[i] *= factor;
        }
        for i in 0..self.poles.len() {
            self.poles[i] *= factor;
        }
        self.gain *= factor.powi(dim as i32);
	}

	
	/// private function
	fn adapt_frequencies_hp(&mut self, factor: f64) {
	
        let dim: usize = self.poles.len() - self.zeros.len();

        // adapt the poles, zeros and gain values to the cutoff frequency
		let mut gain_factor: Complex<f64> = Complex{re: 1., im: 0.};
        for i in 0..self.zeros.len() {
			gain_factor *= -self.zeros[i];
            self.zeros[i] = factor / self.zeros[i];
        }
		self.zeros.append(&mut vec![Complex{re:0., im:0.}; dim]);
        for i in 0..self.poles.len() {
			gain_factor /= -self.poles[i];
            self.poles[i] = factor / self.poles[i];
        }
        self.gain *= gain_factor.re;
	}

	/// private function
	fn adapt_frequencies_bp(&mut self, factor: f64, width: f64) {
	
		let dim: usize = self.poles.len() - self.zeros.len();
		let mut fshift: Complex<f64>;
        
		// adapt the poles, zeros and gain values to the cutoff frequency
        for i in 0..self.zeros.len() {
            self.zeros[i] *= width / 2.;
			fshift = (self.zeros[i] * self.zeros[i] - factor * factor).sqrt();
			self.zeros.push(self.zeros[i] - fshift);
			self.zeros[i] += fshift;
        }
		self.zeros.append(&mut vec![Complex{re:0., im:0.}; dim]);
        
		for i in 0..self.poles.len() {
            self.poles[i] *= width / 2.;
			fshift = (self.poles[i] * self.poles[i] - factor * factor).sqrt();
			self.poles.push(self.poles[i] - fshift);
			self.poles[i] += fshift;
        }
        
		self.gain *= width.powi(dim as i32);
	}

	/// private function
	/// Transform the filter into a lowpass, highpass, bandpass... and the given cutoff frequency(ies)
	fn adapt_frequencies(&mut self, band: BType) {

		match band {
			BType::LowPass(cutoff) => {
				let omega: f64 = 2. * self.fs * (PI * cutoff / self.fs).tan();
				self.adapt_frequencies_lp(omega);
			}
			BType::HighPass(cutoff) => {
				let omega: f64 = 2. * self.fs * (PI * cutoff / self.fs).tan();
				self.adapt_frequencies_hp(omega);
			}
			BType::BandPass(cutoff_1, cutoff_2) => {
				ma::assert_ge!(cutoff_2, cutoff_1);
				let omega_1: f64 = 2. * self.fs * (PI * cutoff_1 / self.fs).tan();
				let omega_2: f64 = 2. * self.fs * (PI * cutoff_2 / self.fs).tan();
				self.adapt_frequencies_bp((omega_1 * omega_2).sqrt(), omega_2 - omega_1);
			}
		}

	}

	/// private function
	fn zero2coef(mut zeros: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
		
		/* Convert the zeros of a polynome into its coefficiants
		 * recursive function !! */

		if zeros.len() < 1 {
			vec![Complex{re: 1., im: 0.}]
		}
		else if zeros.len() == 1 {
			vec![Complex{re: 1., im: 0.}, -zeros[0]]
		}
		else {
			let z: Complex<f64> = zeros.remove(0);
			let temp: Vec<Complex<f64>> = Self::zero2coef(zeros);
			let mut output: Vec<Complex<f64>> = vec![temp[0]];
			
			for i in 0..(temp.len() - 1) {
				output.push(temp[i+1] - z * temp[i]);
			}
			output.push(-z * temp[temp.len()-1]);
			output
		}
	}
	
	/// Transform a filter defined by its poles and zeros, into its polynomial coefficiants
	/// This function is called when the filter is applied to a time series using the [apply_filter]
	/// method.
	pub fn polezero_to_coef(&self) -> (Vec<f64>, Vec<f64>) {
		
		let a_comp: Vec<Complex<f64>> = Self::zero2coef(self.poles.clone());
		let b_comp: Vec<Complex<f64>> = Self::zero2coef(self.zeros.clone());
		let mut a: Vec<f64> = Vec::new();
		let mut b: Vec<f64> = Vec::new();

		for i in 0..b_comp.len() {
			b.push(self.gain * b_comp[i].re);
		}
		for j in 0..a_comp.len() {
			a.push(a_comp[j].re);
		}
		
		(b, a)
	}

	/// Standard filter generators: These methods are Filter object constructors.
	/// So far, Butterworth, types 1 and 2 Chebyshev filters, low pass, high pass and band pass
	/// can be generated.
	/// 
	/// 
	/// Generates butterworth filter
	/// 
	/// # Example
    /// ```
    /// use ::timeseries::{
    ///     filter as flt,
    /// };
    /// 
	/// let fs: f64 = 1e3;
	/// 
    /// // generates an 8th butterworth lowpass filter at 10 Hz
    /// let my_filter: flt::Filter::butterworth(8, flt::BType::LowType(10.), fs);
    ///
    /// ```

	pub fn butterworth(order: usize, band: BType, fs: f64) -> Filter {
		/* ------------------------------------------------------------------------------------- *
		 * Initialize Butterworth filter with the order given in parameter
		 * ------------------------------------------------------------------------------------- */

		// initialize filter
		let mut output: Self = Filter{
			gain: 1.,
			poles: Vec::new(),
			zeros: Vec::new(),
			fs: fs
		};
		
		// compute filter poles
		for i in 0..order {
			output.poles.push(-Complex{
				re: 0.,
				im: PI * (-(order as i32) + 1 + 2 * (i as i32)) as f64 / (2 * order) as f64
			}.exp());
		}
		
		// adapt the poles, zeros and gain values to the cutoff frequency
		output.adapt_frequencies(band);
		output
	}
	/// Generates type 1 Chebyshev filter
	/// 
	/// # Example
    /// ```
    /// use ::timeseries::{
    ///     filter as flt,
    /// };
    /// 
	/// let fs: f64 = 1e3;
	/// 
    /// // generates an 8th butterworth lowpass filter at 10 Hz, with 0.1 dB of ripple in the band pass
    /// let my_filter: flt::Filter::chebyshev_type1(8, 0.1, flt::BType::LowType(10.), fs);
    ///
    /// ```

	pub fn chebyshev_type1(order: usize, ripple: f64, band: BType, fs: f64) -> Filter {
		/* ------------------------------------------------------------------------------------- *
		 * Type 1 Chebyshev filter
		 * the order and the ripple factor (in dB) are given in parameter
		 * ------------------------------------------------------------------------------------- */

		let eps: f64 = ((10_f64).powf(ripple/10.) - 1.).sqrt();
		let mu: f64 = (1. / eps).asinh() / (order as f64);

		// initialize filter
		let mut output: Self = Filter{
			gain: 1.,
			poles: Vec::new(),
			zeros: Vec::new(),
			fs: fs
		};
		
		// compute filter poles
		let mut comp_gain: Complex<f64> = Complex{re:1., im:0.};
		let mut phi: Complex<f64>;
		for i in 0..order {
			phi = Complex{
				re: mu,
				im: PI * (-(order as i32) + 1 + 2 * (i as i32)) as f64 / (2 * order) as f64
			};
			output.poles.push(-phi.sinh());
			comp_gain *= phi.sinh();
		}
		
		// compute filter gain
		output.gain = comp_gain.re;

		if order % 2 == 0 {
			output.gain /= (1. + eps * eps).sqrt();
		}
		// adapt the poles, zeros and gain values to the cutoff frequency
		output.adapt_frequencies(band);
		output
	}
	
	/// Generates type 2 Chebyshev filter
	/// 
	/// # Example
    /// ```
    /// use ::timeseries::{
    ///     filter as flt,
    /// };
    /// 
	/// let fs: f64 = 1e3;
	/// 
    /// // generates an 8th lowpass type 2 Chebyshev filter at 10 Hz, with an attenuation of -40 dB output of the pand pass.
    /// let my_filter: flt::Filter::chebyshev_type2(8, 40., flt::BType::LowType(10.), fs);
    ///
    /// ```
	pub fn chebyshev_type2(order: usize, attenuation: f64, band: BType, fs: f64) -> Filter {
		/* ------------------------------------------------------------------------------------- *
		 * Type 1 Chebyshev filter
		 * the order and the attenuation factor (in dB) are given in parameter
		 * ------------------------------------------------------------------------------------- */

		let eps: f64 = 1. / ((10_f64).powf(attenuation/10.) - 1.).sqrt();
		let mu: f64 = (1. / eps).asinh() / (order as f64);

		// initialize filter
		let mut output: Self = Filter{
			gain: 1.,
			poles: Vec::new(),
			zeros: Vec::new(),
			fs: fs
		};
		

		// compute filter zeros
		let mut zero: Complex<f64>;
		let mut comp_gain: Complex<f64> = Complex{re:1., im:0.};
		let mut index: i32;

		for i in 0..order {
			index = -(order as i32) + 1 + 2 * (i as i32);
			if index != 0 {
				zero = Complex{
					re: 0.,
					im: 1. / (PI * (index as f64) / (2 * order) as f64).sin()
				};
				output.zeros.push(zero);
				comp_gain /= -zero;
			}
		}


		// compute filter poles
		let mut pole: Complex<f64>;
		for i in 0..order {
			pole = -Complex{
				re: 0.,
				im: PI * (-(order as i32) + 1 + 2 * (i as i32)) as f64 / (2 * order) as f64
			}.exp();
			pole = Complex{
				re: pole.re * mu.sinh(),
				im: pole.im * mu.cosh()
			};
			pole = 1. / pole;
			output.poles.push(pole);
			comp_gain *= -pole;
		}

		// compute filter gain
		output.gain = comp_gain.re;

		// adapt the poles, zeros and gain values to the cutoff frequency
		output.adapt_frequencies(band);
		output
	}
	
	/// Initializes an allpass filter with a gain of 1.
	/// This method is used to create a custom filter.
	pub fn init_filter(fs: f64) -> Filter {

		Filter {
			gain: 1.,
			poles: Vec::new(),
			zeros: Vec::new(),
			fs: fs
		}
	}
	/// Multiply the gain by given factor
	/// ```math
	/// H(s) = G
	/// ```
	pub fn gain_factor(&mut self, gain: f64) {

		self.gain *= gain;
	}

	/// Add a first order pole at the given cutoff frequency 
	/// ```math
	/// H(s) = \frac{1}{1 + \frac{s}{2 \pi fc}}
	/// ```
	pub fn add_pole_1(&mut self, fc: f64) {
		
		let omega: Complex<f64> = Complex{re: 2. * PI * fc, im: 0.};
		self.poles.push(omega);
		self.gain *= omega.norm();
	}

	/// Add a first order zero at the given cutoff frequency 
	/// ```math 
	/// H(s) = \frac{1}{1 + \frac{s}{q 2 \pi fc} + \left(\frac{s}{2 \pi fc}\right)^2}
	/// ```
	pub fn add_pole_2(&mut self, fc: f64, q: f64) {
		let sum: f64 = PI * fc / q;
		let root: Complex<f64>;
		let eps: f64 = 1e-20;

		if q > 0.5+eps {
			root = Complex{re: 0., im: (4. * q * q - 1.).sqrt()};
		} else if q < 0.5-eps {
			root = Complex{re: (1. - 4. * q * q).sqrt(), im: 0.};
		} else {
			root = Complex{re: 0., im: 0.};
		}
		let s1: Complex<f64> = sum * (-1. - root);
		let s2: Complex<f64> = sum * (-1. + root);
		self.poles.push(s1);
		self.poles.push(s2);
		self.gain *= (s1 * s2).norm();
	}
	
	/// Add a n-order integrator filter
	/// ```math
	/// H(s) = \frac{1}{s^n}
	/// ```
	pub fn add_integrator(&mut self, order: usize) {
		
		for _i in 0..order {
			self.poles.push(Complex{re: 0., im: 0.});
		}

	}

	/// Add a first order zero at the given cutoff frequency
	/// ```math
	/// H(s) = 1 + \frac{s}{2 \pi fc}
	/// ```
	pub fn add_zero_1(&mut self, fc: f64) {
		
		let omega: Complex<f64> = Complex{re: 2. * PI * fc, im: 0.};
		self.poles.push(omega);
		self.gain /= omega.norm();
	}
	
	/// Add a first order zero at the given cutoff frequency
	/// ```math 
	/// H(s) = 1 + \frac{s}{q 2 \pi fc} + \left(\frac{s}{2 \pi fc}\right)^2
	/// ```
	pub fn add_zeros_2(&mut self, fc: f64, q: f64) {
		
		let sum: f64 = PI * fc / q;
		let root: Complex<f64>;
		let eps: f64 = 1e-20;

		if q > 0.5+eps {
			root = Complex{re: 0., im: (4. * q * q - 1.).sqrt()};
		} else if q < 0.5-eps {
			root = Complex{re: (1. - 4. * q * q).sqrt(), im: 0.};
		} else {
			root = Complex{re: 0., im: 0.};
		}
		let s1: Complex<f64> = sum * (-1. - root);
		let s2: Complex<f64> = sum * (-1. + root);
		self.zeros.push(s1);
		self.zeros.push(s2);
		self.gain /= (s1 * s2).norm();
	}

	/// Add a n-order derivator filter
	/// ```math
	/// H(s) = s^n
	/// ```
	pub fn add_derivator(&mut self, order: usize) {
		
		for _i in 0..order {
			self.zeros.push(Complex{re: 0., im: 0.});
		}

	}

	/* ----------------------------------------------------------------------------------------- *
	 * getter functions
	 * ----------------------------------------------------------------------------------------- */

	pub fn get_poles(&self) -> Vec<Complex<f64>> {
		self.poles.clone()
	}
	pub fn get_zeros(&self) -> Vec<Complex<f64>> {
		self.zeros.clone()
	}
	pub fn get_gain(&self) -> f64 {
		self.gain
	}
	pub fn get_fs(&self) -> f64 {
		self.fs
	}

	/* ----------------------------------------------------------------------------------------- *
	 * print functions
	 * ----------------------------------------------------------------------------------------- */

	pub fn print(&self) {
		println!("gain: {}", self.gain);
		println!("zeros: {:#?}", self.zeros);
		println!("poles: {:#?}", self.poles);
	}

}




