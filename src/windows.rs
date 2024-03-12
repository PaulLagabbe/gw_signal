/* --------------------------------------------------------------------------------------------- *
 * Window functions
 * --------------------------------------------------------------------------------------------- */

use rustfft::num_complex::Complex;
use std::f64::consts::PI;

/// Window object, to be used to compute the Fourrier transform of a signal
/// 
/// # Example
/// 
/// ```
/// // libraries
/// use gw_signal::{
///		timeseries as ts,
/// 	frequencyseries as fs,
/// 	window as win,
/// }
/// // sampling frequency in Hz
/// let frequency: f64 = 1e3; 
/// // define hann window of 1 sec duration overlapping over 0.5 sec.
/// let hann_window: win::Window = win::Window::hann(1., 0.5, frequency);
/// // generate noise signal with 20000 samples, with an ampliture of 0.1.
/// let signal: ts::TimeSeries::white_noise(20000, frequency, 0.1);
///
/// // compute power spectral density
/// let response = signal.psd(&hann_window);
/// ```
#[derive(Clone)]
pub struct Window {
	overlap: usize,
	vector: Vec<f64>,
}


impl Window {

	/// constructor method, create a rectangle window object
	pub fn rectangle(duration: f64, overlap_time: f64, fs: f64) -> Self {

		let size: usize = (duration * fs).round() as usize;
		let overlap: usize = (overlap_time * fs).round() as usize;
		assert!(size > overlap);
		// make window vector
		let mut vector: Vec<f64> = Vec::new();
		for _i in 0..size {
			vector.push(1.);
		}
		Window {
			overlap,
			vector,
		}

	}
	/// constructor method, create a Hann window object
	pub fn hann(duration: f64, overlap_time: f64, fs: f64) -> Self {

		let size: usize = (duration * fs).round() as usize;
		let overlap: usize = (overlap_time * fs).round() as usize;
		assert!(size > overlap);
		// make window vector
		let mut vector: Vec<f64> = Vec::new();
		let mut phi: f64;
		for i in 0..size {
			phi =  2. * PI * (i as f64) / (size as f64);
			vector.push(0.5 * (1. - phi.cos()));
		}

		Window {
			overlap,
			vector,
		}

	}
	/// constructor method, create a Blackman window object
	pub fn blackman(duration: f64, overlap_time: f64, fs: f64) -> Self {

		let size: usize = (duration * fs).round() as usize;
		let overlap: usize = (overlap_time * fs).round() as usize;
		assert!(size > overlap);
		// make window vector
		let mut vector: Vec<f64> = Vec::new();
		let mut phi: f64;
		for i in 0..size {
			phi =  2. * PI * (i as f64) / (size as f64);
			vector.push(0.42 - 0.5 * phi.cos() + 0.08 * (2. * phi).cos());
		}

		Window {
			overlap,
			vector,
		}

	}
	/// constructor method, create a Hamming window object
	pub fn hamming(duration: f64, overlap_time: f64, fs: f64) -> Self {

		let size: usize = (duration * fs).round() as usize;
		let overlap: usize = (overlap_time * fs).round() as usize;
		assert!(size > overlap);
		// make window vector
		let mut vector: Vec<f64> = Vec::new();
		let mut phi: f64;
		for i in 0..size {
			phi =  2. * PI * (i as f64) / ((size-1) as f64);
			vector.push(0.54 - 0.46 * phi.cos());
		}

		Window {
			overlap,
			vector,
		}

	}
	/// constructor method, create a Hamming window object
	pub fn gaussian(sigma: f64, duration: f64, overlap_time: f64, fs: f64) -> Self {

		let size: usize = (duration * fs).round() as usize;
		let overlap: usize = (overlap_time * fs).round() as usize;
		assert!(size > overlap);
		// make window vector
		let mut vector: Vec<f64> = Vec::new();
		let mut phi: f64;
		for i in 0..size {
			phi = -0.5 * ((i - size/2) as f64 / sigma).powi(2);
			vector.push(phi.exp());
		}

		Window {
			overlap,
			vector,
		}

	}
	/// Compute the number of window that can be computed on the data vector given in parameter
	pub fn nb_fft(&self, data_size: usize) -> usize {
		(data_size - self.vector.len()) / (self.vector.len() - self.overlap) + 1
	}

	/// Get size of the window in number of points
	pub fn get_size(&self) -> usize {
		self.vector.len()
	}

	/// return a slive of the parameter data vector that has been windowed
	pub fn get_windowed_data(&self, data: Vec<f64>, k: usize) -> Vec<Complex<f64>> {
		
		assert!(k < self.nb_fft(data.len()));

		// initialize output vector
		let mut output: Vec<Complex<f64>> = Vec::new();
		let start: usize = k * (self.vector.len() - self.overlap);
		for i in 0..self.vector.len() {

			output.push(Complex{
				re:	data[i+start] * self.vector[i],
				im: 0.});
		}

		output
	}

	/// Compute the integral of the square window
	pub fn get_norm_factor(&self) -> f64 {

        let mut alpha: f64 = 0.;
        for i in 0..self.vector.len() {
            alpha += self.vector[i].powi(2);
        };
		alpha
	}
}

	
