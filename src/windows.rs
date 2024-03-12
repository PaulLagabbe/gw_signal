/* --------------------------------------------------------------------------------------------- *
 * Window functions
 * --------------------------------------------------------------------------------------------- */

use rustfft::num_complex::Complex;
use std::f64::consts::PI;

/// Window object iterable
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
		for i in 0..size {
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

	
