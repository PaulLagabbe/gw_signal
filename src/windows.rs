/* --------------------------------------------------------------------------------------------- *
 * Window functions
 * --------------------------------------------------------------------------------------------- */

use rustfft::num_complex::Complex;
use std::f64::consts::PI;

/// Window object iterable
#[derive(Clone)]
pub struct Window {
	size: usize,
	overlap: usize,
	function: fn(usize, usize) -> f64,
}


impl Window {

	/// constructor method, create a window object
	pub fn rectangle(duration: f64, overlap_time: f64, fs: f64) -> Self {

		let size: usize = (duration * fs).round() as usize;
		let overlap: usize = (overlap_time * fs).round() as usize;
		assert!(size > overlap);

		Window {
			size,
			overlap,
			function: rectangle,
		}

	}
	/// constructor method, create a window object
	pub fn hann(duration: f64, overlap_time: f64, fs: f64) -> Self {

		let size: usize = (duration * fs).round() as usize;
		let overlap: usize = (overlap_time * fs).round() as usize;
		assert!(size > overlap);

		Window {
			size,
			overlap,
			function: hann,
		}

	}
	/// Compute the number of window that can be computed on the data vector given in parameter
	pub fn nb_fft(&self, data_size: usize) -> usize {
		(data_size - self.size) / (self.size - self.overlap) + 1
	}

	/// Get size of the window in number of points
	pub fn get_size(&self) -> usize {
		self.size
	}

	/// return a slive of the parameter data vector that has been windowed
	pub fn get_windowed_data(&self, data: Vec<f64>, k: usize) -> Vec<Complex<f64>> {
		
		assert!(k < self.nb_fft(data.len()));

		// initialize output vector
		let mut output: Vec<Complex<f64>> = Vec::new();
		let start: usize = k * (self.size - self.overlap);
		for i in 0..self.size {

			output.push(Complex{
				re:	data[i+start] * (self.function)(self.size, i),
				im: 0.});
		}

		output
	}

	/// Compute the integral of the square window
	pub fn get_norm_factor(&self) -> f64 {

        let mut alpha: f64 = 0.;
        for i in 0..self.size {
            alpha += (self.function)(self.size, i).powi(2);
        };
		alpha
	}
}

/* --------------------------------------------------------------------------------------------- *
 * Window function
 * --------------------------------------------------------------------------------------------- */
fn hann(n: usize, k: usize) -> f64 {
	let phi: f64 = 2. * PI * (k as f64) / (n as f64);
	0.5 * (1. - phi.cos())
}

fn rectangle(_n: usize, _k: usize) -> f64 {
	1.
}

	
