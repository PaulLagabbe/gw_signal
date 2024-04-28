/* --------------------------------------------------------------------------------------------- *
 * Libraries
 * --------------------------------------------------------------------------------------------- */

pub mod timeseries;
pub mod frequencyseries;
pub mod spectrogram;
pub mod filters;
pub mod windows;
//pub mod plot;


use std::{
	str::FromStr,
	string::ToString,
	//f64::consts::PI,
	fmt::Display,
	fs::File,
	io::BufReader,
	io::BufRead,
	io::Write};

use crate::{
	timeseries::TimeSeries,
	frequencyseries::FrequencySeries,
	spectrogram::Spectrogram,
	filters::Filter,
	windows::Window,
//	plot::{RealPlot, ComplexPlot}
};
use rustfft::FftPlanner;
use num::{Complex, complex::ComplexFloat};
use std::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
//use more_asserts as ma;



/* --------------------------------------------------------------------------------------------- *
 * Constructors
 * --------------------------------------------------------------------------------------------- */
/// Spectral analysis methods
impl<D> TimeSeries<D> 
	where D: ComplexFloat + AddAssign + SubAssign + MulAssign + DivAssign,
{

	
	/// Compute the cross spectal density between two signals, using the Welch's method
	/// 
	/// # Example
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	frequencyseries::*,
	/// 	windows::*,
	/// };
	///
	/// // creates two white noise signals
	/// let window: Window = hann(1., 0.5, 1e3);
	/// let mut signal_1: TimeSeries = TimeSeries::white_noise(20000, 1e3, 0f64, 1f64);
	/// let mut signal_2: TimeSeries = signal_1.clone * 2.;
	///
	/// // compute the csd
	/// let csd: FrequencySeries = signal_1.csd(&signal_2, &window);
	/// 
	/// ```
	pub fn csd(
		&self,
		other: &TimeSeries<D>,
		window: &Window) -> FrequencySeries {
	
		assert_eq!(self.get_fs(), other.get_fs());
		// initialize fft
		let mut planner = FftPlanner::new();
		let fft = planner.plan_fft_forward(window.get_size());

		/*
		// compute the mean of each time series
		let mut self_clone = &mut self.clone();
		self_clone = self_clone - self.mean();
		let mut other_clone = &mut other.clone();
		other_clone = other_clone - other.mean();
		*/

		// initialize frequency series
		let mut temp_1: Vec<Complex<f64>>;
		let mut temp_2: Vec<Complex<f64>>;
		let f_max: f64 = self.get_fs() / 2.;
		let mut output: &mut FrequencySeries = &mut FrequencySeries::from_vector(
			f_max, vec![Complex{ re: 0.0f64, im: 0.0f64 }; window.get_size() / 2 + 1]);
		
		let nb_fft: usize = window.nb_fft(self.get_size());
		for i in 0..nb_fft {
			// compute windowed data
			temp_1 = window.get_windowed_data(self.get_data(), i);
			temp_2 = window.get_windowed_data(other.get_data(), i);
			// compute fft
			fft.process(&mut temp_1); fft.process(&mut temp_2);
			// compute product
			output = output + &*(
				&mut FrequencySeries::from_vector(
					f_max, 
					temp_1[0..(window.get_size() / 2 + 1)].to_vec().clone()
				).conj()
				* &FrequencySeries::from_vector(
					f_max,
					temp_2[0..(window.get_size() / 2 + 1)].to_vec().clone()
				)
			);
		}
		output = output / (f_max * window.get_norm_factor() * nb_fft as f64);
		output.clone()

	}
	/// Compute power spectral density using cross spectral density with itself
	/// 
	/// # Example
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	frequencyseries::*,
	/// 	windows::*,
	/// };
	///
	/// // creates two white noise signals
	/// let window: Window = hann(1., 0.5, 1e3);
	/// let mut signal_1: TimeSeries = TimeSeries::white_noise(2000000, 1e3, 0f64, 1f64);
	///
	/// // compute the csd
	/// let psd: FrequencySeries = signal_1.psd(&window);
	/// 
	/// ```
	pub fn psd(
		&self,
		window: &Window) -> FrequencySeries {

		// use csd
		let self_copy: &TimeSeries<D> = &(self.clone());
		self.csd(&self_copy, window)
	}
	
	/// Compute the amplitude spectral density of a signal. Uses the psd function
	/// 
	/// # Example
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	frequencyseries::*,
	/// 	windows::*,
	/// };
	///
	/// // creates two white noise signals
	/// let window: Window = hann(1., 0.5, 1e3);
	/// let mut signal_1: TimeSeries = TimeSeries::white_noise(2000000, 1e3, 0f64, 1f64);
	///
	/// // compute the csd
	/// let asd: FrequencySeries = signal_1.asd(&window);
	/// 
	/// ```
	pub fn asd(
		&self,
		window: &Window) -> FrequencySeries {
		
		self.psd(window).sqrt()
	}
	/// Compute the coherence between two signals.
	/// `\gamma_{1,2}(f) = \frac{|csd_{1,2}(f)|}{psd_1(f) \cdot psd_2(f)}`
	/// 
	/// # Example
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	frequencyseries::*,
	/// 	windows::*,
	/// };
	///
	/// // creates two white noise signals
	/// let window: Window = hann(1., 0.5, 1e3);
	/// let mut signal_1: TimeSeries = TimeSeries::white_noise(2000000, 1e3, 0f64, 1f64);
	/// let mut signal_2: TimeSeries = signal_1.clone() * 2.;
	///
	/// // compute the csd
	/// let coherence: FrequencySeries = signal_1.coherence(&signal_2, &window);
	/// 
	/// ```
	pub fn coherence(
		&self,
		other: &TimeSeries<D>,
		window: &Window) -> FrequencySeries {

		let psd1: FrequencySeries = self.clone().psd(window);
		let psd2: FrequencySeries = other.clone().psd(window);
		let mut csd: FrequencySeries = self.csd(other, window).abs2();
		((&mut csd / &psd1) / &psd2).clone()
	}

	/// Compute the transfer functions between two signals.
	/// `\TF_{1,2}(f) = \frac{csd_{1,2}(f)}{psd_1(f)}`
	/// 
	/// # Example
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	frequencyseries::*,
	/// 	windows::*,
	/// };
	///
	/// // creates two white noise signals
	/// let window: Window = hann(1., 0.5, 1e3);
	/// let mut signal_1: TimeSeries = TimeSeries::white_noise(2000000, 1e3, 0f64, 1f64);
	/// let mut signal_2: TimeSeries = signal_1.clone() * 2.;
	///
	/// // compute the csd
	/// let transfer_function: FrequencySeries = signal_1.transfer_function(&signal_2, &window);
	/// 
	/// ```
	pub fn transfer_function(
		&self,
		other: &TimeSeries<D>,
		window: &Window) -> FrequencySeries {

		let psd: FrequencySeries = self.clone().psd(window);
		let mut csd: FrequencySeries = self.csd(other, window);
		(&mut csd / &psd).clone()
	}

	/* ----------------------------------------------------------------------------------------- *
	 * compute spectrogram
	 * ----------------------------------------------------------------------------------------- */
	/// Compute the cross spectal density between two signals, using the Welch's method
	/// 
	/// # Example
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	frequencyseries::*,
	/// 	windows::*,
	/// };
	///
	/// // creates two white noise signals
	/// let window: Window = hann(1., 0.5, 1e3);
	/// let mut signal_1: TimeSeries = TimeSeries::white_noise(2000000, 1e3, 0f64, 1f64);
	/// let mut signal_2: TimeSeries = signal_1.clone() * 2.;
	///
	/// // compute the csd
	/// let csd: Spectrogram = signal_1.time_csd(&signal_2, &window, 10.);
	/// 
	/// ```
	pub fn time_csd(
		&self,
		other: &TimeSeries<D>,
		window: &Window,
		nb_fft: usize) -> Spectrogram {

		let step: usize = window.get_size() - window.get_overlap();
		let f_max: f64 = self.get_fs() / 2.;
		// check if the number of fft is over 1 and below the maximum number of fft
		assert!(nb_fft > 0);
		assert!(nb_fft <= window.nb_fft(self.get_size()));
		assert_eq!(self.get_fs(), other.get_fs());

		// initialize fft
		let mut planner = FftPlanner::new();
		let fft = planner.plan_fft_forward(window.get_size());

		/*
		// clone the time series and subtract their mean
		let mut self_clone = &mut self.clone();
		self_clone = self_clone - self.mean();
		let mut other_clone = &mut other.clone();
		other_clone = other_clone - other.mean();
		*/
		// compute all frequency series and put them into a vector
		let mut temp_1: Vec<Complex<f64>>;
		let mut temp_2: Vec<Complex<f64>>;

		let mut temp_series: &mut FrequencySeries = &mut FrequencySeries::from_vector(
			f_max, vec![Complex{ re: 0.0f64, im: 0.0f64 }; window.get_size() / 2 + 1]);
		let mut series_vec: Vec<FrequencySeries> = Vec::new();

		for i in 0..window.nb_fft(self.get_size()) {
			temp_series.set_to_zero();
			// compute windowed data
			temp_1 = window.get_windowed_data(self.get_data(), i);
			temp_2 = window.get_windowed_data(other.get_data(), i);
			// compute fft
			fft.process(&mut temp_1); fft.process(&mut temp_2);
			// compute product and add result into the frequency series vector
			temp_series = temp_series +
				&*( &mut FrequencySeries::from_vector(
					f_max, 
					temp_1[0..(window.get_size() / 2 + 1)].to_vec().clone()
				).conj()
				* &FrequencySeries::from_vector(
					f_max,
					temp_2[0..(window.get_size() / 2 + 1)].to_vec().clone()
				) / (f_max * window.get_norm_factor() * nb_fft as f64) );

			series_vec.push(temp_series.clone());

		}

		
		// computes first value and push it into the output vector
		let mut one_series: &mut FrequencySeries = &mut FrequencySeries::from_vector(
			f_max, vec![Complex{ re: 0.0f64, im: 0.0f64 }; window.get_size() / 2 + 1]);
		for i in 0..nb_fft {
			one_series = one_series + &series_vec[i];
		}
		
		// initialize output data vector
		let mut output: Vec<Vec<Complex<f64>>> = Vec::new();
		output.push( one_series.get_data() );

		// roll over the rest of the time series
		for i in nb_fft..window.nb_fft(self.get_size()) {
			one_series = one_series + &series_vec[i];
			one_series = one_series - &series_vec[i-nb_fft];
			output.push( one_series.get_data() );
		}
		
		Spectrogram::from_vector(
			f_max,
			((nb_fft - 1) * step + window.get_size()) as f64 / self.get_fs(),
			step as f64 / self.get_fs(),
			output)
	}
	/// Compute power spectral density using cross spectral density with itself
	/// 
	/// # Example
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	frequencyseries::*,
	/// 	windows::*,
	/// };
	///
	/// // creates two white noise signals
	/// let window: Window = hann(1., 0.5, 1e3);
	/// let mut signal_1: TimeSeries = TimeSeries::white_noise(2000000, 1e3, 0f64, 1f64);
	///
	/// // compute the csd
	/// let psd: FrequencySeries = signal_1.time_psd(&window, 10);
	/// 
	/// ```
	pub fn time_psd(
		&self,
		window: &Window,
		nb_fft: usize) -> Spectrogram {

		// use csd
		let self_copy: &TimeSeries<D> = &(self.clone());
		self.time_csd(&self_copy, window, nb_fft)
	}
	/// Compute the amplitude spectral density of a signal. Uses the psd function
	/// 
	/// # Example
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	frequencyseries::*,
	/// 	windows::*,
	/// };
	///
	/// // creates two white noise signals
	/// let window: Window = hann(1., 0.5, 1e3);
	/// let mut signal_1: TimeSeries = TimeSeries::white_noise(2000000, 1e3, 0f64, 1f64);
	///
	/// // compute the csd
	/// let asd: FrequencySeries = signal_1.time_asd(&window, 10);
	/// 
	/// ```
	pub fn time_asd(
		&self,
		window: &Window,
		nb_fft: usize) -> Spectrogram {
		
		self.time_psd(window, nb_fft).sqrt()
	}
	/// Compute the coherence between two signals.
	/// `\gamma_{1,2}(f) = \frac{|csd_{1,2}(f)|}{psd_1(f) \cdot psd_2(f)}`
	/// 
	/// # Example
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	frequencyseries::*,
	/// 	windows::*,
	/// };
	///
	/// // creates two white noise signals
	/// let window: Window = hann(1., 0.5, 1e3);
	/// let mut signal_1: TimeSeries = TimeSeries::white_noise(2000000, 1e3, 0f64, 1f64);
	/// let mut signal_2: TimeSeries = signal_1.clone() * 2.;
	///
	/// // compute the csd
	/// let coherence: FrequencySeries = signal_1.time_cohe(&signal_2, &window, 10);
	/// 
	/// ```
	pub fn time_cohe(
		&self,
		other: &TimeSeries<D>,
		window: &Window,
		nb_fft: usize) -> Spectrogram {

		let psd1: Spectrogram = self.clone().time_psd(window, nb_fft);
		let psd2: Spectrogram = other.clone().time_psd(window, nb_fft);
		let mut csd: Spectrogram = self.time_csd(other, window, nb_fft).abs2();
		((&mut csd / &psd1) / &psd2).clone()
	}
	/// Compute the transfer functions between two signals.
	/// `\TF_{1,2}(f) = \frac{csd_{1,2}(f)}{psd_1(f)}`
	/// 
	/// # Example
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	frequencyseries::*,
	/// 	windows::*,
	/// };
	///
	/// // creates two white noise signals
	/// let window: Window = hann(1., 0.5, 1e3);
	/// let mut signal_1: TimeSeries = TimeSeries::white_noise(2000000, 1e3, 0f64, 1f64);
	/// let mut signal_2: TimeSeries = signal_1.clone() * 2.;
	///
	/// // compute the csd
	/// let transfer_function: FrequencySeries = signal_1.time_tf(&signal_2, &window 10);
	/// 
	/// ```
	pub fn time_tf(
		&self,
		other: &TimeSeries<D>,
		window: &Window,
		nb_fft: usize) -> Spectrogram {

		let psd: Spectrogram = self.clone().time_psd(window, nb_fft);
		let mut csd: Spectrogram = self.time_csd(other, window, nb_fft);
		(&mut csd / &psd).clone()
	}
}

/* --------------------------------------------------------------------------------------------- */

/// Signal processing methods
impl<D> TimeSeries<D> 
	where D: ComplexFloat,
{

	/// The following method apply an IIR filter to a time series
	/// modify the original time series object.
	/// # Example	
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	filter::*,
	/// };
	///
	/// // creates two white noise signals
	/// let fs: f64 = 1e3;
	/// let mut signal_1: TimeSeries<f64> = TimeSeries::white_noise(20000, fs, 0., 1.);
	/// 
	/// // generates an 8th butterworth lowpass filter at 10 Hz
	/// let butter: Filter::butterworth(8, BType::LowPass(10.), fs);
	/// 
	/// // apply the filter to the signal
	/// let mut signal_2: TimeSeries = signal_1.apply_filter(&butter);
	///
	/// ```
	pub fn apply_filter(&mut self, input_filter: &Filter) {
		
		assert!((1. - self.get_fs() / input_filter.get_fs()).abs() < 1e-10);
		let mut flt: Filter = input_filter.clone();
		// warp frequencies
		flt.adapt_frequencies(true);
		// compute bilinear transform of the filter, the filter is now in the z-space
		flt.bilinear_transform();
		// compute the polynomial coefficiants of the z-transform of the filter
		let (mut b, mut a): (Vec<f64>, Vec<f64>) = flt.polezero_to_coef();

		// complete a or b with 0. so that the two vectors have the same size
		if a.len() < b.len() {
			a.append(&mut vec![0.; b.len()-a.len()]);
		}
		else if a.len() > b.len() {
			b.append(&mut vec![0.; a.len()-b.len()]);
		}
		// number of pole and zeros
		let n: usize = a.len();

		// apply filter to the data vector
		let mut x: Vec<D> = vec![self[0]; n-1];
		x.append(&mut self.get_data());
		let mut y: Vec<D> = x.clone();

		for i in 0..self.get_size() {
			// computed y signal value
			let mut temp_y: D = D::zero();
			for j in 0..b.len() {
				temp_y = temp_y + x[i + b.len()-1 - j] * D::from(b[j]).unwrap();
			}
			for j in 1..a.len() {
				temp_y = temp_y - y[i + a.len()-1 - j] * D::from(a[j]).unwrap();
			}
			temp_y = temp_y / D::from(a[0]).unwrap();
			// apply computed value
			self[i] = temp_y;
			y[i+a.len()-1] = temp_y;
		}
	}


	/// Computes the convolution product of two signals using the fft methods
	/// # Example	
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// };
	///
	/// // creates two white noise signals
	/// let fs: f64 = 1e3;
	/// let signal_1: TimeSeries<f64> = TimeSeries::white_noise(20000, fs, 0., 1.);
	/// let signal_2: TimeSeries<f64> = TimeSeries::white_noise(10000, fs, 0., 1.);
	/// 
	/// // computes convolution
	/// let signal: TimeSeries<f64> = signal_1.convolution(&signal_2);
	/// 
	/// ```
	pub fn convolution(&self, other: &TimeSeries<D>) -> TimeSeries<D> {
		
		// test if the sampling frequencies of the two time series are the same
		assert_eq!(self.get_fs(), other.get_fs());
		let fs: f64 = self.get_fs();
		let t0: f64 = 0.; // starting time set to 0
		// clone data vectors with complex data type
		let complex_self = self.to_c64();
		let complex_other = other.to_c64();

		// zero padding the data vectors
		let mut x: Vec<Complex<f64>> = complex_self.get_data();
		x.append(&mut vec![Complex{re: 0., im: 0.}; other.get_size()-1]);
		let mut y: Vec<Complex<f64>> = complex_other.get_data();
		y.append(&mut vec![Complex{re: 0., im: 0.}; self.get_size()-1]);
		let n: usize = x.len();

		// computes the fft
		let mut planner = FftPlanner::new();
		let fft = planner.plan_fft_forward(x.len());
		fft.process(&mut x); fft.process(&mut y);

		// multply the fft
		for it in y.iter().zip(x.iter_mut()) {
			*it.1 *= *it.0 / Complex{re: n as f64, im: 0.};
		}

		// computes the inverse fft
		planner = FftPlanner::new();
		let ifft = planner.plan_fft_inverse(x.len());
		ifft.process(&mut x);
		// test if the data value is complex or not
		let mut sum: Complex<f64> = Complex::from(0f64);
		for val in x.iter() {
			sum = sum + *val;
		}
		// convert data into data type
		let mut data: Vec<D> = Vec::new();
		for val in x.iter() {
			if sum.im().abs() < sum.re().abs() * 1e-10 {
				data.push(D::from((*val).re).unwrap());
			} else {
				data.push(D::from(*val).unwrap());
			}
		}

		// create data vector
		TimeSeries::from_vector(fs, t0, data)
	}
}


/* --------------------------------------------------------------------------------------------- *
 * SeriesIO trait 
 * --------------------------------------------------------------------------------------------- */

/// This trait is for dedug purpose only.
/// It provides a function to print some samples of the time/frequency series and to print it into a csv file.
pub trait SeriesIO {
	/// 
	/// # Example	
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	Series_IO,
	/// };
	///
	/// // creates a white noise signals
	/// let fs: f64 = 1e3;
	/// let mut signal: TimeSeries = TimeSeries::white_noise(20000, fs, 1.);
	/// 
	/// // print the 10 first values of the time series
	/// signal.print(0, 10);
	///
	/// ```
	fn print(&self, n1: usize, n2: usize);

	/// # Example	
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	filter::*,
	/// };
	/// 
	/// // creates a white noise signals
	/// let fs: f64 = 1e3;
	/// let mut signal: TimeSeries = TimeSeries::white_noise(20000, fs, 1.);
	/// 
	/// // Write csv files
	/// signal.write_csv("TimeSeries.csv");
	/// ```
	fn write_csv(&self, file_name: &str);

	/// # Example	
	/// ```
	/// use gw_signal::{
	/// 	timeseries::*,
	/// 	filter::*,
	/// };
	/// 
	/// // Read time series from csv files
	/// let mut signal: TimeSeries<f64> = TimeSeries::<f64>::read_csv("TimeSeries.csv");
	/// ```
	fn read_csv(file_name: &str) -> Self;
}


impl<D: ComplexFloat + ToString + FromStr + Display> SeriesIO for TimeSeries<D> {
	
	fn print(&self, n1: usize, n2: usize) {
		let mut time: f64;		 
		for i in n1..n2 {
			// compute time
			time = self.get_t0() + (i as f64) / self.get_fs();
			println!("t = {:.3} s: {:.6}", time, self[i]);
		}
	}

	fn write_csv(&self, file_name: &str) {
		let mut w = File::create(file_name).unwrap();
		writeln!(&mut w, "time,value").unwrap();
		let mut time: f64 = self.get_t0();
		for value in self.get_data().iter() {
			// compute time
			time += 1f64 / self.get_fs();
			writeln!(&mut w, "{},{}", time, value.to_string()).unwrap();
		}
	}

	fn read_csv(file_name: &str) -> Self {
		println!("Read file: {}", file_name);
		// read file
		let r = File::open(file_name).expect("The file is not found!");
		let buffer = BufReader::new(r);
		// initialize time and data vectors
		let (mut time, mut data): (Vec<f64>, Vec<D>) = (Vec::new(), Vec::new());
		// make iterator over lines and read file header
		let mut line_iter = buffer.lines();
		// read line, split it over the "," character and make a vector of strings
		let mut line_str = line_iter.next().unwrap().unwrap();
		let mut line_vec: Vec<&str> = line_str.split(",").collect();
		assert_eq!(line_vec[0], "time");
		for line in line_iter {
			line_str = line.expect("Unable to read line");
			line_vec = line_str.split(",").collect();

			time.push(f64::from_str(&line_vec[0]).expect("Unable to read time value"));
			let read_data = D::from_str(&line_vec[1]);
			match read_data {
				Ok(x) => data.push(x),
				Err(_) => panic!("Unable to read data value"),
			}
		}
		let frequency: f64 = (time.len()-1) as f64 / (time[time.len()-1] - time[0]);
		TimeSeries::from_vector(frequency, time[0], data)
	}

}



impl SeriesIO for FrequencySeries {
	
	fn print(&self, n1: usize, n2: usize) {
		let mut freq: f64;
		for i in n1..n2 {
			// compute time
			freq = self.get_f_max() * (i as f64) / ((self.get_size()-1) as f64);
			println!("f = {:.3} Hz: {:.6} + {:.6}i", freq, self[i].re, self[i].im);
		}
	}
	
	fn write_csv(&self, file_name: &str) {
		let mut w = File::create(file_name).unwrap();
		writeln!(&mut w, "frequency,value").unwrap();
		let mut freq: f64 = 0f64;
		for value in self.get_data().iter() {
			// compute time
			writeln!(&mut w, "{},{}", freq, value.to_string()).unwrap();
			freq += self.get_f_max() / ((self.get_size()-1) as f64);
		}
	}

	fn read_csv(file_name: &str) -> Self {
		println!("Read file: {}", file_name);
		// read file
		let r = File::open(file_name).expect("The file is not found!");
		let buffer = BufReader::new(r);
		// initialize time and data vectors
		let (mut freq, mut data): (Vec<f64>, Vec<Complex<f64>>) = (Vec::new(), Vec::new());
		// make iterator over lines and read file header
		let mut line_iter = buffer.lines();
		// read line, split it over the "," character and make a vector of strings
		let mut line_str = line_iter.next().unwrap().unwrap();
		let mut line_vec: Vec<&str> = line_str.split(",").collect();
		assert_eq!(line_vec[0], "frequency");
		for line in line_iter {
			line_str = line.expect("Unable to read line");
			line_vec = line_str.split(",").collect();
			freq.push(f64::from_str(&line_vec[0]).expect("Unable to read frequency value"));
			let read_data = Complex::from_str(&line_vec[1]);
			match read_data {
				Ok(x) => data.push(x),
				Err(_) => panic!("Unable to read data value"),
			}
		}
		FrequencySeries::from_vector(freq[freq.len()-1], data)
	}
}

impl SeriesIO for Spectrogram {
	
	fn print(&self, n1: usize, n2: usize){
		println!("Print the first frequency series of the spectrogram");
		let frequency_series = &self[0];
		let mut freq: f64;
		for i in n1..n2 {
			// compute time
			freq = self.get_f_max() * (i as f64) / ((self.get_size().1-1) as f64);
			println!("f = {:.3} Hz: {:.6} + {:.6}i",
				freq, frequency_series[i].re, frequency_series[i].im);
		}
	}

	fn write_csv(&self, file_name: &str) {
		let (_size_time, size_freq) = self.get_size();
		// create file
		let mut w = File::create(file_name).unwrap();
		// write first line with frequency values
		let mut line = String::from("time\\frequency");
		let mut freq: f64;
		for i in 0..size_freq {
			// compute frequency
			freq = self.get_f_max() * (i as f64) / (size_freq - 1) as f64;
			line.push_str(",");
			line.push_str(&freq.to_string());
		}
		//line.push_str("\n");
		// write firsst line
		writeln!(&mut w, "{}", line).unwrap();
		// write data vector
		let mut time: f64 = self.get_t0();
		for frequency_series in self.get_data().iter() {
			// compute time
			line = time.to_string();
			for value in frequency_series.iter() {
				line.push_str(",");
				line.push_str(&value.to_string());
			}
			//line.push_str("\n");
			writeln!(&mut w, "{}", line).unwrap();
			time += self.get_dt();
		}
	}

	fn read_csv(file_name: &str) -> Self {
		println!("Read file: {}", file_name);
		// read file
		let r = File::open(file_name).expect("The file is not found!");
		let buffer = BufReader::new(r);
		
		// make iterator over lines and read file header
		let mut line_iter = buffer.lines();
		
		// read line, split it over the "," character and make a vector of strings
		let mut line_str = line_iter.next().unwrap().unwrap();
		let mut line_vec: Vec<&str> = line_str.split(",").collect();
		
		// read last value of the frequency vector
		let f_max: f64 = f64::from_str(line_vec[line_vec.len() - 1])
			.expect("Unable to read maximum frequency");
		// initialize data vector and time vector
		let mut time: Vec<f64> = Vec::new();
		let mut data: Vec<Vec<Complex<f64>>> = Vec::new();
		let mut is_time: bool;

		// fill the vector
		for line in line_iter {
			line_str = line.expect("Unable to read line");
			line_vec = line_str.split(",").collect();
			is_time = true;
			let mut line_vector: Vec<Complex<f64>> = Vec::new();
			for str_value in line_vec.iter() {
				if is_time {
					time.push(f64::from_str(str_value).expect("unable to read time value"));
				} else {
					let read_data = Complex::from_str(str_value);
					match read_data {
						Ok(x) => line_vector.push(x),
						Err(_) => panic!("Unable to read data value"),
					}
				}
				is_time = false;
			}
			data.push(line_vector);
		}
		Spectrogram::from_vector(f_max, time[0], time[1]-time[0], data)
	}
}


/* --------------------------------------------------------------------------------------------- *
 * Plot time series
 * --------------------------------------------------------------------------------------------- */
/*
impl RealPlot {
	/// Add one time series to the plot
	pub fn add_timeseries<D, F>(&mut self, series: &TimeSeries<D>) 
		where D: ComplexFloat<Real = F>, F: Float,
	{
		
		// build time axis
		let mut time: Vec<f64> = Vec::new();
		for i in 0..series.get_size() {
			time.push(series.get_t0() + (i as f64) / series.get_fs());
		}
		self.add_data_vector(time, series.to_f64().get_data());
	}

	pub fn add_frequencyseries(&mut self, series: &FrequencySeries) {
		
		// build time axis
		let mut time: Vec<f64> = Vec::new();
		let mut real_data: Vec<f64> = Vec::new();
		for i in 1..series.get_size() {
			time.push((i as f64) / ((series.get_size() - 1) as f64) * series.get_f_max());
			real_data.push(series[i].re);
		}
		self.add_data_vector(time, real_data);
		self.set_x_scale_to_log(true);
	}

}
impl ComplexPlot {

	/// Add one time series to the plot
	pub fn add_timeseries<D, F>(&mut self, series: &TimeSeries<D>) 
		where D: ComplexFloat<Real = F>, F: Float,
	{
		
		// build time axis
		let mut time: Vec<f64> = Vec::new();
		for i in 0..series.get_size() {
			time.push(series.get_t0() + (i as f64) / series.get_fs());
		}
		self.add_data_vector(time, series.to_c64().get_data());
	}

	pub fn add_frequencyseries(&mut self, series: &FrequencySeries) {
		
		// build time axis
		let mut time: Vec<f64> = Vec::new();
		let mut data: Vec<Complex<f64>> = Vec::new();
		for i in 1..series.get_size() {
			time.push((i as f64) / ((series.get_size() - 1) as f64) * series.get_f_max());
			data.push(series[i]);
		}
		self.add_data_vector(time, data);
		self.set_x_scale_to_log(true);
		self.set_y_scale_to_log(true);
	}

}
*/
/* --------------------------------------------------------------------------------------------- *
 * Filter methods
 * --------------------------------------------------------------------------------------------- */

/// 
impl Filter {
	/// Compute frequency response of the filter
	pub fn frequency_response(&self, size: usize) -> FrequencySeries {

		// warp frequencies
		let mut clone_filter = self.clone();
		clone_filter.adapt_frequencies(false);
		// initialize frequency series
		let mut response: FrequencySeries = FrequencySeries::from_vector(
			clone_filter.get_fs()/2., vec![Complex{re: clone_filter.get_gain(), im: 0.}; size]);
		let mut frequency: f64;

		for i in 0..size {
			frequency = clone_filter.get_fs() / 2. * (i as f64) / ((size-1) as f64);

			// apply zeros
			for z in clone_filter.get_zeros().iter() {
				response[i] *= Complex{re: 0., im: frequency} - z
			}
			// apply poles
			for p in clone_filter.get_poles().iter() {
				response[i] /= Complex{re: 0., im: frequency} - p
			}
		}
		response
	}

}



