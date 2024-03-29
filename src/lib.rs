/* --------------------------------------------------------------------------------------------- *
 * Libraries
 * --------------------------------------------------------------------------------------------- */

pub mod timeseries;
pub mod frequencyseries;
pub mod filters;
pub mod windows;

use std::{
	str::FromStr,
	f64::consts::PI,
	fmt::Display,
	fs::File,
	io::BufReader,
	io::BufRead,
	io::Write};

use crate::{
	timeseries::TimeSeries,
	timeseries::data::Data,
	frequencyseries::FrequencySeries,
	filters::Filter,
	windows::Window,
};
use rustfft::FftPlanner;
use num::{Complex, complex::ComplexFloat, NumCast};
use std::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use plotters::prelude::*;
//use more_asserts as ma;



/* --------------------------------------------------------------------------------------------- *
 * Constructors
 * --------------------------------------------------------------------------------------------- */

impl<D> TimeSeries<D> 
	where D: ComplexFloat + AddAssign + SubAssign + MulAssign + DivAssign + Data,
{

	
	/// Compute the cross spectal density between two signals, using the Welch's method
	/// 
	/// # Example
	/// ```
	/// use gw_signal::{
	/// 	timeseries as ts,
	/// 	windows as win,
	/// };
	///
	/// // creates two white noise signals
	/// let mut signal_1: ts::TimeSeries = ts::TimeSeries::white_noise(20000, 1e3, 1.);
	/// let mut signal_2: ts::TimeSeries = signal_1.clone * 2.;
	///
	/// // compute the csd
	/// let csv: fs::FrequencySeries = signal_1.csd(signal_2, 1., 0.5, win::hann);
	/// 
	/// ```
	pub fn csd(
		&self,
		other: &TimeSeries<D>,
		window: &Window) -> FrequencySeries {
	
		// initialize fft
		let mut planner = FftPlanner::new();
		let fft = planner.plan_fft_forward(window.get_size());

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
	/// 	timeseries as ts,
	/// 	windows as win,
	/// };
	///
	/// // creates two white noise signals
	/// let mut signal_1: ts::TimeSeries = ts::TimeSeries::white_noise(20000, 1e3, 1.);
	///
	/// // compute the csd
	/// let psv: fs::FrequencySeries = signal_1.psd(1., 0.5, win::hann);
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
	/// 	timeseries as ts,
	/// 	windows as win,
	/// };
	///
	/// // creates two white noise signals
	/// let mut signal_1: ts::TimeSeries = ts::TimeSeries::white_noise(20000, 1e3, 1.);
	///
	/// // compute the csd
	/// let asd: fs::FrequencySeries = signal_1.asd(1., 0.5, win::hann);
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
	/// 	timeseries as ts,
	/// 	windows as win,
	/// };
	///
	/// // creates two white noise signals
	/// let mut signal_1: ts::TimeSeries = ts::TimeSeries::white_noise(20000, 1e3, 1.);
	/// let mut signal_2: ts::TimeSeries = signal_1.clone() * 2.;
	///
	/// // compute the csd
	/// let coherence: fs::FrequencySeries = signal_1.coherence(signal_2, 1., 0.5, win::hann);
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
	/// 	timeseries as ts,
	/// 	windows as win,
	/// };
	///
	/// // creates two white noise signals
	/// let mut signal_1: ts::TimeSeries = ts::TimeSeries::white_noise(20000, 1e3, 1.);
	/// let mut signal_2: ts::TimeSeries = signal_1.clone() * 2.;
	///
	/// // compute the csd
	/// let transfer_function: fs::FrequencySeries = signal_1.transfer_function(signal_2, 1., 0.5, win::hann);
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


/* --------------------------------------------------------------------------------------------- */
	/// The following method apply an IIR filter to a time series
	/// modify the original time series object.
	/// # Example	
	/// ```
	/// use gw_signal::{
	/// 	timeseries as ts,
	/// 	filter as flt,
	/// };
	///
	/// // creates two white noise signals
	/// let fs: f64 = 1e3;
	/// let mut signal_1: ts::TimeSeries = ts::TimeSeries::white_noise(20000, f64, 1.);
	/// 
	/// // generates an 8th butterworth lowpass filter at 10 Hz
	/// let butter: flt::Filter::butterworth(8, flt::BType::LowType(10.), fs);
	/// 
	/// // apply the filter to the signal
	/// let mut signal_2: ts::TimeSeries = signal_1.apply_filter(butter);
	///
	/// ```

	pub fn apply_filter(&mut self, mut flt: Filter) {
		
		assert!((1. - self.get_fs() / flt.get_fs()).abs() < 1e-10);
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
			let mut temp_y: D = NumCast::from(0).unwrap();
			for j in 0..b.len() {
				temp_y += x[i + b.len()-1 - j].scale(b[j]);
			}
			for j in 1..a.len() {
				temp_y -= y[i + a.len()-1 - j].scale(a[j]);
			}
			temp_y = temp_y.scale(1./a[0]);
			self[i] = temp_y;
			y[i+a.len()-1] = temp_y;
		}

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
	/// 	timeseries as ts,
	/// 	Series_IO,
	/// };
	///
	/// // creates a white noise signals
	/// let fs: f64 = 1e3;
	/// let mut signal: ts::TimeSeries = ts::TimeSeries::white_noise(20000, fs, 1.);
	/// 
	/// // print the 10 first values of the time series
	/// signal.print(0, 10);
	///
	/// ```
    fn print(&self, n1: usize, n2: usize);

	/// # Example	
	/// ```
	/// use gw_signal::{
	/// 	timeseries as ts,
	/// 	filter as flt,
	/// };
	/// 
	/// // creates a white noise signals
	/// let fs: f64 = 1e3;
	/// let mut signal: ts::TimeSeries = ts::TimeSeries::white_noise(20000, fs, 1.);
	/// 
	/// // Write csv files
	/// signal.write_csv("TimeSeries.csv");
	/// ```
	fn write_csv(&self, file_name: &str);

	/// # Example	
	/// ```
	/// use gw_signal::{
	/// 	timeseries as ts,
	/// 	filter as flt,
	/// };
	/// 
	/// // Read time series from csv files
	/// let mut signal: ts::TimeSeries<f64> = ts::TimeSeries::<f64>::read_csv("TimeSeries.csv");
	/// ```
	fn read_csv(file_name: &str) -> Self;
}


impl<D: ComplexFloat + Display + Data> SeriesIO for TimeSeries<D> {
    fn print(&self, n1: usize, n2: usize){
		let mut time: f64;         
		for i in n1..n2 {
            // compute time
            time = self.get_t0() + (i as f64) / self.get_fs();
            println!("t = {:.3} s: {:.6}", time, self[i]);
        }
    }

	fn write_csv(&self, file_name: &str){

		let mut w = File::create(file_name).unwrap();
		writeln!(&mut w, "time,value").unwrap();
		let mut time: f64;
        for i in 0..self.get_size() {
            // compute time
            time = self.get_t0() + (i as f64) / self.get_fs();
            writeln!(&mut w, "{}", self[i].to_str(time)).unwrap();
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
		let line_str = line_iter.next().unwrap().unwrap();
		let line_vec: Vec<&str> = line_str.split(",").collect();
		assert_eq!(line_vec[0], "time");

		for line in line_iter {
			let sample = D::from_string(line.expect("Unable to read line")).unwrap();
			time.push(sample.0);
			data.push(sample.1);
		}
		let frequency: f64 = (time.len()-1) as f64 / (time[time.len()-1] - time[0]);
		TimeSeries::from_vector(frequency, time[0], data)
	}

}



impl SeriesIO for FrequencySeries {
    
	fn print(&self, n1: usize, n2: usize){
        let mut freq: f64;

        for i in n1..n2 {
            // compute time
            freq = self.get_f_max() * (i as f64) / ((self.get_size()-1) as f64);
            println!("f = {:.3} Hz: {:.6} + {:.6}i", freq, self[i].re, self[i].im);
        }
    }
	
	fn write_csv(&self, file_name: &str){

		let mut w = File::create(file_name).unwrap();
		writeln!(&mut w, "frequency,real,imaginary").unwrap();
		let mut freq: f64;

        for i in 0..self.get_size() {
            // compute time
            freq = self.get_f_max() * (i as f64) / ((self.get_size()-1) as f64);
            writeln!(&mut w, "{},{},{}", freq, self[i].re(), self[i].im()).unwrap();
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
			line_str = line.unwrap();
			line_vec = line_str.split(",").collect();
			freq.push(f64::from_str(line_vec[0]).unwrap());
			data.push(Complex{
				re: f64::from_str(line_vec[1]).unwrap(),
				im: f64::from_str(line_vec[2]).unwrap()
			});
		}
		FrequencySeries::from_vector(freq[freq.len()-1], data)
	}

}


/* --------------------------------------------------------------------------------------------- *
 * Plot time series
 * --------------------------------------------------------------------------------------------- */

/// Draw time and frequency series
/// If the series is real, draw one canvas
/// If the series is complex, the canvas is splitted into the modulus and phase
pub trait Plot {
	/// Plot time series
	/// 
	/// # Example
	/// ```
	/// use gw_signal::{
	/// 	timeseries as ts,
	/// 	Plot
	/// };
	/// 
	/// let sampling: f64 = 1e3;
	/// let mut signal: ts::TimeSeries = ts::TimeSeries::white_noise(20000, sampling, 0f64, 1f64);
	/// signal.plot("signal.png");
	/// ```
	fn plot(&self, name: &str);
}

impl Plot for TimeSeries<f32> {
	fn plot(&self, name: &str) {
		
		let (x_min, x_max) = (self.get_t0(), self.get_t0() + self.get_size() as f64 / self.get_fs());
		let (y_min, y_max) = (self.real().min(), self.real().max());
		
		// define white canvas
		let drawing_area = BitMapBackend::new(name, (1365, 768)).into_drawing_area();
		drawing_area.fill(&WHITE).unwrap();
		
		// define canvas features
		let mut ctx = ChartBuilder::on(&drawing_area)
			.margin(20)
			// axis
			.set_label_area_size(LabelAreaPosition::Left, 50)
			.set_label_area_size(LabelAreaPosition::Bottom, 50)
			// grid
			.build_cartesian_2d(x_max..x_min, y_min..y_max)
			.unwrap();

		ctx.configure_mesh().draw().unwrap();
		
		ctx.draw_series(LineSeries::new( (0..self.get_size()).map( |index| (
			self.get_t0() + index as f64 / self.get_fs(),
			self[index]
		)), &RED )).unwrap();

	}
}
impl Plot for TimeSeries<f64> {
	fn plot(&self, name: &str) {
		
		let (x_min, x_max) = (self.get_t0(), self.get_t0() + self.get_size() as f64 / self.get_fs());
		let (y_min, y_max) = (self.real().min(), self.real().max());
	
		// define white canvas
		let drawing_area = BitMapBackend::new(name, (1365, 768)).into_drawing_area();
		drawing_area.fill(&WHITE).unwrap();
		
		// define canvas features
		let mut ctx = ChartBuilder::on(&drawing_area)
			.margin(20)
			// axis
			.set_label_area_size(LabelAreaPosition::Left, 100)
			.set_label_area_size(LabelAreaPosition::Bottom, 100)
			// caption
			.caption("Timeseries",("sans-serif", 40))
			// grid
			.build_cartesian_2d(x_min..x_max, y_min..y_max)
			.unwrap();

		ctx.configure_mesh().draw().unwrap();
		
		ctx.draw_series(LineSeries::new( (0..self.get_size()).map( |index| (
			self.get_t0() + index as f64 / self.get_fs(),
			self[index]
		)), &RED )).unwrap();

	}
}
impl Plot for TimeSeries<Complex<f32>> {
	fn plot(&self, name: &str) {
		
		// define white canvas
		let drawing_area = BitMapBackend::new(name, (1365, 768)).into_drawing_area();
		drawing_area.fill(&WHITE).unwrap();
		let (top, bottom) = drawing_area.split_vertically(128);
		
		// draw real part
		let (x_min, x_max) = (self.get_t0(), self.get_t0() + self.get_size() as f64 / self.get_fs());
		let (y_min, y_max) = (self.real().min(), self.real().max());
		
		// draw canvas
		let mut ctx = ChartBuilder::on(&top)
			.margin(20)
			// axis
			.set_label_area_size(LabelAreaPosition::Left, 64)
			.set_label_area_size(LabelAreaPosition::Bottom, 64)
			// grid
			.build_cartesian_2d(x_max..x_min, y_min..y_max)
			.unwrap();

		ctx.configure_mesh().draw().unwrap();
		
		ctx.draw_series(LineSeries::new( (0..self.get_size()).map( |index| (
			self.get_t0() + index as f64 / self.get_fs(),
			self.real()[index]
		)), &RED )).unwrap();

		// draw imaginary part
		let (y_min, y_max) = (self.imag().min(), self.imag().max());
		
		// draw canvas
		let mut ctx = ChartBuilder::on(&bottom)
			.margin(20)
			// axis
			.set_label_area_size(LabelAreaPosition::Left, 64)
			.set_label_area_size(LabelAreaPosition::Bottom, 64)
			// grid
			.build_cartesian_2d(x_max..x_min, y_min..y_max)
			.unwrap();

		ctx.configure_mesh().draw().unwrap();
		
		ctx.draw_series(LineSeries::new( (0..self.get_size()).map( |index| (
			self.get_t0() + index as f64 / self.get_fs(),
			self.imag()[index]
		)), &RED )).unwrap();

	}
}
impl Plot for TimeSeries<Complex<f64>> {
	fn plot(&self, name: &str) {
		
		// define white canvas
		let drawing_area = BitMapBackend::new(name, (1365, 768)).into_drawing_area();
		drawing_area.fill(&WHITE).unwrap();
		let (top, bottom) = drawing_area.split_vertically(384);
		
		// draw real part
		let (x_min, x_max) = (self.get_t0(), self.get_t0() + self.get_size() as f64 / self.get_fs());
		let (y_min, y_max) = (self.abs().min(), self.abs().max());
		

		// draw canvas
		let mut ctx = ChartBuilder::on(&top)
			.margin(20)
			// axis
			.set_label_area_size(LabelAreaPosition::Left, 64)
			.set_label_area_size(LabelAreaPosition::Bottom, 64)
			// grid
			.build_cartesian_2d(x_max..x_min, y_min..y_max)
			.unwrap();

		ctx.configure_mesh().draw().unwrap();
		
		ctx.draw_series(LineSeries::new( (0..self.get_size()).map( |index| (
			self.get_t0() + index as f64 / self.get_fs(),
			self.abs()[index]
		)), &RED )).unwrap();

		// draw imaginary part
		let (y_min, y_max) = (self.arg().min(), self.arg().max());
		
		// draw canvas
		let mut ctx = ChartBuilder::on(&bottom)
			.margin(20)
			// axis
			.set_label_area_size(LabelAreaPosition::Left, 64)
			.set_label_area_size(LabelAreaPosition::Bottom, 64)
			// grid
			.build_cartesian_2d(x_max..x_min, y_min..y_max)
			.unwrap();

		ctx.configure_mesh().draw().unwrap();
		
		ctx.draw_series(LineSeries::new( (0..self.get_size()).map( |index| (
			self.get_t0() + index as f64 / self.get_fs(),
			self.arg()[index]
		)), &RED )).unwrap();

	}
}

impl Plot for FrequencySeries {
	fn plot(&self, name: &str) {
		
		// define white canvas
		let drawing_area = BitMapBackend::new(name, (1365, 768)).into_drawing_area();
		drawing_area.fill(&WHITE).unwrap();
		let (top, bottom) = drawing_area.split_vertically(384);
		
		// draw real part
		let (x_min, x_max) = (self.get_f_max() / (self.get_size() - 1) as f64, self.get_f_max());
		let (y_min, y_max) = (self.min_abs(), self.max_abs());
		

		// draw canvas
		let mut ctx = ChartBuilder::on(&top)
			.margin(20)
			// axis
			.set_label_area_size(LabelAreaPosition::Left, 64)
			.set_label_area_size(LabelAreaPosition::Bottom, 64)
			// grid
			.build_cartesian_2d((x_min..x_max).log_scale(), (y_min..y_max).log_scale())
			.unwrap();

		ctx.configure_mesh().draw().unwrap();
		
		ctx.draw_series(LineSeries::new( (1..self.get_size()).map( |index| (
			(index as f64 / (self.get_size() - 1) as f64) * self.get_f_max(),
			self.clone().abs()[index].re
		)), &RED )).unwrap();

		// draw imaginary part
		let (y_min, y_max) = (-PI, PI);
		
		// draw canvas
		let mut ctx = ChartBuilder::on(&bottom)
			.margin(20)
			// axis
			.set_label_area_size(LabelAreaPosition::Left, 64)
			.set_label_area_size(LabelAreaPosition::Bottom, 64)
			// grid
			.build_cartesian_2d((x_min..x_max).log_scale(), y_min..y_max)
			.unwrap();

		ctx.configure_mesh().draw().unwrap();
		
		ctx.draw_series(LineSeries::new( (1..self.get_size()).map( |index| (
			(index as f64 / (self.get_size()-1) as f64) * self.get_f_max(),
			self.clone().arg()[index].re
		)), &RED )).unwrap();

	}
}


/* --------------------------------------------------------------------------------------------- *
 * Filter methods
 * --------------------------------------------------------------------------------------------- */

/// 
impl Filter {


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




