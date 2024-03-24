use std::str::FromStr;
use num::{Float, Complex, complex::ComplexFloat};

/// Trait common to every data types: real or complex, 4 or 8 bytes.
/// This trait contains functions to read and write one data sample from a String,
/// and functions multiply by scalar (f64)
pub trait Data: ComplexFloat {
	type Base;

	/// Multiply by a f64 factor
    fn scale(self, slope: f64) -> Self;
	/// Cast into a String
	fn to_str(self, time: f64) -> String;
	/// Read a String and return the value in the appropriate type
	fn from_string(input: String) -> Result<(f64, Self), String>
		where Self: Sized;
}

impl Data for f32 {
	type Base = f32;

    fn scale(self, slope: f64) -> Self {
        self * (slope as f32)
    }
	fn to_str(self, time: f64) -> String {
		format!("{},{}", time, self).to_string()
	}
	fn from_string(input: String) -> Result<(f64, Self), String> {
		let string_vec: Vec<&str> = input.split(",").collect();
		if string_vec.len() == 2 {
			Ok((
				f64::from_str(string_vec[0]).unwrap(),
				f32::from_str(string_vec[1]).unwrap()
			))
		} else {
			Err(format!("Error while reading string: {}", input))
		}
	}
}



impl Data for f64 {
	type Base = f64;

    fn scale(self, slope: f64) -> Self {
        self * slope
    }
	fn to_str(self, time: f64) -> String {
		format!("{},{}", time, self).to_string()
	}
	fn from_string(input: String) -> Result<(f64, Self), String> {
		let string_vec: Vec<&str> = input.split(",").collect();
		if string_vec.len() == 2 {
			Ok((
				f64::from_str(string_vec[0]).unwrap(),
				f64::from_str(string_vec[1]).unwrap()
			))
		} else {
			Err(format!("Error while reading string: {}", input))
		}
	}
}




impl Data for Complex<f32> {
	type Base = f32;

    fn scale(self, slope: f64) -> Self {
        self * Complex{re: slope as f32, im: 0f32}
    }
	fn to_str(self, time: f64) -> String {
		format!("{},{},{}", time, self.re, self.im).to_string()
	}
	fn from_string(input: String) -> Result<(f64, Self), String> {
		let string_vec: Vec<&str> = input.split(",").collect();
		if string_vec.len() == 3 {
			Ok((
				f64::from_str(string_vec[0]).unwrap(),
				Complex{
					re: f32::from_str(string_vec[1]).unwrap(), 
					im: f32::from_str(string_vec[2]).unwrap()}
			))
		} else if string_vec.len() == 2 {
			Ok((
				f64::from_str(string_vec[0]).unwrap(),
				Complex{re: f32::from_str(string_vec[1]).unwrap(), im: 0f32}
			))
		} else {
			Err(format!("Error while reading string: {}", input))
		}
	}

}




impl Data for Complex<f64> {
	type Base = f64;

    fn scale(self, slope: f64) -> Self {
        self * Complex{re: slope, im: 0f64}
    }
	fn to_str(self, time: f64) -> String {
		format!("{},{},{}", time, self.re, self.im).to_string()
	}
	fn from_string(input: String) -> Result<(f64, Self), String> {
		let string_vec: Vec<&str> = input.split(",").collect();
		if string_vec.len() == 3 {
			Ok((
				f64::from_str(string_vec[0]).unwrap(),
				Complex{
					re: f64::from_str(string_vec[1]).unwrap(),
					im: f64::from_str(string_vec[2]).unwrap()}
			))
		} else if string_vec.len() == 2 {
			Ok((
				f64::from_str(string_vec[0]).unwrap(),
				Complex{re: f64::from_str(string_vec[1]).unwrap(), im: 0f64}
			))
		} else {
			Err(format!("Error while reading string: {}", input))
		}
	}
}
/* --------------------------------------------------------------------------------------------- */

/// Trait containing specific function for float data 
pub trait FloatData: Data + Float {
	/// Return an array of data type sample
	fn linspace(start: Self, end: Self, size: usize) -> Vec<Self>;

}

impl FloatData for f32 {

	fn linspace(start: Self, end: Self, size: usize) -> Vec<Self> {
		
		// initialize output vector
		let mut output: Vec<f32> = Vec::new();
		let step: f32 = (end - start) / (size - 1) as f32;

		for i in 0..size {
			output.push(i as f32 * step);
		}
		output
	}
}
impl FloatData for f64 {

	fn linspace(start: Self, end: Self, size: usize) -> Vec<Self> {
		
		// initialize output vector
		let mut output: Vec<f64> = Vec::new();
		let step: f64 = (end - start) / (size - 1) as f64;

		for i in 0..size {
			output.push(i as f64 * step);
		}
		output
	}
}

