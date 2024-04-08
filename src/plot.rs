/* --------------------------------------------------------------------------------------------- *
 * import libraries
 * --------------------------------------------------------------------------------------------- */

use plotters::prelude::*;
use num::Complex;
use std::f64::consts::PI;

/* --------------------------------------------------------------------------------------------- *
 * define structures
 * --------------------------------------------------------------------------------------------- */

pub struct RealPlot {
	x_values: Vec<Vec<f64>>,
	y_values: Vec<Vec<f64>>,
	x_scale_log: bool,
	y_scale_log: bool,
}


pub struct ComplexPlot {
	x_values: Vec<Vec<f64>>,
	y_values: Vec<Vec<Complex<f64>>>,
	x_scale_log: bool,
	y_scale_log: bool,
}



/* --------------------------------------------------------------------------------------------- *
 * define methods
 * --------------------------------------------------------------------------------------------- */

impl RealPlot {
	/// Initialize real plot canvas
	pub fn new() -> Self {

		RealPlot{
			x_values: Vec::new(),
			y_values: Vec::new(),
			x_scale_log: false,
			y_scale_log: false,
		}
	}
	
	/// Get number of data sample
	pub fn get_size(&self) -> Vec<usize> {
		let mut output: Vec<usize> = Vec::new();
		for x in self.x_values.iter() {
			output.push(x.len());
		}
		output
	}
	/// Add one data vector to data list
	pub fn add_data_vector(&mut self, x: Vec<f64>, y: Vec<f64>) {
		assert_eq!(x.len(), y.len());
		self.x_values.push(x);
		self.y_values.push(y);
	}
	/// Set axis to log scale
	pub fn set_x_scale_to_log(&mut self, x_log: bool) {
		self.x_scale_log = x_log;
	}
	pub fn set_y_scale_to_log(&mut self, y_log: bool) {
		self.y_scale_log = y_log;
	}
	/// Draw plot
	pub fn draw(self, name: &str) {

		// define x axis
		let (mut xmin, mut xmax): (f64, f64) = (f64::MAX, f64::MIN);
		for x_vec in self.x_values.iter() {
			let (temp_min, temp_max) = (x_vec[0], x_vec[x_vec.len()-1]);
			if temp_min < xmin { xmin = temp_min; }
			if temp_max > xmax { xmax = temp_max; }
		}
		
		// define y axis
		let (mut ymin, mut ymax): (f64, f64) = (f64::MAX, f64::MIN);
		for y_vec in self.y_values.iter() {
			for y in y_vec.iter() {
				if *y < ymin { ymin = *y; }
				if *y > ymax { ymax = *y; }
			}
		}
		
		// define white canvas
        let drawing_area = SVGBackend::new(name, (1365, 768)).into_drawing_area();
        drawing_area.fill(&WHITE).unwrap();

        // define plot context
		let size: Vec<usize> = self.get_size();
		if self.x_scale_log & self.y_scale_log {
			
			let mut ctx = ChartBuilder::on(&drawing_area)
            	.margin(20)
            	.set_left_and_bottom_label_area_size(64)
				.build_cartesian_2d((xmin..xmax).log_scale(), (ymin..ymax).log_scale())
				.unwrap();
			ctx.configure_mesh().draw().unwrap();
			
			// draw data
			for i in 0..self.y_values.len() {
				ctx.draw_series(LineSeries::new( (0..size[i]).map( |index| (
					self.x_values[i][index], self.y_values[i][index]
				)), &Palette99::pick(i) )).unwrap();
			}

		} else if self.x_scale_log & !self.y_scale_log {

			let mut ctx = ChartBuilder::on(&drawing_area)
            	.margin(20)
            	.set_left_and_bottom_label_area_size(64)
				.build_cartesian_2d((xmin..xmax).log_scale(), ymin..ymax)
				.unwrap();
			ctx.configure_mesh().draw().unwrap();

			// draw data
			for i in 0..self.y_values.len() {
				ctx.draw_series(LineSeries::new( (0..size[i]).map( |index| (
					self.x_values[i][index], self.y_values[i][index]
				)), &Palette99::pick(i) )).unwrap();
			}

		} else if !self.x_scale_log & self.y_scale_log {
			let mut ctx = ChartBuilder::on(&drawing_area)
            	.margin(20)
            	.set_left_and_bottom_label_area_size(64)
				.build_cartesian_2d(xmin..xmax, (ymin..ymax).log_scale())
				.unwrap();
			ctx.configure_mesh().draw().unwrap();

			// draw data
			for i in 0..self.y_values.len() {
				ctx.draw_series(LineSeries::new( (0..size[i]).map( |index| (
					self.x_values[i][index], self.y_values[i][index]
				)), &Palette99::pick(i) )).unwrap();
			}

		} else {

			let mut ctx = ChartBuilder::on(&drawing_area)
            	.margin(20)
            	.set_left_and_bottom_label_area_size(64)
				.build_cartesian_2d(xmin..xmax, ymin..ymax)
				.unwrap();
			ctx.configure_mesh().draw().unwrap();

			// draw data
			for i in 0..self.y_values.len() {
				ctx.draw_series(LineSeries::new( (0..size[i]).map( |index| (
					self.x_values[i][index], self.y_values[i][index]
				)), &Palette99::pick(i) )).unwrap();
			}

		}
	}

}



impl ComplexPlot {
	/// Initialize real plot canvas
	pub fn new() -> Self {

		ComplexPlot{
			x_values: Vec::new(),
			y_values: Vec::new(),
			x_scale_log: false,
			y_scale_log: false,
		}
	}
	
	/// Get number of data sample
	pub fn get_size(&self) -> Vec<usize> {
		let mut output: Vec<usize> = Vec::new();
		for x in self.x_values.iter() {
			output.push(x.len());
		}
		output
	}
	/// Add one data vector to data list
	pub fn add_data_vector(&mut self, x: Vec<f64>, y: Vec<Complex<f64>>) {
		assert_eq!(x.len(), y.len());
		self.x_values.push(x);
		self.y_values.push(y);
	}
	/// Set axis to log scale
	pub fn set_x_scale_to_log(&mut self, x_log: bool) {
		self.x_scale_log = x_log;
	}
	pub fn set_y_scale_to_log(&mut self, y_log: bool) {
		self.y_scale_log = y_log;
	}
	/// Draw plot
	pub fn draw(self, name: &str) {

		// define x axis
		let (mut xmin, mut xmax): (f64, f64) = (f64::MAX, f64::MIN);
		for x_vec in self.x_values.iter() {
			let (temp_min, temp_max) = (x_vec[0], x_vec[x_vec.len()-1]);
			if temp_min < xmin { xmin = temp_min; }
			if temp_max > xmax { xmax = temp_max; }
		}
		
		// define y axis
		let (mut ymin, mut ymax): (f64, f64) = (f64::MAX, f64::MIN);
		for y_vec in self.y_values.iter() {
			for y in y_vec.iter() {
				if y.norm() < ymin { ymin = y.norm(); }
				if y.norm() > ymax { ymax = y.norm(); }
			}
		}
		
		// define white canvas
        let drawing_area = SVGBackend::new(name, (1365, 768)).into_drawing_area();
        drawing_area.fill(&WHITE).unwrap();
		
		let (top, bottom) = drawing_area.split_vertically(384);

        // define top plot context
		let size: Vec<usize> = self.get_size();
		if self.x_scale_log & self.y_scale_log {
			
			let mut ctx = ChartBuilder::on(&top)
            	.margin(20)
            	.set_left_and_bottom_label_area_size(64)
				.build_cartesian_2d((xmin..xmax).log_scale(), (ymin..ymax).log_scale())
				.unwrap();
			ctx.configure_mesh().draw().unwrap();
			
			// draw data
			for i in 0..self.y_values.len() {
				ctx.draw_series(LineSeries::new( (0..size[i]).map( |index| (
					self.x_values[i][index], self.y_values[i][index].norm()
				)), &Palette99::pick(i) )).unwrap();
			}

		} else if self.x_scale_log & !self.y_scale_log {

			let mut ctx = ChartBuilder::on(&top)
            	.margin(20)
            	.set_left_and_bottom_label_area_size(64)
				.build_cartesian_2d((xmin..xmax).log_scale(), ymin..ymax)
				.unwrap();
			ctx.configure_mesh().draw().unwrap();

			// draw data
			for i in 0..self.y_values.len() {
				ctx.draw_series(LineSeries::new( (0..size[i]).map( |index| (
					self.x_values[i][index], self.y_values[i][index].norm()
				)), &Palette99::pick(i) )).unwrap();
			}

		} else if !self.x_scale_log & self.y_scale_log {
			let mut ctx = ChartBuilder::on(&top)
            	.margin(20)
            	.set_left_and_bottom_label_area_size(64)
				.build_cartesian_2d(xmin..xmax, (ymin..ymax).log_scale())
				.unwrap();
			ctx.configure_mesh().draw().unwrap();

			// draw data
			for i in 0..self.y_values.len() {
				ctx.draw_series(LineSeries::new( (0..size[i]).map( |index| (
					self.x_values[i][index], self.y_values[i][index].norm()
				)), &Palette99::pick(i) )).unwrap();
			}

		} else {

			let mut ctx = ChartBuilder::on(&top)
            	.margin(20)
            	.set_left_and_bottom_label_area_size(64)
				.build_cartesian_2d(xmin..xmax, ymin..ymax)
				.unwrap();
			ctx.configure_mesh().draw().unwrap();

			// draw data
			for i in 0..self.y_values.len() {
				ctx.draw_series(LineSeries::new( (0..size[i]).map( |index| (
					self.x_values[i][index], self.y_values[i][index].norm()
				)), &Palette99::pick(i) )).unwrap();
			}

		}

		// fill bottom plot
		(ymin, ymax) = (-PI, PI);

		if self.x_scale_log {
			let mut ctx = ChartBuilder::on(&bottom)
            	.margin(20)
            	.set_left_and_bottom_label_area_size(64)
				.build_cartesian_2d((xmin..xmax).log_scale(), ymin..ymax)
				.unwrap();
			ctx.configure_mesh().draw().unwrap();

			// draw data
			for i in 0..self.y_values.len() {
				ctx.draw_series(LineSeries::new( (0..size[i]).map( |index| (
					self.x_values[i][index], self.y_values[i][index].arg()
				)), &Palette99::pick(i) )).unwrap();
			}
		} else {
			let mut ctx = ChartBuilder::on(&bottom)
            	.margin(20)
            	.set_left_and_bottom_label_area_size(64)
				.build_cartesian_2d(xmin..xmax, ymin..ymax)
				.unwrap();
			ctx.configure_mesh().draw().unwrap();

			// draw data
			for i in 0..self.y_values.len() {
				ctx.draw_series(LineSeries::new( (0..size[i]).map( |index| (
					self.x_values[i][index], self.y_values[i][index].arg()
				)), &Palette99::pick(i) )).unwrap();
			}
		}

	}

}








