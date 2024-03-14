I started this project to train myself on the rust programming language. If you have questions or suggestions, or if you want to report a bug, please write an email on plagabbe@gmail.com !!

This library provides tools for signal processing and analysis. It is composed by four classes:


# TimeSeries
The `TimeSeries` object describes time series data, which consists in time indexed vector.

## Attributes
 - fs: sampling frequency
 - t0: time of the first sample
 - data: data vector

## Methods

### Constructor methods
 - white_noise: generates a white noise signal with a given amplitude
 - constant: generates a constant signal at a given value
 - wave: generates a sinusoidal signal withe the given frequency, amplitude and phase at t0
 - from_vector: transform a data vector into a `TimeSeries` object
### Methods for spectral analysis.
These methods create a `FrequencySeries` object.
 - csd: cross spectral density
 - psd: power spectral density
 - asd: amplitude spectral density
 - tf: transfer function
 - cohe: coherence
### Getter methods
 - get_fs: get sampling frequency
 - get_t0: get start time
 - get_size: get data vector size
 - get_data: get data vector
### Math methods
 - abs: computes modulus of the data
 - sqrt: computes the square root of the data
 - inv: computes the inverse of the data
### Operator overloading
 - (+) add
 - (-) subtract
 - (*) multiply
 - (/) divide
 - ([]) index
### Others
 - apply_filter: apply the `Filter` object to the time series
 - print: print the chosen elements (debug function)
 - write: wrte time and data vectors into an ascii file


# FrequencySeries
The `FrequencySeries` object describes the frequency indexed data. The vector contains complex data type.

## Attributes
 - f_max: sampling frequency
 - data: data vector

## Methods

### Constructor methods
 - from_vector: transform a data vector into a `FrequencySeries` object.
### Getter methods
 - get_f_max: get maximum frequency
 - get_size: get data vector size
 - get_data: get data vector
### Math methods
 - re: computes real part of the data
 - im: computes imaginary part of the data
 - abs: computes modulus of the data
 - abs2: computes the square modulus of the data
 - conj: computes the complex conjugates of the data
 - sqrt: computes the square root of the data
 - inv: computes the inverse of the data
### Operator overloading
 - (+) add
 - (-) subtract
 - (*) multiply
 - (/) divide
 - ([]) index
### Others
 - print: print the chosen elements (debug function)
 - write: wrte time and data vectors into an ascii file


# Window

## Attributes
 - overlap: number of samples two successives windows overlap over
 - vector: window vector

## Methods

### Constructors
 - rectangle
 - hann
### Other
 - nb_fft: computes the number of windows
 - get_size: return window size
 - get_windowed_data: compute the window data
 - get_norm_factor; computes the integral of the squared window


# Filter

## Attributes
 - gain: gain of the filter
 - poles: pole frequencies list
 - zeros: zero frequencies list
 - fs: sampling frequency

## Methods

### filter application methods
 - bilinear_transform
 - adapt_frequencies
 - polezero_to_coef
### Generic filter constructor
 - butterworth
 - chebyshev_type1
 - chebyshev_type2
### Custom filter contructor
 - init_filter
 - gain_factor
 - add_pole_1
 - add_pole_2
 - add_zero_1
 - add_zero_2
 - add_integrator
 - add_derivator
### getter function
 - get_poles
 - get_zeros
 - get_gain
 - get_fs
### Other
 - frequency_response

