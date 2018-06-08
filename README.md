# Time-Regularized Blind Deconvolution Approach for Radio Interferometry

**Description:** Matlab codes associated with the method described in 

>P.-A. Thouvenin, A. Repetti, A. Dabbech, and Y. Wiaux - 
<strong>Non-convex optimization for self-calibration of direction-dependent effects in radio interferometric imaging</strong>, <em>Proc. IEEE Sensor Array and Multichannel Signal Process. Workshop (SAM)</em>, to appear, Sheffield, United-Kingdom, 8--11 July 2018.

**Author:** P.-A. Thouvenin, p.thouvenin@hw.ac.uk

**Experiments:** to reproduce the experiments reported in the article, configure and successively run the following files:

- main_generate_synth_data.m
- main_imaging_true_dde.m
- main_calibration_imaging_reg.m
- main_calibration_imaging_no_reg.m

**Dependencies:** the present codes includes a slightly modified version of the MATLAB NUFFT algorithm available at http://web.eecs.umich.edu/~fessler/irt/fessler.tgz, described in

> J. A. Fessler and B. P. Sutton - 
<strong>Nonuniform Fast Fourier Transforms Using Min-Max Interpolation</strong>, <em>IEEE Trans. Image Process.</em>, vol. 51, n. 2, pp. 560--574, Feb. 2003.
