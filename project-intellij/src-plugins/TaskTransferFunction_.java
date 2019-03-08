// ImageJ plugin to calculate the task transfer function of a circular object
// present in an image.
//
// The routine will run on a user-selected oval region of interest (ROI). If no
// ROI is selected then the routine exits.
//
//package ij.plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.measure.CurveFitter;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.measure.Calibration;
import ij.process.ImageStatistics;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.differentiation.*;
import org.apache.commons.math3.analysis.interpolation.*;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.MathArrays;
import ij.process.FHT;
import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class TaskTransferFunction_ implements PlugInFilter {
    private ImagePlus imp;
	private double obj_dia_mm = 25.0;
    private double pixel_reduction_factor = 5.0;
    private double loess_bandwidth = 0.02;
    private int loess_robustness = 10;
    private Roi main_roi;
    private Boolean find_com = Boolean.TRUE;

	private DecimalFormat two_dp = new DecimalFormat("0.00");
	private DecimalFormat three_dp = new DecimalFormat("0.000");

	public int setup(String arg, ImagePlus imp) {
    	this.imp = imp;

    	if(this.imp == null) {
			IJ.error("Please open an image");
			return(0);
		}

    	// Remove any overlays that are present
		imp.setOverlay(null);

		// If there is no ROI set then exit
		main_roi = imp.getRoi();
		if(main_roi == null) {
			IJ.error("Please draw an ROI around an object");
			return(0);
		}

		if ((arg != null) && (arg.length() > 0)) {
			pixel_reduction_factor = Float.parseFloat(arg);
		}
		else {
    		 GenericDialog gd = new GenericDialog("Processing options");
    		 gd.addNumericField("Object diameter (mm)", obj_dia_mm, 0);

    		 gd.addNumericField("Pixel reduction factor", pixel_reduction_factor, 0);
    		 gd.addNumericField("LOESS bandwidth", loess_bandwidth, 2);
    		 gd.addNumericField("LOESS robustness", loess_robustness, 0);

    		 gd.addCheckbox("Find centre of mass", find_com);

    		 gd.showDialog();
    		 if(gd.wasCanceled()) {
    			 IJ.error("Plugin cancelled");
    			 return(0);
    		 }

    		 obj_dia_mm = gd.getNextNumber();
    		 pixel_reduction_factor = gd.getNextNumber();
    		 loess_bandwidth = gd.getNextNumber();
    		 loess_robustness = (int) gd.getNextNumber();
    		 find_com = gd.getNextBoolean();
    	}
    	return DOES_ALL;
    }

    public void run(ImageProcessor ip) {
		//---------------------------------------------------------------------
    	// Obtain the pixel size of the image. This will return 1.0 if ImageJ
		// isn't able to read the pixel size from the DICOM tags.
		Calibration cal = imp.getCalibration();
		double pixel_size_in_mm = cal.pixelWidth;
		double obj_dia_pix = obj_dia_mm/pixel_size_in_mm;
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Determine the centre of mass (Com). Initialise the CoM to the centre
		// of the main_roi.
		double[] com = new double[2];
		com[0] = cal.getX(main_roi.getXBase() + (main_roi.getFloatWidth()/2.0));
		com[1] = cal.getY(main_roi.getYBase() + (main_roi.getFloatHeight()/2.0));

		double inner_dia = obj_dia_pix + (10.0 / pixel_size_in_mm);
		Roi inner_roi = new OvalRoi(cal.getRawX(com[0]) - inner_dia / 2.0, cal.getRawY(com[1]) - inner_dia / 2.0, inner_dia, inner_dia);
		Roi outer_roi = new OvalRoi(cal.getRawX(com[0]) - obj_dia_pix, cal.getRawY(com[1]) - obj_dia_pix, obj_dia_pix * 2.0, obj_dia_pix * 2.0);
		Roi annular_roi = new ShapeRoi(outer_roi).xor(new ShapeRoi(inner_roi));

		if (find_com) {
			//---------------------------------------------------------------------
			// Obtain the centre of mass of the ROI set by the user.
			imp.setRoi(main_roi);
			com = centreOfMass(ip, main_roi, cal, 0.0);
			IJ.log("Initial CoM: " + three_dp.format(com[0]) + ", " + three_dp.format(com[1]));
			//---------------------------------------------------------------------


			//---------------------------------------------------------------------
			// Recalculate the CoM after subtracting the background using an
			// annular region of interest centred on the centre of mass that has an
			// internal diameter of object diameter + 10 mm and an external
			// diameter of 2 x object diameter.
			inner_roi = new OvalRoi(cal.getRawX(com[0]) - inner_dia / 2.0, cal.getRawY(com[1]) - inner_dia / 2.0, inner_dia, inner_dia);
			outer_roi = new OvalRoi(cal.getRawX(com[0]) - obj_dia_pix, cal.getRawY(com[1]) - obj_dia_pix, obj_dia_pix * 2.0, obj_dia_pix * 2.0);
			annular_roi = new ShapeRoi(outer_roi).xor(new ShapeRoi(inner_roi));


			// Obtain the mean pixel value of the annular ROI
			imp.setRoi(annular_roi);
			ImageStatistics annular_stats = imp.getAllStatistics();
			double bgd_mean = annular_stats.mean;

			// Recalculate the CoM with background subtracted
			imp.setRoi(outer_roi);
			com = centreOfMass(ip, outer_roi, cal, -bgd_mean);
			IJ.log("Adjusted CoM: " + three_dp.format(com[0]) + ", " + three_dp.format(com[1]));

			// Recentre the annualar roi on the new CoM and recalculate CoM  for a final time
			annular_roi.setLocation(cal.getRawX(com[0]) - obj_dia_pix, cal.getRawY(com[1]) - obj_dia_pix);
			imp.setRoi(annular_roi);
			annular_stats = imp.getAllStatistics();
			bgd_mean = annular_stats.mean;
			imp.setRoi(outer_roi);
			com = centreOfMass(ip, outer_roi, cal, -bgd_mean);
		}

		IJ.log("Final CoM: " + three_dp.format(com[0]) + ", " + three_dp.format(com[1]));

		// Add a marker to the image showing the location of the CoM
		Roi com_marker = new PointRoi(cal.getRawX(com[0]), cal.getRawY(com[1]));
		Overlay com_overlay = new Overlay(com_marker);
		com_overlay.setStrokeColor(Color.red);
		imp.setOverlay(com_overlay);

		// Recentre the outer_roi using the centre of mass.
		outer_roi.setLocation(cal.getRawX(com[0])-obj_dia_pix, cal.getRawY(com[1])-obj_dia_pix);

		// Create an overlay to display the annular ROI and CoM on the image.
		Overlay overlay = new Overlay();
		imp.setOverlay(overlay);
		overlay.add(annular_roi);
		overlay.add(com_marker);
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Calculate the CNR of the object
		double obj_roi_dia = obj_dia_pix - (3.0 / pixel_size_in_mm);
		Roi obj_roi = new OvalRoi(cal.getRawX(com[0])-obj_roi_dia/2.0, cal.getRawY(com[1])-obj_roi_dia/2.0, obj_roi_dia, obj_roi_dia);
		CNR_result cnr_results = contrastToNoiseRatio(imp, obj_roi, annular_roi);

		IJ.log("CNR is: " + three_dp.format(cnr_results.cnr));
		IJ.log("Contrast is: " + three_dp.format(cnr_results.contrast));
		IJ.log("Noise is: " + three_dp.format(cnr_results.noise));

		overlay.add(obj_roi);
		overlay.setStrokeColor(Color.green);
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Obtain the raw edge spread function (ESF).
		int i=0;
		int num_points = outer_roi.getContainedPoints().length;
		double[] raw_esf_pos = new double[num_points];
		double[] raw_esf_val = new double[num_points];
		for (Point p : outer_roi.getContainedPoints()) {
		    raw_esf_pos[i] = calculateDistanceBetweenPoints(com[0], com[1], cal.getX(p.x), cal.getY(p.y));
			raw_esf_val[i] = ip.getPixelValue(p.x, p.y); // getPixelValue includes calibration
            i++;
		}

		// There's a chance that two or more elements are at exactly the same
		// position. This violates the monotonic requirement of the
		// interpolator used later on. So find duplicate position entries and
		// replace with the mean value at that position.
		// http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/util/MathArrays.html#unique(double[])
		double[] unique_pos = MathArrays.unique(raw_esf_pos);

		// If the length of unique_pos is less than raw_esf_pos then there must
		// be multiple data points at the same distance from the CoM.
		if (unique_pos.length < raw_esf_pos.length) {
			double[] unique_pos_vals = new double[unique_pos.length];

			for (i = 0; i < unique_pos.length; i++) {
				int count = 0;
				double sum = 0.0;
				for (int j = 0; j < raw_esf_pos.length; j++) {
					if (raw_esf_pos[j] == unique_pos[i]) {
						count++;
						sum += raw_esf_val[j];
					}
				}
				unique_pos_vals[i] = sum / count;
			}

			// Overwrite the initial position and value data with the new set
			// of data at unique positions.
			raw_esf_pos = unique_pos;
			raw_esf_val = unique_pos_vals;
		}

		// Sort the raw ESF arrays into ascending position order using Apache
		// Commons Math 3.6 API. See:
		// http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/util/MathArrays.html#sortInPlace(double[],%20double[]...)
		MathArrays.sortInPlace(raw_esf_pos, raw_esf_val);

		// Plot the raw esf
		Plot plot_esf = new Plot("Edge spread function", "Distance (mm)", "Pixel value");
		plot_esf.add("dot", raw_esf_pos, raw_esf_val);
		String esf_labels = "Raw ESF\t";
		plot_esf.show();
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Interpolate the raw ESF to a regularly-spaced array with sub-pixel
		// sampling using a local regression algorithm.  See:
		// http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/analysis/interpolation/LoessInterpolator.html
		// and https://en.wikipedia.org/wiki/Local_regression

		// Work out the distance scale increment.
		double rebinned_sample_inc = pixel_size_in_mm / pixel_reduction_factor;

		// Define the starting position of the resampled ESF, with an offset of
		// one sample increment, and work out the number of samples.
		double start_pos = StatUtils.min(raw_esf_pos) + rebinned_sample_inc;
		double num_samples = Math.floor( (StatUtils.max(raw_esf_pos) - start_pos) / rebinned_sample_inc );

		// Populate an array containing the regularly-spaced locations for the ESF.
		double current_pos = start_pos;
		double[] esf_rebinned_pos = new double[(int) num_samples];
		for (i=0; i<esf_rebinned_pos.length; i++) {
			esf_rebinned_pos[i] = current_pos;
			current_pos += rebinned_sample_inc;
		}

		//.....................................................................
		// Calculate resampled ESF values using local regression with the
		// Apache Commons Math 3.6 API. See:
		// http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/analysis/interpolation/LoessInterpolator.html
		// and https://en.wikipedia.org/wiki/Local_regression
		UnivariateInterpolator loess_interpolator = new LoessInterpolator(loess_bandwidth, loess_robustness);
		double[] smoothed_esf = new LoessInterpolator(loess_bandwidth, loess_robustness).smooth(raw_esf_pos, raw_esf_val);
		UnivariateFunction loess_function = loess_interpolator.interpolate(raw_esf_pos, smoothed_esf);
		double[] esf_rebinned_val = new double[esf_rebinned_pos.length];
		for (i=0; i<esf_rebinned_val.length; i++) {
			esf_rebinned_val[i] = loess_function.value(esf_rebinned_pos[i]);
		}

		// Add the ESF to the existing ESF plot as a red line.
		plot_esf.setColor("Red");
		plot_esf.add("line", esf_rebinned_pos, esf_rebinned_val);
		esf_labels += "Local regression ESF\t";
		//.....................................................................


		//.....................................................................
		// Create a monotonic ESF from the local regression ESF using brute
		// force.
		double[] monotonic_esf = new double[esf_rebinned_val.length];
		// Check to see if ESF starts high and goes low
		if (esf_rebinned_val[0] > esf_rebinned_val[esf_rebinned_val.length-1]) {
			// Need to reverse the values array
			MathArrays.sortInPlace(esf_rebinned_pos, MathArrays.OrderDirection.DECREASING, esf_rebinned_val);
		}
		monotonic_esf[0] = esf_rebinned_val[0];
		for (i=1; i<esf_rebinned_val.length; i++) {
			if (esf_rebinned_val[i] >= monotonic_esf[i-1]) {
				monotonic_esf[i] = esf_rebinned_val[i];
			}
			else {
				monotonic_esf[i] = monotonic_esf[i-1];
			}
		}

		// Sort the esf_rebinned_pos and associated arrays back to ascending order
		MathArrays.sortInPlace(esf_rebinned_pos, esf_rebinned_val, monotonic_esf);

		// Plot the monotonic ESF
		plot_esf.setColor("Blue");
		plot_esf.add("line", esf_rebinned_pos, monotonic_esf);
		esf_labels += "Monotonic local regression ESF\t";
		//.....................................................................


		//.....................................................................
		// Calculate resampled ESF values using linear interpolation
		UnivariateInterpolator linear_interpolator = new LinearInterpolator();
		UnivariateFunction linear_function = linear_interpolator.interpolate(raw_esf_pos, smoothed_esf);
		double[] esf_rebinned_val_lin = new double[esf_rebinned_pos.length];
		for (i=0; i<esf_rebinned_val_lin.length; i++) {
			esf_rebinned_val_lin[i] = linear_function.value(esf_rebinned_pos[i]);
		}

		// Add the rebinned linear ESF to the existing ESF plot as a green line.
		plot_esf.setColor("Green");
		plot_esf.add("line", esf_rebinned_pos, esf_rebinned_val_lin);
		esf_labels += "Linear regression ESF\t";
		plot_esf.addLegend(esf_labels);
		//.....................................................................


/*		//---------------------------------------------------------------------
		// TEST: Rebin the raw ESF into a regularly-spaced array using the
		// method described by Samei et at here: http://doi.org/10.1118/1.598165
		start_pos = StatUtils.min(raw_esf_pos);
		num_samples = Math.floor( (StatUtils.max(raw_esf_pos) - start_pos) / rebinned_sample_inc );
		double[] esf_binned_pos = new double[(int) num_samples];
		double[] esf_binned_val = new double[esf_binned_pos.length];
		current_pos = start_pos;
		for (i=0; i<esf_binned_pos.length; i++) {
			int count = 0;
			double sum = 0.0;

			esf_binned_pos[i] = current_pos;

			for (int j=0; j<raw_esf_pos.length; j++) {
				if ( (raw_esf_pos[j] >= esf_binned_pos[i] - rebinned_sample_inc/2.0) && (raw_esf_pos[j] <= esf_binned_pos[i] + rebinned_sample_inc/2.0) ) {
					count ++;
					sum += raw_esf_val[j];
				}
			}

			if (count > 0) {
				esf_binned_val[i] = sum / count;
			}
			else esf_binned_val[i] = esf_binned_val[i-1];

			current_pos += rebinned_sample_inc;
		}

		//plot_esf.setColor("Blue");
		//plot_esf.add("line", esf_binned_pos, esf_binned_val);
		//esf_labels += "Rebinned ESF";
		//plot_esf.addLegend(esf_labels);
		// TEST end
		//---------------------------------------------------------------------*/

		// End of creating ESFs
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Differentiate the ESFs to obtain line spread functions (LSFs)

		//.....................................................................
		// First using an Apache Commons Math 3.6 API differentiator. See:
		// http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/analysis/differentiation/FiniteDifferencesDifferentiator.html
		// I have configured the FiniteDifferencesDifferentiator to use two
		// points with a step size of rebinned_sample_inc to provide the
		// gradient between each pair of points in the rebinned ESF.
		FiniteDifferencesDifferentiator differentiator = new FiniteDifferencesDifferentiator(2, rebinned_sample_inc);
		UnivariateDifferentiableFunction completeF = differentiator.differentiate(loess_function);
		double[] lsf_val = differentiate(esf_rebinned_pos, completeF);

		// Plot the LSF.
		Plot plot_lsf = new Plot("Line spread function", "Distance (mm)", "Value");
		plot_lsf.setColor("Red");
		plot_lsf.add("line", esf_rebinned_pos, lsf_val);
		String lsf_labels = "LSF from local regression ESF\t";
		plot_lsf.show();
		//.....................................................................


		//.....................................................................
		// Differentiate the monotonic ESF to produce a LSF
		//UnivariateInterpolator linear_interpolator = new LinearInterpolator();
		linear_function = linear_interpolator.interpolate(esf_rebinned_pos, monotonic_esf);
		completeF = differentiator.differentiate(linear_function);
		double[] lsf_val_monotonic = differentiate(esf_rebinned_pos, completeF);

		// Plot the LSF calculated from the monotonic ESF.
		plot_lsf.setColor("Blue");
		plot_lsf.add("line", esf_rebinned_pos, lsf_val_monotonic);
		lsf_labels += "LSF from monotonic ESF\t";
		//.....................................................................


		//.....................................................................
		// Differentiate the linear ESF to produce a LSF
		linear_function = linear_interpolator.interpolate(esf_rebinned_pos, esf_rebinned_val_lin);
		completeF = differentiator.differentiate(linear_function);
		double[] lsf_val_linear = differentiate(esf_rebinned_pos, completeF);

		// Plot the LSF calculated from the monotonic ESF.
		plot_lsf.setColor("Green");
		plot_lsf.add("line", esf_rebinned_pos, lsf_val_linear);
		lsf_labels += "LSF from linear ESF\t";
		//.....................................................................


		//.....................................................................
		// Carry out a Gaussian fit to the linear LSF result
		CurveFitter gaussian_fitter = new CurveFitter(esf_rebinned_pos, lsf_val_linear);
		gaussian_fitter.doFit(CurveFitter.GAUSSIAN);
		double[] p = gaussian_fitter.getParams();
		double[] lsf_val_gaussian = new double[esf_rebinned_pos.length];
		for (i=0; i<esf_rebinned_pos.length; i++) {
			double x = esf_rebinned_pos[i];
			lsf_val_gaussian[i] = p[0] + (p[1] - p[0]) * Math.exp(-(x - p[2]) * (x - p[2]) / (2.0 * p[3] * p[3]));
		}

		// Plot the LSF calculated from the Gaussian fit of the linear ESF.
		plot_lsf.setColor("Orange");
		plot_lsf.add("line", esf_rebinned_pos, lsf_val_gaussian);
		lsf_labels += "LSF from Gaussian fit of linear ESF\t";
		plot_lsf.addLegend(lsf_labels);
		//.....................................................................

		// End of LSF calculation
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Fourier transform the LSFs to obtain TTFs
		TTF_result ttf_data = TTF(lsf_val, rebinned_sample_inc, pixel_size_in_mm);
		TTF_result ttf_data_linear = TTF(lsf_val_linear, rebinned_sample_inc, pixel_size_in_mm);
		TTF_result ttf_data_gaussian = TTF(lsf_val_gaussian, rebinned_sample_inc, pixel_size_in_mm);
		TTF_result ttf_data_monotonic = TTF(lsf_val_monotonic, rebinned_sample_inc, pixel_size_in_mm);

		// Plot the nTTFs
		Plot plot_ttf = new Plot("Normalised TTF", "Frequency (per mm)", "nTTF");

		plot_ttf.setColor("Red");
		plot_ttf.add("line", ttf_data.freq, ttf_data.val);
		String ttf_labels = "TTF from local regression ESF\t";

		plot_ttf.setColor("Green");
		plot_ttf.add("line", ttf_data_linear.freq, ttf_data_linear.val);
		ttf_labels += "TTF from linear ESF\t";

		plot_ttf.setColor("Orange");
		plot_ttf.add("line", ttf_data_gaussian.freq, ttf_data_gaussian.val);
		ttf_labels += "TTF from Gaussian LSF fit\t";

		plot_ttf.setColor("Blue");
		plot_ttf.add("line", ttf_data_monotonic.freq, ttf_data_monotonic.val);
		ttf_labels += "TTF from monotonic ESF\t";

		plot_ttf.show();
		plot_ttf.addLegend(ttf_labels);
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Work out TTF 50 and TTF 10 for each TTF
		double ttf_50 = freqAtSpecificTTF(ttf_data, 0.5);
		double ttf_10 = freqAtSpecificTTF(ttf_data, 0.1);

		double ttf_50_linear = freqAtSpecificTTF(ttf_data_linear, 0.5);
		double ttf_10_linear = freqAtSpecificTTF(ttf_data_linear, 0.1);

		double ttf_50_gaussian = freqAtSpecificTTF(ttf_data_gaussian, 0.5);
		double ttf_10_gaussian = freqAtSpecificTTF(ttf_data_gaussian, 0.1);

		double ttf_50_mono = freqAtSpecificTTF(ttf_data_monotonic, 0.5);
		double ttf_10_mono = freqAtSpecificTTF(ttf_data_monotonic, 0.1);

		IJ.log("MTF 50 and 10 are: " + three_dp.format(ttf_50) + ", " + three_dp.format(ttf_10));
		IJ.log("MTF linear 50 and 10 are: " + three_dp.format(ttf_50_linear) + ", " + three_dp.format(ttf_10_linear));
		IJ.log("MTF Gaussian 50 and 10 are: " + three_dp.format(ttf_50_gaussian) + ", " + three_dp.format(ttf_10_gaussian));
		IJ.log("MTF mono 50 and 10 are: " + three_dp.format(ttf_50_mono) + ", " + three_dp.format(ttf_10_mono));
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Resample the TTF data to 0.05 lp/mm increments
		TTF_result resampled_ttf = resampleTTF(ttf_data, 0.02);
		TTF_result resampled_ttf_linear = resampleTTF(ttf_data_linear, 0.02);
		TTF_result resampled_ttf_gaussian = resampleTTF(ttf_data_gaussian, 0.02);
		TTF_result resampled_ttf_mono = resampleTTF(ttf_data_monotonic, 0.02);
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Write the resampled TTF data to the ImageJ log
		logTTF(resampled_ttf, "TTF (local regression ESF)");
		logTTF(resampled_ttf_linear, "TTF (linear ESF)");
		logTTF(resampled_ttf_gaussian, "TTF (Gaussian fit to LSF)");
		logTTF(resampled_ttf_mono, "TTF (monotonic ESF)");
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Remove the current ROI
		imp.deleteRoi();
		//---------------------------------------------------------------------
	}


	private double calculateDistanceBetweenPoints(double x1, double y1, double x2,double y2) {
		// Calculate the distance between two points in cartesian space.
        return Math.sqrt( (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }


	private double[] centreOfMass(ImageProcessor ip, Roi roi, Calibration cal, double pv_offset) {
		// Calculate the centre of mass of pixels in a Roi positioned on an
		// ImageProcessor with a specific calibration. The user can apply a
		// pixel value offset to subtract background.
		double x_wt=0;
		double y_wt=0;
		double total_wt=0;
		for (Point p : roi.getContainedPoints()) {
			double pixel_val = ip.getPixelValue(p.x, p.y);
			pixel_val += pv_offset;
			x_wt += cal.getX(p.x) * pixel_val;
			y_wt += cal.getY(p.y) * pixel_val;
			total_wt += pixel_val;
		}
		double x_com = x_wt / total_wt;
		double y_com = y_wt / total_wt;
		return new double[] {x_com, y_com};
	}


	private CNR_result contrastToNoiseRatio(ImagePlus imp, Roi obj_roi, Roi bgd_roi) {
		// Calculate the contrast to noise ratio between the two regions of interest.
		// CNR calculated as (obj_mean - bgd_mean) / bgd_stdDev as show in equation
		// 1 of Christianson et al, https://doi.org/10.1148/radiol.15132091
		imp.setRoi(obj_roi);
		ImageStatistics obj_stats = imp.getAllStatistics();

		imp.setRoi(bgd_roi);
		ImageStatistics bgd_stats = imp.getAllStatistics();

		CNR_result result = new CNR_result();
		result.contrast = obj_stats.mean - bgd_stats.mean;
		result.noise = bgd_stats.stdDev;
		result.cnr = Math.abs(result.contrast) / result.noise;

		return result;
	}


	private TTF_result TTF(double[] lsf, double lsf_inc_mm, double px_size_mm) {
		// This uses the FHT routine within ImageJ. See:
		// https://imagej.nih.gov/ij/developer/api/ij/process/FHT.html#fourier1D-float:A-int-

		// The FHT routine requires the values to be of type "float".
		float[] lsf_float = new float[lsf.length];
		for (int i=0; i<lsf.length; i++) {
			lsf_float[i] = (float) lsf[i];
		}

		// Calculate the TTF of the (float) LSF values; apply a Hann window
		// before the FFT to match ImaQuest software - see paper by Chen et al:
		// http://dx.doi.org/10.1118/1.4881519]
		// Also see http://download.ni.com/evaluation/pxi/Understanding%20FFTs%20and%20Windowing.pdf
		FHT fht = new FHT();
		float[] ttf_val = fht.fourier1D(lsf_float, FHT.HANN);

		// Work out the Nyquist frequency and how many TTF elements to plot in
		// order to get to 2 x Nyquist.
		double nyquist_freq = 1.0 / (2.0 * px_size_mm);
		int ttf_elements_to_plot = (int)Math.floor((2.0 * nyquist_freq) * (2.0 * ttf_val.length * lsf_inc_mm));

		// Calculate the normalised TTF (nTTF). Using double data type so they
		// can be used directly with Plot.
		double[] nttf_val = new double[ttf_elements_to_plot-1];
		double[] freq_scale = new double[ttf_elements_to_plot-1];
		// Start at element 1 of ttf_val as we don't want to use the DC
		// component value that is stored in element 0.
		for (int i=1; i<ttf_elements_to_plot; i++) {
			// ttf_val[0] is the DC component so don't use it to normalise,
			// use ttf_val[1] instead.
			nttf_val[i-1] = (double) (ttf_val[i] / ttf_val[1]);

			// The frequencies are calculated according to the documentation
			// under the "Returns" section of this page:
			// https://imagej.nih.gov/ij/developer/api/ij/process/FHT.html#fourier1D-float:A-int-
			freq_scale[i-1] = (double) i / (2.0 * ttf_val.length * lsf_inc_mm);
		}

		TTF_result results = new TTF_result();
		results.freq = freq_scale;
		results.val = nttf_val;

		return results;
	}


	private double[] differentiate(double[] locations, UnivariateDifferentiableFunction diff_func) {
		double[] deriv_val = new double[locations.length];
		for (int i=0; i<locations.length; i++) {
			try {
				DerivativeStructure xDS = new DerivativeStructure(1, 1, 0, locations[i]);
				DerivativeStructure yDS = diff_func.value(xDS);
				deriv_val[i] = yDS.getPartialDerivative(1);
			}
			catch (OutOfRangeException e) {
				// This is thrown for the first and last locations if a linear
				// interpolator is used because the interpolator doesn't have a
				// data point at either side of the required position to use
				// for its calculation.
			}
		}
		return deriv_val;
	}


	private TTF_result resampleTTF(TTF_result ttf, double freq_inc) {
		UnivariateInterpolator linear_interpolator = new LinearInterpolator();
		UnivariateFunction fn_ttf_lin = linear_interpolator.interpolate(ttf.freq, ttf.val);

		ArrayList<Double> ttf_freq = new ArrayList<Double>();
		ArrayList<Double> ttf_val = new ArrayList<Double>();

		double freq = 0.0;
		while (freq < StatUtils.max(ttf.freq)) {
			try {
				ttf_val.add(fn_ttf_lin.value(freq));
				ttf_freq.add(freq);
			}
			catch (OutOfRangeException e) {
				// If the interpolator doesn't have data on either side of the
				// required position then it throws an error.
			}
			freq += freq_inc;
		}

		double[] freq_array = new double[ttf_freq.size()];
		double[] val_array = new double[ttf_val.size()];
		for (int i=0; i<freq_array.length; i++) {
			freq_array[i] = ttf_freq.get(i);
			val_array[i] = ttf_val.get(i);
		}

		TTF_result result = new TTF_result();
		result.freq = freq_array;
		result.val = val_array;

		return result;
    }


	private void logTTF(TTF_result ttf, String ttf_label) {
		IJ.log("Freq (mm^-1)," + ttf_label);

		for (int i=0; i<ttf.freq.length; i++) {
			IJ.log("" + two_dp.format(ttf.freq[i]) + ", " + three_dp.format(ttf.val[i]));
		}
	}


	private double freqAtSpecificTTF(TTF_result ttf, double req_ttf_val) {
		UnivariateInterpolator linear_interpolator = new LinearInterpolator();
		UnivariateFunction fn_ttf_lin = linear_interpolator.interpolate(ttf.freq, ttf.val);

		double current_ttf = 100.0;
		double prev_ttf = 100.0;

		double freq_step = 0.01;

		double current_freq = ttf.freq[0];
		double prev_freq = ttf.freq[0];

		while (current_ttf > req_ttf_val) {
			prev_ttf = current_ttf;
			prev_freq = current_freq - freq_step;

			try {
				current_ttf = fn_ttf_lin.value(current_freq);
			}
			catch (OutOfRangeException e) {
				// If the interpolator doesn't have data on either side of the
				// required position then it throws an error. This may also
				// be thrown if there is a problem with the TTF data that means
				// this while loop never reaches the required ttf value.
				return 999.0;
			}
			current_freq += freq_step;
		}
		double[] ttf_vals = {prev_ttf, current_ttf};
		double[] ttf_freqs = {prev_freq, current_freq-freq_step};

		MathArrays.sortInPlace(ttf_vals, ttf_freqs);
		fn_ttf_lin = linear_interpolator.interpolate(ttf_vals, ttf_freqs);

		return fn_ttf_lin.value(req_ttf_val);
	}


	private class TTF_result {
    	private double[] freq;
    	private double[] val;
	}


	private class CNR_result {
		private double cnr;
		private double contrast;
		private double noise;
	}
}
