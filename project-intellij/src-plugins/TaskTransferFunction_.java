// ImageJ plugin to calculate the task transfer function of a circular object
// present in an image.
//
// The routine will run on a user-selected oval region of interest (ROI). If no
// ROI is selected then the user is returned out of the routine.
//
//package ij.plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.measure.Calibration;
import ij.process.ImageStatistics;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.differentiation.*;
import org.apache.commons.math3.analysis.interpolation.*;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.MathArrays;
import ij.process.FHT;
import java.awt.*;

public class TaskTransferFunction_ implements PlugInFilter {
    private ImagePlus imp;
	private double obj_dia_mm = 25.0;
    private double pixel_reduction_factor = 5.0;
    private double loess_bandwidth = 0.02;
    private int loess_robustness = 10;
    private Roi main_roi;
    private Boolean find_com = Boolean.TRUE;

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
			IJ.log("Initial CoM: " + com[0] + ", " + com[1]);
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
			IJ.log("Adjusted CoM: " + com[0] + ", " + com[1]);

			// Recentre the annualar roi on the new CoM and recalculate CoM  for a final time
			annular_roi.setLocation(cal.getRawX(com[0]) - obj_dia_pix, cal.getRawY(com[1]) - obj_dia_pix);
			imp.setRoi(annular_roi);
			annular_stats = imp.getAllStatistics();
			bgd_mean = annular_stats.mean;
			imp.setRoi(outer_roi);
			com = centreOfMass(ip, outer_roi, cal, -bgd_mean);
		}

		IJ.log("Final CoM: " + com[0] + ", " + com[1]);

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

		IJ.log("CNR is: " + cnr_results.cnr);
		IJ.log("Contrast is: " + cnr_results.contrast);
		IJ.log("Noise is: " + cnr_results.noise);

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
		double[] unique_pos_vals = new double[unique_pos.length];

		// If the length of unique_pos is less than raw_esf_pos then there must
		// be multiple data points at the same distance from the CoM.
		if (unique_pos.length < raw_esf_pos.length) {
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
		plot_esf.show();
		String esf_labels = "Raw ESF\t";
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Interpolate the raw ESF to a regularly-spaced array.
		// Work out the distance scale increment.
		double rebinned_sample_inc = pixel_size_in_mm / pixel_reduction_factor;

		// Define the starting position, with an offset of one sample increment
		// and work out the number of samples.
		double start_pos = StatUtils.min(raw_esf_pos) + rebinned_sample_inc;
		double num_samples = Math.floor( (StatUtils.max(raw_esf_pos) - start_pos) / rebinned_sample_inc );

		// Populate an array containing the regularly-spaced locations for the ESF.
		double current_pos = start_pos;
		double[] esf_rebinned_pos = new double[(int) num_samples];
		for (i=0; i<esf_rebinned_pos.length; i++) {
			esf_rebinned_pos[i] = current_pos;
			current_pos += rebinned_sample_inc;
		}

		// Work out the value of the ESF at the rebinned positions using the
		// Apache Commons Math 3.6 API local regression algorithm. See:
		// http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/analysis/interpolation/LoessInterpolator.html
		// and https://en.wikipedia.org/wiki/Local_regression
		UnivariateInterpolator loess_interpolator = new LoessInterpolator(loess_bandwidth, loess_robustness);
		double[] smoothed_esf = new LoessInterpolator(loess_bandwidth, loess_robustness).smooth(raw_esf_pos, raw_esf_val);
		UnivariateFunction loess_function = loess_interpolator.interpolate(raw_esf_pos, smoothed_esf);
		double[] esf_rebinned_val = new double[esf_rebinned_pos.length];
		for (i=0; i<esf_rebinned_val.length; i++) {
			esf_rebinned_val[i] = loess_function.value(esf_rebinned_pos[i]);
		}

		// Add the rebinned ESF to the existing ESF plot as a red line.
		plot_esf.setColor("Red");
		plot_esf.add("line", esf_rebinned_pos, esf_rebinned_val);
		esf_labels += "Local regression ESF\t";
		plot_esf.addLegend(esf_labels);
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
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
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Differentiate the ESF to obtain a line spread function (LSF) using
		// a Apache Commons Math 3.6 API differentiator. See:
		// http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/analysis/differentiation/FiniteDifferencesDifferentiator.html
		// I have configured the FiniteDifferencesDifferentiator to use two
		// points with a step size of rebinned_sample_inc to provide the
		// gradient between each pair of points in the rebinned ESF.
		FiniteDifferencesDifferentiator differentiator = new FiniteDifferencesDifferentiator(2, rebinned_sample_inc);
		UnivariateDifferentiableFunction completeF = differentiator.differentiate(loess_function);
		double[] lsf_val = new double[esf_rebinned_pos.length];
		for (i=0; i<esf_rebinned_pos.length; i++) {
			DerivativeStructure xDS = new DerivativeStructure(1, 1, 0, esf_rebinned_pos[i]);
			DerivativeStructure yDS = completeF.value(xDS);
			lsf_val[i] = yDS.getPartialDerivative(1);
		}

		// Plot the LSF.
		Plot plot_lsf = new Plot("Line spread function", "Distance (mm)", "Value");
		plot_lsf.setColor("Red");
		plot_lsf.add("line", esf_rebinned_pos, lsf_val);
		plot_lsf.show();
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Fourier transform the LSF to obtain the MTF (TTF) using the FHT
		// routine within ImageJ. See:
		// https://imagej.nih.gov/ij/developer/api/ij/process/FHT.html#fourier1D-float:A-int-

		// The FHT routine requires the values to be of type "float".
		float[] lsf_val_float = new float[lsf_val.length];
		for (i=0; i<lsf_val.length; i++) {
			lsf_val_float[i] = (float) lsf_val[i];
		}

		// Calculate the TTF of the (float) LSF values; apply a Hann window
		// before the FFT to match ImaQuest software - see paper by Chen et al:
		// http://dx.doi.org/10.1118/1.4881519]
		// Also see http://download.ni.com/evaluation/pxi/Understanding%20FFTs%20and%20Windowing.pdf
		FHT fht = new FHT();
		float[] ttf_val = fht.fourier1D(lsf_val_float, FHT.HANN);

		// Work out the Nyquist frequency and how many TTF elements to plot in
		// order to get to 2 x Nyquist.
		double nyquist_freq = 1.0 / (2.0 * pixel_size_in_mm);
		int ttf_elements_to_plot = (int)Math.floor((2.0 * nyquist_freq) * (2.0 * ttf_val.length * rebinned_sample_inc));

		// Calculate the normalised TTF (nTTF).
		// Using double data type so they can be used directly with Plot.
		double[] nttf_val = new double[ttf_elements_to_plot];
		double[] freq_scale = new double[ttf_elements_to_plot];
		for (i=0; i<ttf_elements_to_plot; i++) {
			// ttf_val[0] is the DC component so don't use it to normalise,
			// use ttf_val[1] instead.
			nttf_val[i] = (double) (ttf_val[i] / ttf_val[1]);

			// The frequencies are calculated according to the documentation
			// under the "Returns" section of this page:
			// https://imagej.nih.gov/ij/developer/api/ij/process/FHT.html#fourier1D-float:A-int-
			freq_scale[i] = (double) i / (2.0 * ttf_val.length * rebinned_sample_inc);
		}

		// Plot the nTTF
		Plot plot_ttf = new Plot("Normalised TTF", "Frequency (per mm)", "nTTF");
		plot_ttf.setColor("Red");
		plot_ttf.add("line", freq_scale, nttf_val);
		plot_ttf.show();
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


	public class CNR_result {
		public double cnr;
		public double contrast;
		public double noise;
	}
}
