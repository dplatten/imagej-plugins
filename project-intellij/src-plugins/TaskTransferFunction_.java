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
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
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

    public int setup(String arg, ImagePlus imp) {
    	this.imp = imp;

    	if(this.imp == null) {
			IJ.error("Please open an image");
			return(0);
		}

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

    		 gd.showDialog();
    		 if(gd.wasCanceled()) {
    			 IJ.error("Plugin cancelled");
    			 return(0);
    		 }

    		 obj_dia_mm = gd.getNextNumber();
    		 pixel_reduction_factor = gd.getNextNumber();
    		 loess_bandwidth = gd.getNextNumber();
    		 loess_robustness = (int) gd.getNextNumber();
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
    	// Obtain the centre of mass of the ROI set by the user.
		imp.setRoi(main_roi);
		double[] com = centreOfMass(ip, main_roi, cal, 0.0);
		IJ.log("Initial CoM: " + com[0] + ", " + com[1]);
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Recalculate the CoM after subtracting the background using an
		// annular region of interest centred on the centre of mass that has an
		// internal diameter of object diameter + 10 mm and an external
		// diameter of 2 x object diameter.
		double inner_dia = obj_dia_pix + (10.0 / pixel_size_in_mm);
		Roi inner_roi = new OvalRoi(cal.getRawX(com[0])-inner_dia/2.0, cal.getRawY(com[1])-inner_dia/2.0, inner_dia, inner_dia);
		Roi outer_roi = new OvalRoi(cal.getRawX(com[0])-obj_dia_pix, cal.getRawY(com[1])-obj_dia_pix, obj_dia_pix*2.0, obj_dia_pix*2.0);
		Roi annular_roi = new ShapeRoi(outer_roi).xor(new ShapeRoi(inner_roi));

		// Set the annular ROI to be active on the image
		imp.setRoi(annular_roi);

		// Get the (uncalibrated) statistics for the ROI
		ImageStatistics annular_stats = annular_roi.getStatistics();

		// Obtain the calibrated mean pixel value
		double bgd_mean = cal.getCValue(annular_stats.mean);

		// Recalculate the CoM with background subtracted
		imp.setRoi(outer_roi);
		com = centreOfMass(ip, outer_roi, cal, -bgd_mean);
		IJ.log("New CoM: " + com[0] + ", " + com[1]);

		// Add a marker to the image showing the location of the CoM
		Roi com_marker = new PointRoi(cal.getRawX(com[0]), cal.getRawY(com[1]));

		// Recentre the outer_roi using the new centre of mass.
		outer_roi.setLocation(cal.getRawX(com[0])-obj_dia_pix, cal.getRawY(com[1])-obj_dia_pix);

		// Display the annular ROI on the image
		Overlay overlay = new Overlay(annular_roi);
		overlay.setStrokeColor(Color.green);
		overlay.add(outer_roi);
		overlay.add(com_marker);
		imp.setOverlay(overlay);
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Obtain the raw edge spread function (ESF)
		int i=0;
		int num_points = outer_roi.getContainedPoints().length;
		double[] raw_esf_pos = new double[num_points];
		double[] raw_esf_val = new double[num_points];
		for (Point p : outer_roi.getContainedPoints()) {
		    raw_esf_pos[i] = calculateDistanceBetweenPoints(com[0], com[1], cal.getX(p.x), cal.getY(p.y));
			raw_esf_val[i] = ip.getPixelValue(p.x, p.y); // getPixelValue includes calibration
            i++;
		}
		// Sort the raw ESF arrays into ascending position order using Apache Commons Math 3.6 API
		// http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/util/MathArrays.html#sortInPlace(double[],%20double[]...)
		MathArrays.sortInPlace(raw_esf_pos, raw_esf_val);

		// Plot the raw esf
		Plot plot_esf = new Plot("Raw ESF", "Distance", "Pixel value");
		plot_esf.add("dot", raw_esf_pos, raw_esf_val);
		plot_esf.show();
		//---------------------------------------------------------------------


		//---------------------------------------------------------------------
		// Interpolate the raw ESF to a regularly-spaced array.
		// Work out the distance scale increment.
		double rebinned_sample_inc = pixel_size_in_mm / pixel_reduction_factor;

		// Define the starting position, with an offset of one sample increment
		// and work out the number of samples.
		double start_pos = StatUtils.min(raw_esf_pos) + rebinned_sample_inc;
		double num_samples = Math.floor( (StatUtils.max(raw_esf_pos) - start_pos) / rebinned_sample_inc );

		// Populate an array containing the regularly-spaced locations for the
		// ESF.
		double current_pos = start_pos;
		double[] esf_rebinned_pos = new double[(int) num_samples];
		for (i=0; i<num_samples; i++) {
			esf_rebinned_pos[i] = current_pos;
			current_pos += rebinned_sample_inc;
		}

		// Work out the value of the ESF at the rebinned positions using the Apache Commons Math 3.6 API
		// local regression algorithm. See:
		// http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/analysis/interpolation/LoessInterpolator.html
		UnivariateInterpolator loess_interpolator = new LoessInterpolator(loess_bandwidth, loess_robustness);
		UnivariateFunction loess_function = loess_interpolator.interpolate(raw_esf_pos, raw_esf_val);
		double[] esf_rebinned_val = new double[(int) num_samples];
		for (i=0; i<num_samples; i++) {
			esf_rebinned_val[i] = loess_function.value(esf_rebinned_pos[i]);
		}

		// Add the rebinned ESF to the existing ESF plot as a red line.
		plot_esf.setColor("Red");
		plot_esf.add("line", esf_rebinned_pos, esf_rebinned_val);
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
		double[] lsf_val = new double[(int) num_samples];
		for (i=0; i<num_samples; i++) {
			DerivativeStructure xDS = new DerivativeStructure(1, 1, 0, esf_rebinned_pos[i]);
			DerivativeStructure yDS = completeF.value(xDS);
			lsf_val[i] = yDS.getPartialDerivative(1);
		}

		// Plot the LSF
		Plot plot_lsf = new Plot("LSF", "Distance", "Value");
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

		// Calculate the TTF of the (float) LSF values
		FHT fht = new FHT();
		float[] ttf_val = fht.fourier1D(lsf_val_float, FHT.NO_WINDOW);

		// Calculate the normalised TTF (nTTF).
		// Using double data type so they can be used directly with Plot.
		double[] nttf_val = new double[ttf_val.length];
		double[] freq_scale = new double[ttf_val.length];
		for (i=0; i<ttf_val.length; i++) {
			// ttf_val[0] is the DC component so don't use it to normalise,
			// use ttf_val[1] instead.
			nttf_val[i] = (double) (ttf_val[i] / ttf_val[1]);

			// The frequencies are calculated according to the documentation
			// under the "Returns" section of this page:
			// https://imagej.nih.gov/ij/developer/api/ij/process/FHT.html#fourier1D-float:A-int-
			freq_scale[i] = (double) i / (2.0 * ttf_val.length * rebinned_sample_inc);
		}

		// Plot the nTTF
		Plot plot_ttf = new Plot("nTTF", "Frequency (mm-1)", "nTTF");
		plot_ttf.setColor("Red");
		plot_ttf.add("line", freq_scale, nttf_val);
		plot_ttf.show();
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
}
