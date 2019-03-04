// ImageJ plugin to calculate a variance or signal-to-noise map of a
// digital x-ray image.
//
// The routine will run on a user-selected rectangular region of
// interest (ROI). If no ROI is selected then the whole image is
// selected and analysed.
//
// It breaks the large ROI into smaller square ROIs, and calculates
// the chosen statistical parameter for each region. The user is
// asked to supply the side-length of the smaller ROIs (in mm).
//
// The routine calculates variance as:
// 		sum[x - mean(x)] / n
//
// The routine calculates standard deviation as:
// 		sqrt(variance)
//
// The routine calculates signal-to-noise as:
// 		mean(x) / standard deviation(x)
//
//package ij.plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
//import ij.gui.NewImage;
import ij.gui.Roi;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.measure.Calibration;
import ij.gui.Plot;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import flanagan.math.Gradient;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.MathArrays;
import ij.process.FHT;

import java.awt.*;
//import java.util.Arrays;

public class TaskTransferFunction_ implements PlugInFilter {
    protected ImagePlus imp;
    private double pixel_reduction_factor = 5.0;
    private double loess_bandwidth = 0.25;
    private int loess_robustness = 4;

    public int setup(String arg, ImagePlus imp) {
    	this.imp = imp;

		if ((arg != null) && (arg.length() > 0)) {
			pixel_reduction_factor = Float.parseFloat(arg);
		}
		else {
    		 GenericDialog gd = new GenericDialog("Processing options");
    		 gd.addNumericField("Pixel reduction factor", 10.0, 0);
    		 gd.addNumericField("LOESS bandwidth", 0.02, 2);
    		 gd.addNumericField("LOESS robustness", 10, 0);
    		 gd.showDialog();
    		 if(gd.wasCanceled()) {
    			 IJ.error("Plugin cancelled");
    			 return(0);
    		 }
    		 pixel_reduction_factor = gd.getNextNumber();
    		 loess_bandwidth = gd.getNextNumber();
    		 loess_robustness = (int) gd.getNextNumber();
    	}
    	return DOES_ALL;
    }

    public void run(ImageProcessor ip) {
    	float progress = 0;	// Used to update the ImageJ progress bar

    	// Obtain the pixel size of the image. This will return 1.0 if
    	// ImageJ isn't able to read the pixel size from the DICOM tags.
		Calibration cal = imp.getCalibration();
		double pixel_size_in_mm = cal.pixelWidth;

    	// If the user has set a rectangular ROI then obtain it's properties.
    	// If there is no ROI set then the whole image is used.
    	Roi main_roi = imp.getRoi();
    	if(main_roi == null) {
    		imp.setRoi(0, 0, ip.getWidth(), ip.getHeight());
			main_roi = imp.getRoi();
    	}

    	// Obtain the centre of mass of the roi
		ResultsTable rt = new ResultsTable();
		Analyzer an = new Analyzer(imp, Analyzer.CENTER_OF_MASS, rt);
		an.measure();
		double x_com = rt.getValue("XM", rt.size()-1) - pixel_size_in_mm/2.0;
		double y_com = rt.getValue("YM", rt.size()-1) - pixel_size_in_mm/2.0;
		IJ.log("CoM: " + x_com + ", " + y_com);

		// Obtain the raw edge spread function
		int num_points = main_roi.getContainedPoints().length;
		double[] raw_esf_pos = new double[num_points];
		double[] raw_esf_val = new double[num_points];
		int i = 0;
		for (Point p : main_roi.getContainedPoints()) {
		    raw_esf_pos[i] = calculateDistanceBetweenPoints(x_com, y_com, cal.getX(p.x), cal.getY(p.y));
			raw_esf_val[i] = ip.getPixelValue(p.x, p.y); // getPixelValue includes calibration
            i++;
			//IJ.log("x, y: " + cal.getX(p.x) + ", " + cal.getY(p.y));
		}
		// Sort the raw esf arrays into ascending position order
		MathArrays.sortInPlace(raw_esf_pos, raw_esf_val);

		// Plot the raw esf
		Plot plot = new Plot("Raw ESF", "Distance", "Pixel value");
		plot.add("dot", raw_esf_pos, raw_esf_val);
		plot.show();

		// Interpolate the raw esf to regular spacing
		double rebinned_sample_inc = pixel_size_in_mm / pixel_reduction_factor;
		double start_pos = StatUtils.min(raw_esf_pos) + rebinned_sample_inc;
		double num_samples = Math.floor( (StatUtils.max(raw_esf_pos) - start_pos) / rebinned_sample_inc );
		double current_pos = start_pos;
		double[] esf_rebinned_val = new double[(int) num_samples];
		double[] esf_rebinned_pos = new double[(int) num_samples];
		for (i=0; i<num_samples; i++) {
			esf_rebinned_pos[i] = current_pos;
			current_pos += rebinned_sample_inc;
		}

		UnivariateInterpolator interpolator = new LinearInterpolator();
		UnivariateFunction function = interpolator.interpolate(raw_esf_pos, raw_esf_val);
		for (i=0; i<num_samples; i++) {
			esf_rebinned_val[i] = function.value(esf_rebinned_pos[i]);
		}
		// Plot the interpolated esf
		plot.setColor("Red");
		plot.add("line", esf_rebinned_pos, esf_rebinned_val);

		UnivariateInterpolator loess_interpolator = new LoessInterpolator(loess_bandwidth, loess_robustness);
		UnivariateFunction loess_function = loess_interpolator.interpolate(raw_esf_pos, raw_esf_val);
		double[] esf_rebinned_loess_val = new double[(int) num_samples];
		for (i=0; i<num_samples; i++) {
			esf_rebinned_loess_val[i] = loess_function.value(esf_rebinned_pos[i]);
		}
		// Plot the interpolated loess esf
		plot.setColor("Green");
		plot.add("line", esf_rebinned_pos, esf_rebinned_loess_val);


		// Differentiate the esf to obtain the line spread function
		Gradient gg = new Gradient(esf_rebinned_pos, esf_rebinned_val);
		double[] lsf_val = gg.splineDeriv_1D_array();
		// Plot the lsf
		Plot plot_lsf = new Plot("LSF", "Distance", "Value");
		plot_lsf.add("line", esf_rebinned_pos, lsf_val);
		plot_lsf.show();

		// Differentiate the loess esf to obtain the line spread function
		Gradient gg_loess = new Gradient(esf_rebinned_pos, esf_rebinned_loess_val);
		double[] lsf_loess_val = gg_loess.splineDeriv_1D_array();
		// Plot the loess lsf
		plot_lsf.setColor("Red");
		plot_lsf.add("line", esf_rebinned_pos, lsf_loess_val);


		// Fourier transform the lsf to provide a ttf (mtf)
		float[] lsf_val_float = new float[lsf_val.length];
		for (i=0; i<lsf_val.length; i++) {
			lsf_val_float[i] = (float) lsf_val[i];
		}


		FHT fht = new FHT();
		float[] ttf_val = fht.fourier1D(lsf_val_float, FHT.NO_WINDOW);
		// The first element of ttf_val is the DC component - ignore
		double[] nttf_val_double = new double[ttf_val.length-1];

		// Work out frequency increment of MTF
		double freq_inc = 1.0 / (StatUtils.max(esf_rebinned_pos) - StatUtils.min(esf_rebinned_pos));
		double[] freq_scale = new double[nttf_val_double.length];

		for (i=1; i<nttf_val_double.length; i++) {
			nttf_val_double[i] = (double) (ttf_val[i] / ttf_val[1]);
			freq_scale[i] = ((double)i-1.0) * freq_inc;
		}


		// Plot the ttf
		Plot plot_ttf = new Plot("nTTF", "Frequency", "Value");
		plot_ttf.add("line", freq_scale, nttf_val_double);
		plot_ttf.show();

		// Fourier transform the loess lsf to provide a ttf (mtf)
		float[] lsf_loess_val_float = new float[lsf_loess_val.length];
		for (i=0; i<lsf_loess_val.length; i++) {
			lsf_loess_val_float[i] = (float) lsf_loess_val[i];
		}

		float[] ttf_loess_val = fht.fourier1D(lsf_loess_val_float, FHT.NO_WINDOW);
		// The first element of ttf_val is the DC component - ignore
		double[] nttf_loess_val_double = new double[ttf_loess_val.length-1];
		for (i=1; i<nttf_loess_val_double.length; i++) {
			nttf_loess_val_double[i] = (double) (ttf_loess_val[i] / ttf_loess_val[1]);
			IJ.log("TTF: " + ttf_loess_val[i] + ", " + nttf_loess_val_double[i]);
		}

		// Plot the loess ttf
		plot_ttf.setColor("Red");
		plot_ttf.add("line", freq_scale, nttf_loess_val_double);
	}

    public double calculateDistanceBetweenPoints(
            double x1,
            double y1,
            double x2,
            double y2) {
        return Math.sqrt( (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }
}
