// ImageJ plugin to resample an image to a lower number of counts.
//
// David Platten, 2019
//
//package ij.plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.NewImage;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

public class AdjustSD_ implements PlugInFilter {
    private ImagePlus imp;
    private float p;

    public int setup(String arg, ImagePlus imp) {
    	this.imp = imp;

        if ((arg != null) && (arg.length() > 0)) {
            p = Float.parseFloat(arg);
        }
        else {
            p = (float) IJ.getNumber("Enter resampling fraction", 0.5);
        }

    	return DOES_ALL;
    }

    public void run(ImageProcessor ip) {
    	// Initiate some variables for the main loop
    	int x, y; //, count, n;
    	double mean; //, stdev;

        float progress;	// Used to update the ImageJ progress bar

        // Obtain the calibration properties of the original image
        Calibration cal = imp.getCalibration();
        double[] cal_coeffs = cal.getCoefficients();

        // If the user has set a region of interest then obtain it's properties.
        // If there is no ROI set then the routine is exited.
        Roi main_roi = imp.getRoi();
        if(main_roi != null) {
            ImageStatistics stats = main_roi.getStatistics();
            mean = stats.umean;
            //stdev = stats.stdDev;
        }
        else {
            return;
        }

    	// Obtain the width and height of the image.
    	int width  = ip.getWidth();
    	int height = ip.getHeight();

    	// Create a float array to contain the resampled data
        float[][] results_pixels = new float[width][height];

    	// Main loop that steps through each pixel in the image
    	for(y=0; y<height; y++) {
            for(x=0; x<width; x++) {

                // Set this pixel location to be the value of count
                results_pixels[x][y] = (float)cal_coeffs[0] + (float)cal_coeffs[1] * ((float)mean + ((ip.getPixel(x, y) - (float)mean) * p));
    		}
            // Update the progress bar
    		progress = (float)y / (float)height;
        	IJ.showStatus("Calculating column " + y + " of " + height);
    		IJ.showProgress(progress);
    	}

        // Create a new image and set the pixel values to be the contents of results_pixels
    	String image_title = "Resampled to " + p + " of original counts";
    	ImagePlus results_image = NewImage.createFloatImage(image_title, width, height, 1, NewImage.FILL_BLACK);
    	ImageProcessor results_ip = results_image.getProcessor();
        results_ip.setFloatArray(results_pixels);

    	// Display the results image.
    	results_image.show();
    	results_image.updateAndDraw();

    	// Auto window and level a bit.
    	IJ.run("Enhance Contrast", "saturated=0.35");
    }
}