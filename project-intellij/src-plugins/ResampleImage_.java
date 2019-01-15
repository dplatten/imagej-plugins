// ImageJ plugin to resample an image to a lower number of counts.
//
// David Platten, 2019
//
//package ij.plugin;
import ij.*;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;

public class ResampleImage_ implements PlugInFilter {
    protected ImagePlus imp;
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
    	int x, y, count, n;

        float progress = 0;	// Used to update the ImageJ progress bar

    	// Obtain the width and height of the image.
    	int width  = ip.getWidth();
    	int height = ip.getHeight();

    	// Obtain the pixel values of the image as an array of floating point numbers
    	float[][] pixels = ip.getFloatArray();

    	// Create a float array to contain the resampled data
        float[][] results_pixels = new float[width][height];

    	// Main loop that steps through each pixel in the image
    	for(y=0; y<height; y++) {
            for(x=0; x<width; x++) {
                count = 0;
                for(n=0; n<pixels[x][y]; n++) {
                    if (Math.random() < p) count++;
                }
                // Set this pixel location to be the value of count
                results_pixels[x][y] = (float) count;
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