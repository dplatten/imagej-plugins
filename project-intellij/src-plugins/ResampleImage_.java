// ImageJ plugin to resample an image to a lower number of counts.
//
// David Platten, 2019
//
//package ij.plugin;
import ij.*;
import ij.gui.*;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;

public class ResampleImage_ implements PlugInFilter {
    protected ImagePlus imp;
    private Calibration cal;
    private int mask_value = -2048;
    private float p = 0.5f;

    public int setup(String arg, ImagePlus imp) {
    	this.imp = imp;
    	this.cal = imp.getCalibration();

        if ((arg != null) && (arg.length() > 0)) {
            this.p = Float.parseFloat(arg);
            this.mask_value = (int) IJ.getNumber("Enter mask value", this.mask_value);
        }
        else {
            this.p = (float) IJ.getNumber("Enter resampling fraction", this.p);
            this.mask_value = (int) IJ.getNumber("Enter mask value", this.mask_value);
        }

    	return DOES_ALL;
    }

    public void run(ImageProcessor ip) {
    	// Obtain the width and height of the image.
    	int width  = ip.getWidth();
    	int height = ip.getHeight();

        // Create a duplicate of the image
    	ImagePlus img_copy = new Duplicator().run(this.imp);
        img_copy.show();

    	// Main loop that steps through each pixel in the duplicate image
        ImageProcessor img_copy_ip = img_copy.getProcessor();
    	for(int y=0; y<height; y++) {
            for(int x=0; x<width; x++) {
                int current_val = img_copy_ip.getPixel(x, y);
                if (this.cal.getCValue(current_val) != this.mask_value) {
                    int count = 0;
                    for (int n = 0; n < current_val; n++) {
                        if (Math.random() < this.p) count++;
                    }
                    // Set this pixel location to be the value of count
                    img_copy_ip.putPixel(x, y, (int) (count / this.p));
                }
    		}
            img_copy.updateAndDraw();

            // Update the progress bar
    		float progress = (float)y / (float)height;
        	IJ.showStatus("Calculating column " + y + " of " + height);
    		IJ.showProgress(progress);
    	}

    	// Auto window and level a bit.
    	IJ.run("Enhance Contrast", "saturated=0.35");
    }
}