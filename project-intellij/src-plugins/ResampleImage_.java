

/*
 * MIT License
 *
 * Copyright (c) 2019 David Platten
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

import ij.*;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;

/**
 * ImageJ plugin to resample an image to a lower number of counts. Assumes that the pixel
 * values in the image correspond to the number of photons detected at that pixel.
 * <br><br>
 * David Platten, 2019
 */
public class ResampleImage_ implements PlugInFilter {
    private ImagePlus imp;
    private Calibration cal;
    private int mask_value = -2048;
    private float p = 0.5f;
    private boolean rescale_values = false;

    /**
     * This method is automatically run by ImageJ when the plugin is called.
     *
     * @param arg Not used by this plugin
     * @param imp Passed to the routine automatically by ImageJ.
     * @return DOES_ALL, indicating that the plugin will run on any image type
     */
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        this.cal = imp.getCalibration();

        GenericDialog gd = new GenericDialog("Processing options");
        gd.addNumericField("Resampling factor (0.0 - 1.0", this.p, 2);
        gd.addNumericField("Enter mask value (pixels with this value will be ignored)", this.mask_value, 0);
        gd.addCheckbox("Re-scale results to maintain mean values?", this.rescale_values);

        gd.showDialog();
        if(gd.wasCanceled()) {
            return(DONE);
        }

        this.p = (float) gd.getNextNumber();
        this.mask_value = (int) gd.getNextNumber();
        this.rescale_values = gd.getNextBoolean();

        return DOES_ALL;
    }

    /**
     * This method is called automatically by ImageJ once the <i>setup</i> method has completed.
     *
     * @param ip The ImageProcessor that is passed to this method automatically by ImageJ
     */
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

                    // Rescale count back to the original value if required
                    if (this.rescale_values) count /= this.p;

                    // Set this pixel location to be the value of count
                    img_copy_ip.putPixel(x, y, (int) count);
                }
            }
            // Update the image each time a row is complete
            img_copy.updateAndDraw();

            // Update the progress bar
            float progress = (float)y / (float)height;
            IJ.showStatus("Calculating row " + y + " of " + height);
            IJ.showProgress(progress);
        }
    }
}