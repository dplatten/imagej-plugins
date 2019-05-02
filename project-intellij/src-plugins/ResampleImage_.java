// ImageJ plugin to resample an image to a lower number of counts.
//
// David Platten, 2019
//
import ij.*;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;

public class ResampleImage_ implements PlugInFilter {
    private ImagePlus imp;
    private Calibration cal;
    private int mask_value = -2048;
    private float p = 0.5f;
    private boolean rescale_values = false;

    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        this.cal = imp.getCalibration();

        GenericDialog gd = new GenericDialog("Processing options");
        gd.addNumericField("Resampling factor (0.0 - 1.0", this.p, 2);
        gd.addNumericField("Enter mask value (pixels with this value will be ignored)", this.mask_value, 0);
        gd.addCheckbox("Re-scale results to maintain mean values?", this.rescale_values);

        gd.showDialog();
        if(gd.wasCanceled()) {
            IJ.error("Plugin cancelled");
            return(0);
        }

        this.p = (float) gd.getNextNumber();
        this.mask_value = (int) gd.getNextNumber();
        this.rescale_values = gd.getNextBoolean();

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