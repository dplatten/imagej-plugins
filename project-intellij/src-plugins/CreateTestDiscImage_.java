// ImageJ plugin to create a disc image of a user-defined diameter
// within an image of user-defined dimensions with a user-defined
// pixel size.

import ij.IJ;
import ij.plugin.PlugIn;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import static ij.plugin.filter.PlugInFilter.DOES_ALL;

/**
 *
 */
public class CreateTestDiscImage_ implements PlugIn {

    /** Placeholder method at the moment
     *
     * @return DOES_ALL
     */
    public int setup() {
        return DOES_ALL;
    }


    /** Method to run the plugin
     *
     * @param arg not used at the moment
     */
    public void run(String arg) {

        GenericDialog gd = new GenericDialog("Create test disc image");
        gd.addNumericField("Disc diameter (mm)", 3.0, 1);
        gd.addNumericField("Image width and height (pixels)", 512, 1);
        gd.addNumericField("Image pixel size (mm)", 0.2, 2);
        gd.showDialog();
        if(gd.wasCanceled()) {
            return; // Return to the calling method
        }

        double disc_diameter = gd.getNextNumber();
        int image_pixel_width = (int) gd.getNextNumber();
        double image_pixel_size_mm = gd.getNextNumber();

        DetectabilityIndex_ DI_methods = new DetectabilityIndex_();
        Calibration cal = new Calibration();
        ImagePlus testDiscImage = DI_methods.createIdealDiscImage(disc_diameter, image_pixel_size_mm, cal, image_pixel_width);
        IJ.run(testDiscImage, "Properties...", "unit=mm pixel_width=" + image_pixel_size_mm + " pixel_height=" + image_pixel_size_mm);
        testDiscImage.show();
    }
}