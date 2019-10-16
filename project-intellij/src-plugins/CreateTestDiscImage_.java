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