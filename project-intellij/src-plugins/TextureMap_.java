

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
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;
import java.util.Arrays;

/**
 * ImageJ plugin to calculate a variance or signal-to-noise map of a
 * digital x-ray image.
 * <br><br>
 * The routine will run on a user-selected rectangular region of
 * interest (ROI). If no ROI is selected then the whole image is
 * selected and analysed.
 * <br><br>
 * It breaks the large ROI into smaller square ROIs, and calculates
 * the chosen statistical parameter for each region. The user is
 * asked to supply the side-length of the smaller ROIs (in mm).
 * <br><br>
 * The routine calculates variance as:
 * <br><br>
 * sum[x - mean(x)] / n
 * <br><br>
 * The routine calculates standard deviation as:
 * <br><br>
 * sqrt(variance)
 * <br><br>
 * The routine calculates signal-to-noise as:
 * <br><br>
 * mean(x) / standard deviation(x)
 */
public class TextureMap_ implements PlugInFilter {
	private ImagePlus imp;
	// It is important that the constants below are in alphabetical order
	// so that the binarySearch that is used later works correctly.
	private static final int  SNR=0, STANDARD_DEV=1, VARIANCE=2;
	private int processingType = VARIANCE;
	private String processingText;


	/**
	 * Setup method that is automatically called if the plugin is triggered
	 * from an ImageJ menu.
	 *
	 * @param arg Used if calling the plugin from another command to set the texture type. Options are <i>variance</i>,
	 *            <i>snr</i> and <i>standard_dev</i>. If <i>arg</i> is empty then the user is presented with a dialogue
	 *            box to choose the required processing type.
	 * @param imp The ImagePlus of the currently selected image. This is automatically passed to the plugin by ImageJ.
	 * @return A value informing ImageJ what type of images this plugin can handle. Set to <i>DOES_ALL</i>
	 */
	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		if 		(arg.equals("variance")) 		processingType = VARIANCE;
		else if (arg.equals("snr")) 			processingType = SNR;
		else if (arg.equals("standard_dev")) 	processingType = STANDARD_DEV;
		else {
			GenericDialog gd = new GenericDialog("Which type of processing?");
			// It is important that the items in the "items" array below
			// appear in alphabetical order (the same order as the constants
			// of the same name that are defined earlier).
			String[] items = {"snr", "std_dev", "variance"};
			gd.addChoice("Processing:", items, "variance");
			gd.showDialog();
			if(gd.wasCanceled()) {
				return(DONE);
			}
			processingType = Arrays.binarySearch(items, gd.getNextChoice());
		}
		return DOES_ALL;
	}


	/**
	 * Method that is automatically called after the <i>setup</i> method when the plugin is triggered from an ImageJ
	 * menu. The user is prompted for a sub-sample size (mm). The routine then calculates the chosen statistic in a
	 * square ROI of the supplied size positioned in the top left hand corner of the image. The calculated value is
	 * used as the pixel value of the corresponding position in a results image. The ROI is then moved in a raster
	 * fashion over the whole image, generating the statistic for each location and assigning this value to the
	 * corresponding location in the results image. The results image values are scaled to enhance contrast and
	 * displayed on the screen.
	 *
	 * @param ip The ImageProcessor of the currently selected image.
	 *           This is automatically passed to the plugin by ImageJ.
	 */
	public void run(ImageProcessor ip) {
		float progress;	// Used to update the ImageJ progress bar

		// Ask the user what size blocks they want the image broken into.
		float sub_sample_size_in_mm = (float) IJ.getNumber("Enter sub-sample size in mm", 2);

		// Obtain the pixel size of the image. This will return 1.0 if
		// ImageJ isn't able to read the pixel size from the DICOM tags.
		float pixel_size_in_mm = (float) imp.getFileInfo().pixelWidth;

		// Work out the block size in terms of pixels rather than mm. This
		// is rounded down to the nearest whole number.
		int sub_sample_size = (int) Math.floor( sub_sample_size_in_mm / pixel_size_in_mm);

		// If the user has set a rectangular ROI then obtain it's properties.
		// If there is no ROI set then the whole image is used.
		Roi main_roi = imp.getRoi();
		if(main_roi == null) {
			imp.setRoi(0, 0, ip.getWidth(), ip.getHeight());
			main_roi = imp.getRoi();
		}

		// Obtain the width and height of the ROI.
		int width  = main_roi.getBounds().width;
		int height = main_roi.getBounds().height;

		// Assign dimensions to the texture_map results array. These
		// dimensions are smaller than width and height to allow for
		// the sub_sample_size.
		double[][] texture_map;
		int x_size = width  - sub_sample_size;
		int y_size = height - sub_sample_size;
		texture_map = new double[x_size][y_size];

		// Obtain the start and end x and y positions of the ROI
		int start_x = main_roi.getBounds().x;
		int start_y = main_roi.getBounds().y;
		int end_x   = start_x + x_size;
		int end_y   = start_y + y_size;

		// Obtain the pixel values of the ROI as an array of floating point
		// numbers.
		float[][] pixels = ip.getFloatArray();

		// Initiate some variables for the main loop
		float total_texture=0;
		long start_time = System.currentTimeMillis();

		// Main loop that steps through the large ROI
		for(int x=start_x; x<end_x; x++) {

			int sub_end_x = x + sub_sample_size;

			for(int y=start_y; y<end_y; y++) {

				int sub_end_y = y + sub_sample_size;
				float mean, variance;

				if(y==start_y) {
					// Calculate the mean of the current block
					float total_value = 0;
					for(int sub_x=x; sub_x<sub_end_x; sub_x++) {
						for(int sub_y=y; sub_y<sub_end_y; sub_y++) {
							total_value += pixels[sub_x][sub_y];
						}
					}
					mean = total_value / (sub_sample_size*sub_sample_size);

					// Calculate the variance of the current block
					float total_diff_from_mean_sqr = 0;
					for(int sub_x=x; sub_x<sub_end_x; sub_x++) {
						for(int sub_y=y; sub_y<sub_end_y; sub_y++) {
							float value = pixels[sub_x][sub_y];
							total_diff_from_mean_sqr += (value - mean) * (value - mean);
						}
					}
					variance = total_diff_from_mean_sqr / (sub_sample_size*sub_sample_size);
				}

				else {
					float total_value = 0;
					for(int sub_x=x; sub_x<sub_end_x; sub_x++) {
						total_value -= pixels[sub_x][y-1];
						total_value += pixels[sub_x][sub_end_y-1];
					}
					mean = total_value / (sub_sample_size*sub_sample_size);

					// Calculate the variance of the current block
					float total_diff_from_mean_sqr = 0;
					for(int sub_x=x; sub_x<sub_end_x; sub_x++) {
						for(int sub_y=y; sub_y<sub_end_y; sub_y++) {
							float value = pixels[sub_x][sub_y];
							total_diff_from_mean_sqr += (value - mean) * (value - mean);
						}
					}
					variance = total_diff_from_mean_sqr / (sub_sample_size*sub_sample_size);
				}

				// Calculate the appropriate statistic
				double current_texture=0;
				if(processingType == VARIANCE) {
					processingText = "Variance";
					current_texture = variance;
				}
				else if(processingType == SNR) {
					processingText = "SNR";
					if(variance != 0.0) current_texture = mean / Math.sqrt(variance);
				}
				else if(processingType == STANDARD_DEV) {
					processingText = "Standard dev";
					current_texture = Math.sqrt(variance);
				}

				// Write the calculated statistic into the results array.
				texture_map[x-start_x][y-start_y] = current_texture;
				total_texture += current_texture;
			}
			progress = (float)(x-start_x) / (float)x_size;
			IJ.showStatus("Calculating row " + (x-start_x) + " of " + x_size);
			IJ.showProgress(progress);
		}

		IJ.log("Elapsed time is " + (System.currentTimeMillis() - start_time) / 1000.0 + " s");

		// Create new Short (16-bit signed integer) image to display the results. All numerical variables are
		// signed in Java. I'd prefer to use a Byte (8-bit signed integer), but in Java Bytes
		// can have values from -128 to +127 which is no good to me. So I'm using a
		// Short instead. This can have values from -32768 to +32767.
		String image_title = processingText + " " + sub_sample_size + " pixel sub-samples";
		ImagePlus results_image = NewImage.createShortImage(image_title, x_size, y_size, 1, NewImage.FILL_BLACK);
		ImageProcessor results_ip = results_image.getProcessor();
		short[] results_pixels = (short[]) results_ip.getPixels();

		// Step through the results array, scaling the values using
		// the known maximum and minimum figures and then allocating the
		// value to the appropriate pixel in the results image.
		double mean_texture = total_texture / (double) ( x_size * y_size);
		for(int y=0; y<=y_size-1; y++) {
			int offset = y * x_size;
			for(int x=0; x<=x_size-1; x++) {
				int i = offset + x;
				results_pixels[i] = (short) (texture_map[x][y] / (mean_texture*3) * (double) Short.MAX_VALUE);
			}
		}

		// Display the results image.
		results_image.show();
		results_image.updateAndDraw();

		// Auto window and level a bit.
		IJ.run("Enhance Contrast", "saturated=0.35");
	}
}
