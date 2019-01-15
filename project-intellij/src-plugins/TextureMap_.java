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
import ij.*;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;
import java.util.Arrays;

public class TextureMap_ implements PlugInFilter {
    protected ImagePlus imp;
    // It is important that the constants below are in alphabetical order
    // so that the binarySearch that is used later works correctly.
    public static final int  SNR=0, STANDARD_DEV=1, VARIANCE=2;
    private int processingType = VARIANCE;
    private String processingText;

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
    			 IJ.error("Plugin cancelled");
    			 return(0);
    		 }
    		 int index = Arrays.binarySearch(items, gd.getNextChoice());
    		 processingType = index;
    	}
    	return DOES_ALL;
    }

    public void run(ImageProcessor ip) {
    	float progress = 0;	// Used to update the ImageJ progress bar

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
    	double current_texture=0, total_texture=0, total_value=0, mean=0, previous_mean=0;
    	int x, y, sub_x, sub_y, sub_end_x, sub_end_y;
    	double total_diff_from_mean_sqr=0, value=0, variance=0;

    	// Main loop that steps through the large ROI
    	for(x=start_x; x<end_x; x++) {

    		sub_end_x = x + sub_sample_size;

			for(y=start_y; y<end_y; y++) {

    			sub_end_y = y + sub_sample_size;

    			if(y==start_y) {
    				// Calculate the mean of the current block
    				total_value = 0;
    				for(sub_x=x; sub_x<sub_end_x; sub_x++) {
    					for(sub_y=y; sub_y<sub_end_y; sub_y++) {
    						total_value += (double) pixels[sub_x][sub_y];
    					}
    				}
    				mean = total_value / ((double)sub_sample_size*(double)sub_sample_size);

    				// Calculate the variance of the current block
    				total_diff_from_mean_sqr = 0;
    				value = 0;
    				for(sub_x=x; sub_x<sub_end_x; sub_x++) {
    					for(sub_y=y; sub_y<sub_end_y; sub_y++) {
    						value = (double) pixels[sub_x][sub_y];
    						total_diff_from_mean_sqr += (value - mean) * (value - mean);
    					}
    				}
    				variance = total_diff_from_mean_sqr / ((double)sub_sample_size*(double)sub_sample_size);
    			}

    			else {
    				previous_mean = mean;

    				for(sub_x=x; sub_x<sub_end_x; sub_x++) {
    					total_value -= (double) pixels[sub_x][y-1];
    					total_value += (double) pixels[sub_x][sub_end_y-1];
    				}
    				mean = total_value / ((double)sub_sample_size*(double)sub_sample_size);

    				// Calculate the variance of the current block
    				total_diff_from_mean_sqr = 0;
    				value = 0;
    				for(sub_x=x; sub_x<sub_end_x; sub_x++) {
    					for(sub_y=y; sub_y<sub_end_y; sub_y++) {
    						value = (double) pixels[sub_x][sub_y];
    						total_diff_from_mean_sqr += (value - mean) * (value - mean);
    					}
    				}
    				variance = total_diff_from_mean_sqr / ((double)sub_sample_size*(double)sub_sample_size);
    			}


    			// Calculate the appropriate statistic
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
    		progress = ((float)x-(float)start_x) / (float)x_size;
        	IJ.showStatus("Calculating row " + (x-start_x) + " of " + x_size);
    		IJ.showProgress(progress);
    	}

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
		int offset, i;
		for(y=0; y<=y_size-1; y++) {
    		offset = y * x_size;
    		for(x=0; x<=x_size-1; x++) {
    			i = offset + x;
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
