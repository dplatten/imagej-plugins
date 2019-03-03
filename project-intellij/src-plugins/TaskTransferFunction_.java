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
import ij.gui.NewImage;
import ij.gui.Roi;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.measure.Calibration;

import java.awt.*;
import java.util.Arrays;

public class TaskTransferFunction_ implements PlugInFilter {
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
		double x_com = rt.getValue("XM", rt.size()-1);
		double y_com = rt.getValue("YM", rt.size()-1);
		IJ.log("CoM: " + x_com + ", " + y_com);


		int num_points = main_roi.getContainedPoints().length;


		double[][] raw_esf = new double[num_points][2];
		int i = 0;
		for (Point p : main_roi.getContainedPoints()) {
			raw_esf[i][1] = ip.getPixelValue(p.x, p.y); // getPixelValue includes calibration
			IJ.log("x, y: " + cal.getX(p.x) + ", " + cal.getY(p.y));
		}
		//IJ.log("mean,count: "+IJ.d2s(sum/count,4)+" "+count);
    }
}
