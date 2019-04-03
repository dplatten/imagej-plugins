package QuantitativeIQ;
// ImageJ plugin to calculate the task transfer function of a circular object present in an image.
//
// The routine will run on a user-selected oval region of interest (ROI). If no
// ROI is selected then the routine exits.
//
import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.plugin.ZProjector;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.measure.Calibration;
import ij.process.ImageStatistics;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.differentiation.*;
import org.apache.commons.math3.analysis.interpolation.*;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.complex.Complex;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.lang.System;

/**
 *
 */
public class TaskTransferFunction_ implements PlugInFilter {
    private ImagePlus imp;
    private Calibration cal;

    private double obj_dia_mm = 25.0;
    private double pixel_reduction_factor = 5.0;
    private double rebin_f_inc = 0.05; // in mm^-1
    private double loess_bandwidth = 0.02;
    private int loess_robustness = 10;
    private Roi main_roi;
    private Boolean find_com = Boolean.TRUE;

    private String[] esf_proc_choices = {"Linear regression", "Local regression", "Local regression (monotonic)"};
    private static int LINEAR_REGRESSION = 0;
    private static int LOCAL_REGRESSION = 1;
    private static int MONOTONIC_LOCAL_REGRESSION = 2;
    private int esf_proc_choice_idx = MONOTONIC_LOCAL_REGRESSION;
    private String esf_proc_choice = esf_proc_choices[esf_proc_choice_idx];

    private TTF_result ttf_result;
    private CNR_result cnr_result;

    private double ttf_50;
    private double ttf_10;

    private boolean use_average_stack = false;
    private int min_stack_index;
    private int max_stack_index;

    private boolean log_results = true;

    private DecimalFormat two_dp = new DecimalFormat("0.00");
    private DecimalFormat three_dp = new DecimalFormat("0.000");


    /** Method to enable or disable using an average stack rather than a single image
     *
     * @param average_stack boolean value to set whether to use the stack average or not
     */
    void setAverageStack(boolean average_stack) {
        this.use_average_stack = average_stack;
    }


    /** Method to enable or disable writing to the ImageJ log
     *
     * @param log_contition Boolean value to set logging on or off
     */
    void setLogging(boolean log_contition) {
        this.log_results = log_contition;
    }


    // Defined to enable other plugins to access the TTF results
    TTF_result getTTFResult() {
        return ttf_result;
    }


    // Defined to enable other plugins to access the TTF results
    double getTTF50() {
        return this.ttf_50;
    }


    // Defined to enable other plugins to access the TTF results
    double getTTF10() {
        return this.ttf_10;
    }


    // Defined to enable other plugins to access the CNR results
    CNR_result getCNRResult() {
        return cnr_result;
    }


    /**
     *
     * @return
     */
    double getFreqBinInc() {
        return rebin_f_inc;
    }


    /**
     *
     * @param arg
     * @param imp
     * @return
     */
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;

        if(this.imp == null) {
            return(0);
        }

        // Obtain the calibration for the ImagePlus, as this is lost if we use a z-projection of a stack
        this.cal = this.imp.getCalibration();

        // Remove any overlays that are present
        this.imp.setOverlay(null);

        // If there is no ROI set then exit
        this.main_roi = imp.getRoi();
        if(this.main_roi == null) {
            IJ.error("Please draw an ROI around an object");
            return(0);
        }

        // Make the main_roi square with a side length that is a power of two
        // Ensure the ROI is square with a side length that is a power of two
        double roi_x = this.main_roi.getXBase();
        double roi_y = this.main_roi.getYBase();

        double N = this.main_roi.getFloatWidth();
        int M = (int) (Math.log(N) / Math.log(2.0));
        int N_power_two = (int) Math.pow(2, M);
        if (N_power_two != N) { // Make the Roi larger using the next power of 2 bigger than the current width
            N_power_two = (int) Math.pow(2, M+1);
            Roi new_roi = new Roi(this.main_roi.getXBase() - (N_power_two-N)/2.0, this.main_roi.getYBase() - (N_power_two-N)/2.0, N_power_two, N_power_two);
            this.main_roi = new_roi;
            this.imp.setRoi(this.main_roi);
            this.imp.updateAndDraw();
        }

        GenericDialog gd = new GenericDialog("Task transfer function options");
        gd.addNumericField("Object diameter (mm)", this.obj_dia_mm, 0);
        gd.addNumericField("Pixel reduction factor", this.pixel_reduction_factor, 0);
        gd.addNumericField("Frequency rebin increment (mm^-1)", this.rebin_f_inc, 2);
        gd.addCheckbox("Find centre of mass", this.find_com);
        gd.addRadioButtonGroup("ESF processing to use", this.esf_proc_choices, 3, 1, this.esf_proc_choices[this.esf_proc_choice_idx]);
        gd.addMessage("Local regression options");
        gd.addNumericField("Bandwidth", this.loess_bandwidth, 2);
        gd.addNumericField("Robustness", this.loess_robustness, 0);
        if (this.imp.isStack()) {
            this.use_average_stack = true;
            gd.addCheckbox("Run on the average stack", this.use_average_stack);
            gd.addNumericField("Min slice index", 1, 0);
            gd.addNumericField("to max slice index", this.imp.getImageStackSize(), 0);
        }
        gd.showDialog();
        if(gd.wasCanceled()) {
            return(0);
        }

        this.obj_dia_mm = gd.getNextNumber();
        this.pixel_reduction_factor = gd.getNextNumber();
        this.rebin_f_inc = gd.getNextNumber();
        this.find_com = gd.getNextBoolean();
        this.esf_proc_choice = gd.getNextRadioButton();
        this.esf_proc_choice_idx = indexOf(this.esf_proc_choices, this.esf_proc_choice);
        this.loess_bandwidth = gd.getNextNumber();
        this.loess_robustness = (int) gd.getNextNumber();
        if (this.imp.isStack()) {
            this.use_average_stack = gd.getNextBoolean();
            this.min_stack_index = (int) gd.getNextNumber();
            this.max_stack_index = (int) gd.getNextNumber();
        }

        if (this.use_average_stack) {
            ZProjector z_projector = new ZProjector(this.imp);
            z_projector.setMethod(ZProjector.AVG_METHOD);
            z_projector.setStartSlice(240);
            z_projector.setStopSlice(273);
            z_projector.doProjection();
            this.imp = z_projector.getProjection();
            this.imp.setCalibration(this.cal);
            this.imp.show();
            this.imp.setRoi(this.main_roi);
        }

        return DOES_ALL;
    }

    /**
     *
     * @param ip
     */
    public void run(ImageProcessor ip) {
        if (this.use_average_stack) {
            ip = this.imp.getProcessor();
        }


        //---------------------------------------------------------------------
        // Obtain the pixel size of the image. This will return 1.0 if ImageJ
        // isn't able to read the pixel size from the DICOM tags.
        //Calibration cal = this.imp.getCalibration();
        double pixel_size_in_mm = this.cal.pixelWidth;
        double obj_dia_pix = this.obj_dia_mm/pixel_size_in_mm;
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Determine the centre of mass (Com). Initialise the CoM to the centre
        // of the main_roi.
        double[] com = new double[2];
        com[0] = this.cal.getX(this.main_roi.getXBase() + (this.main_roi.getFloatWidth()/2.0));
        com[1] = this.cal.getY(this.main_roi.getYBase() + (this.main_roi.getFloatHeight()/2.0));

        double[] com_new = com;

        double inner_dia = obj_dia_pix + (10.0 / pixel_size_in_mm);
        Roi inner_roi = new OvalRoi(this.cal.getRawX(com[0]) - inner_dia / 2.0, this.cal.getRawY(com[1]) - inner_dia / 2.0, inner_dia, inner_dia);
        Roi outer_roi = new OvalRoi(this.cal.getRawX(com[0]) - obj_dia_pix, this.cal.getRawY(com[1]) - obj_dia_pix, obj_dia_pix * 2.0, obj_dia_pix * 2.0);
        Roi annular_roi = new ShapeRoi(outer_roi).xor(new ShapeRoi(inner_roi));

        if (this.find_com) {
            //---------------------------------------------------------------------
            // Obtain the centre of mass of the ROI set by the user.
            this.imp.setRoi(this.main_roi);
            com = centreOfMass(ip, this.main_roi, this.cal, 0.0);
            if (this.log_results) IJ.log("Initial CoM:," + this.three_dp.format(com[0]) + "," + this.three_dp.format(com[1]));
            //---------------------------------------------------------------------


            //---------------------------------------------------------------------
            // Recalculate the CoM after subtracting the background using an
            // annular region of interest centred on the centre of mass that has an
            // internal diameter of object diameter + 10 mm and an external
            // diameter of 2 x object diameter.
            inner_roi = new OvalRoi(this.cal.getRawX(com[0]) - inner_dia / 2.0, this.cal.getRawY(com[1]) - inner_dia / 2.0, inner_dia, inner_dia);
            outer_roi = new OvalRoi(this.cal.getRawX(com[0]) - obj_dia_pix, this.cal.getRawY(com[1]) - obj_dia_pix, obj_dia_pix * 2.0, obj_dia_pix * 2.0);
            annular_roi = new ShapeRoi(outer_roi).xor(new ShapeRoi(inner_roi));


            // Obtain the mean pixel value of the annular ROI
            this.imp.setRoi(annular_roi);
            ImageStatistics annular_stats = this.imp.getAllStatistics();
            double bgd_mean = annular_stats.mean;

            // Recalculate the CoM with background subtracted
            this.imp.setRoi(outer_roi);
            com = centreOfMass(ip, outer_roi, this.cal, -bgd_mean);
            if (this.log_results) IJ.log("Adjusted CoM:," + this.three_dp.format(com[0]) + "," + this.three_dp.format(com[1]));


            // Recentre the annualar roi on the new CoM and recalculate CoM  for a final time
            //annular_roi.setLocation(this.cal.getRawX(com[0]) - obj_dia_pix, this.cal.getRawY(com[1]) - obj_dia_pix);
            annular_roi.setLocation(this.cal.getRawX(com_new[0]) - obj_dia_pix, this.cal.getRawY(com_new[1]) - obj_dia_pix);
            this.imp.setRoi(annular_roi);
            annular_stats = this.imp.getAllStatistics();
            bgd_mean = annular_stats.mean;
            this.imp.setRoi(outer_roi);
            com = centreOfMass(ip, outer_roi, cal, -bgd_mean);


            main_roi.setLocation(this.cal.getRawX(com[0]) - this.main_roi.getFloatWidth()/2.0, this.cal.getRawY(com[1]) - this.main_roi.getFloatWidth()/2.0);
            com_new = centreOfMassLowFreqFilter(ip, this.main_roi, this.cal, this.obj_dia_mm, -bgd_mean);
            IJ.log("Low pass com is:," + com_new[0] + "," + com_new[1]);
        }

        if (this.log_results) IJ.log("Final CoM:," + this.three_dp.format(com[0]) + "," + this.three_dp.format(com[1]));
        com = com_new;

        // Add a marker to the image showing the location of the CoM
        Roi com_marker = new PointRoi(this.cal.getRawX(com[0]), this.cal.getRawY(com[1]));
        Overlay com_overlay = new Overlay(com_marker);
        com_overlay.setStrokeColor(Color.red);
        this.imp.setOverlay(com_overlay);

        // Recentre the outer_roi and annular_roi using the centre of mass.
        outer_roi.setLocation(this.cal.getRawX(com[0])-obj_dia_pix, this.cal.getRawY(com[1])-obj_dia_pix);
        annular_roi.setLocation(this.cal.getRawX(com_new[0]) - obj_dia_pix, this.cal.getRawY(com_new[1]) - obj_dia_pix);

        // Create an overlay to display the annular ROI and CoM on the image.
        Overlay overlay = new Overlay();
        this.imp.setOverlay(overlay);
        overlay.add(annular_roi);
        overlay.add(com_marker);
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Calculate the CNR of the object
        double obj_roi_dia = obj_dia_pix - (3.0 / pixel_size_in_mm);
        Roi obj_roi = new OvalRoi(this.cal.getRawX(com[0])-obj_roi_dia/2.0, this.cal.getRawY(com[1])-obj_roi_dia/2.0, obj_roi_dia, obj_roi_dia);
        CNR_result cnr_results = contrastToNoiseRatio(this.imp, obj_roi, annular_roi);

        if (this.log_results) IJ.log("CNR is: " + this.three_dp.format(cnr_results.cnr));
        if (this.log_results) IJ.log("Contrast is: " + this.three_dp.format(cnr_results.contrast));
        if (this.log_results) IJ.log("Noise is: " + this.three_dp.format(cnr_results.noise));

        overlay.add(obj_roi);
        overlay.setStrokeColor(Color.green);

        cnr_result = cnr_results;
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Obtain the raw edge spread function (ESF).
        int i=0;
        int num_points = outer_roi.getContainedPoints().length;
        double[] raw_esf_pos = new double[num_points];
        double[] raw_esf_val = new double[num_points];
        for (Point p : outer_roi.getContainedPoints()) {
            raw_esf_pos[i] = calculateDistanceBetweenPoints(com[0], com[1], this.cal.getX(p.x), this.cal.getY(p.y));
            raw_esf_val[i] = ip.getPixelValue(p.x, p.y); // getPixelValue includes calibration
            i++;
        }

        // There's a chance that two or more elements are at exactly the same position. This violates the monotonic
        // requirement of the interpolator used later on. So find duplicate position entries and replace with the mean
        // value at that position.
        // http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/util/MathArrays.html#unique(double[])
        double[] unique_pos = MathArrays.unique(raw_esf_pos);

        // If the length of unique_pos is less than raw_esf_pos then there must
        // be multiple data points at the same distance from the CoM.
        if (unique_pos.length < raw_esf_pos.length) {
            double[] unique_pos_vals = new double[unique_pos.length];

            for (i = 0; i < unique_pos.length; i++) {
                int count = 0;
                double sum = 0.0;
                for (int j = 0; j < raw_esf_pos.length; j++) {
                    if (raw_esf_pos[j] == unique_pos[i]) {
                        count++;
                        sum += raw_esf_val[j];
                    }
                }
                unique_pos_vals[i] = sum / count;
            }

            // Overwrite the initial position and value data with the new set
            // of data at unique positions.
            raw_esf_pos = unique_pos;
            raw_esf_val = unique_pos_vals;
        }

        // Sort the raw ESF arrays into ascending position order using Apache Commons Math 3.6 API. See:
        // http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/util/MathArrays.html#sortInPlace(double[],%20double[]...)
        MathArrays.sortInPlace(raw_esf_pos, raw_esf_val);

        // Plot the raw esf
        Plot plot_esf = new Plot("Edge spread function", "Distance (mm)", "Pixel value");
        plot_esf.add("dot", raw_esf_pos, raw_esf_val);
        String esf_labels = "Raw ESF\t";
        plot_esf.show();
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Rebin the ESF
        double rebinned_sample_inc = pixel_size_in_mm / this.pixel_reduction_factor;

        // Define the starting position of the resampled ESF, with an offset of
        // one sample increment, and work out the number of samples.
        double start_pos = StatUtils.min(raw_esf_pos) + rebinned_sample_inc;
        double num_samples = Math.floor((StatUtils.max(raw_esf_pos) - start_pos) / rebinned_sample_inc);

        // Populate an array containing the regularly-spaced locations for the ESF.
        double current_pos = start_pos;
        double[] esf_rebinned_pos = new double[(int) num_samples];
        double[] esf_rebinned_val = new double[(int) num_samples];
        for (i = 0; i < esf_rebinned_pos.length; i++) {
            esf_rebinned_pos[i] = current_pos;
            current_pos += rebinned_sample_inc;
        }

        UnivariateInterpolator interpolator;
        UnivariateFunction function;

        if (this.esf_proc_choice_idx == LINEAR_REGRESSION) {
            // Calculate resampled ESF values using linear interpolation
            interpolator = new LinearInterpolator();
            function = interpolator.interpolate(raw_esf_pos, raw_esf_val);
            //double[] esf_rebinned_val = new double[esf_rebinned_pos.length];
            for (i = 0; i < esf_rebinned_val.length; i++) {
                esf_rebinned_val[i] = function.value(esf_rebinned_pos[i]);
            }
        }
        else {
            // Calculate resampled ESF values using local regression with the Apache Commons Math 3.6 API. See:
            // http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/analysis/interpolation/LoessInterpolator.html
            // and https://en.wikipedia.org/wiki/Local_regression
            interpolator = new LoessInterpolator(this.loess_bandwidth, this.loess_robustness);
            double[] smoothed_esf = new LoessInterpolator(this.loess_bandwidth, this.loess_robustness).smooth(raw_esf_pos, raw_esf_val);
            function = interpolator.interpolate(raw_esf_pos, smoothed_esf);

            for (i = 0; i < esf_rebinned_val.length; i++) {
                esf_rebinned_val[i] = function.value(esf_rebinned_pos[i]);
            }

            if (this.esf_proc_choice_idx == MONOTONIC_LOCAL_REGRESSION) {
                // Create a monotonic ESF from the local regression ESF using brute force.

                // Check to see if ESF starts high and goes low
                if (esf_rebinned_val[0] > esf_rebinned_val[esf_rebinned_val.length - 1]) {
                    // Need to reverse the values array
                    MathArrays.sortInPlace(esf_rebinned_pos, MathArrays.OrderDirection.DECREASING, esf_rebinned_val);
                }

                for (i = 1; i < esf_rebinned_val.length; i++) {
                    if (esf_rebinned_val[i] < esf_rebinned_val[i - 1]) {
                        esf_rebinned_val[i] = esf_rebinned_val[i - 1];
                    }
                }

                // Sort the esf_rebinned_pos and associated arrays back to ascending order
                MathArrays.sortInPlace(esf_rebinned_pos, esf_rebinned_val);
            }
        }

        // Add the rebinned ESF to the existing ESF plot as a line.
        plot_esf.setColor("Red");
        plot_esf.add("line", esf_rebinned_pos, esf_rebinned_val);
        esf_labels += this.esf_proc_choice;
        plot_esf.setColor("Black");
        plot_esf.addLegend(esf_labels);
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Differentiate the ESF to obtain line spread function (LSF)
        FiniteDifferencesDifferentiator differentiator = new FiniteDifferencesDifferentiator(2, rebinned_sample_inc);;
        function = interpolator.interpolate(esf_rebinned_pos, esf_rebinned_val);
        UnivariateDifferentiableFunction complete_f = differentiator.differentiate(function);;
        double[] lsf_val = differentiate(esf_rebinned_pos, complete_f);

        // Plot the LSF.
        Plot plot_lsf = new Plot("Line spread function", "Distance (mm)", "Value");
        plot_lsf.setColor("Red");
        plot_lsf.add("line", esf_rebinned_pos, lsf_val);
        plot_lsf.setColor("Black");
        plot_lsf.addLegend(this.esf_proc_choice);
        plot_lsf.show();
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Fourier transform the LSF to obtain TTF
        TTF_result ttf_data = TTF(lsf_val, rebinned_sample_inc, pixel_size_in_mm);

        // Plot the nTTFs
        Plot plot_ttf = new Plot("Normalised TTF", "Frequency (per mm)", "nTTF");
        plot_ttf.setColor("Red");
        plot_ttf.add("line", ttf_data.freq, ttf_data.val);
        plot_ttf.setColor("Black");
        plot_ttf.addLegend(this.esf_proc_choice);
        plot_ttf.show();
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Work out TTF 50 and TTF 10
        double ttf_50 = freqAtSpecificTTF(ttf_data, 0.5);
        double ttf_10 = freqAtSpecificTTF(ttf_data, 0.1);
        if (this.log_results) IJ.log("MTF 50 and 10 are: " + this.three_dp.format(ttf_50) + ", " + this.three_dp.format(ttf_10));
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Resample the TTF data
        TTF_result resampled_ttf = resampleTTF(ttf_data, rebin_f_inc);
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Write the resampled TTF data to the ImageJ log
        if (this.log_results) logTTF(resampled_ttf, "TTF ("+this.esf_proc_choice+")");
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Remove the current ROI
        this.imp.deleteRoi();
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Set the final_results property so it can be obtained for further use
        this.ttf_result = ttf_data;
        this.ttf_50 = ttf_50;
        this.ttf_10 = ttf_10;
        //---------------------------------------------------------------------
    }


    /**
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    private double calculateDistanceBetweenPoints(double x1, double y1, double x2,double y2) {
        // Calculate the distance between two points in cartesian space.
        return Math.sqrt( (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }


    /**
     *
     * @param ip
     * @param roi
     * @param cal
     * @param pv_offset
     * @return
     */
    private double[] centreOfMass(ImageProcessor ip, Roi roi, Calibration cal, double pv_offset) {
        // Calculate the centre of mass of pixels in a Roi positioned on an
        // ImageProcessor with a specific calibration. The user can apply a
        // pixel value offset to subtract background.
        double x_wt=0;
        double y_wt=0;
        double total_wt=0;
        for (Point p : roi.getContainedPoints()) {
            double pixel_val = ip.getPixelValue(p.x, p.y);
            pixel_val += pv_offset;
            x_wt += cal.getX(p.x) * pixel_val;
            y_wt += cal.getY(p.y) * pixel_val;
            total_wt += pixel_val;
        }
        double x_com = x_wt / total_wt;
        double y_com = y_wt / total_wt;
        return new double[] {x_com, y_com};
    }


    private double[] centreOfMassLowFreqFilter(ImageProcessor ip, Roi roi, Calibration cal, double obj_diam, double pv_offset) {

        double[][] roi_pixels = getCalibratedPixelValues(ip, roi);

/*
        // ---------------------------------------------------------------------
        // Display the ROI
        ImagePlus temp = NewImage.createFloatImage("Roi pixels", roi_pixels.length, roi_pixels[0].length, 1, NewImage.FILL_BLACK);
        temp.getProcessor().setFloatArray(convertDoubleToFloat(roi_pixels));
        //IJ.run(nps_2d, "Swap Quadrants", "");
        temp.show();
        temp.updateAndDraw();
        IJ.run(temp, "Enhance Contrast", "saturated=0.35");
        IJ.run("Set... ", "zoom=800");
        // ---------------------------------------------------------------------
*/

        TwoDFFT twod_fourier_transformer = new TwoDFFT();
        Complex[][] fft = twod_fourier_transformer.transform(roi_pixels, DftNormalization.STANDARD);

/*
        // ---------------------------------------------------------------------
        // Display the ROI
        temp = NewImage.createFloatImage("Roi FFT", fft.length, fft[0].length, 1, NewImage.FILL_BLACK);
        temp.getProcessor().setFloatArray(convertComplexToFloat(fft));
        //IJ.run(nps_2d, "Swap Quadrants", "");
        temp.show();
        temp.updateAndDraw();
        IJ.run(temp, "Enhance Contrast", "saturated=0.35");
        IJ.run("Set... ", "zoom=800");
        // ---------------------------------------------------------------------
*/


        double freq_inc = 1.0 / (roi_pixels.length * cal.pixelWidth);


        // Make frequencies higher than freq_cut_off zero in the FFT.
        //double low_freq_cut_off = 0.5 / obj_diam;
        double high_freq_cut_off = 1.0 / obj_diam;

        Complex complexZero = new Complex(0.0, 0.0);

        // First the top-left quadrant, +ve i and +ve j
        for(int i=0; i<roi_pixels.length/2.0; i++) {
            for(int j=0; j<roi_pixels[0].length/2.0; j++) {
                // Calculate distance from the top-left hand corner (0,0)
                double dist = calculateDistanceBetweenPoints(0,0, i, j);
                double radial_freq = dist * freq_inc;
                if (radial_freq > high_freq_cut_off) fft[i][j] = complexZero;
            }
        }

        // Next the -i and +j quadrant
        int i_start = (int) (roi_pixels.length/2.0);
        for(int i=i_start; i<roi_pixels.length; i++) {
            for(int j=0; j<roi_pixels[0].length/2.0; j++) {
                // Calculate distance from the top-right hand corner (roi_pixels.length, 0)
                double dist = calculateDistanceBetweenPoints(roi_pixels.length, 0, i, j);
                double radial_freq = dist * freq_inc;
                if (radial_freq > high_freq_cut_off) fft[i][j] = complexZero;
            }
        }

        // Next the bottom-left quadrant, +ve i and -ve j
        int j_start = (int) (roi_pixels[0].length/2.0);
        for(int i=0; i<roi_pixels.length/2.0; i++) {
            for(int j=j_start; j<roi_pixels[0].length; j++) {
                // Calculate the distance from the bottom-left hand corner (0, roi_pixels[0].length)
                double dist = calculateDistanceBetweenPoints(0, roi_pixels[0].length, i, j);
                double radial_freq = dist * freq_inc;
                if (radial_freq > high_freq_cut_off) fft[i][j] = complexZero;
            }
        }

        // Finally the bottom-right quadrant, -ve i and -ve j
        i_start = (int) (roi_pixels.length/2.0);
        j_start = (int) (roi_pixels[0].length/2.0);
        for(int i=i_start; i<roi_pixels.length; i++) {
            for(int j=j_start; j<roi_pixels[0].length; j++) {
                // Calculate the distance from the bottom-right hand corner (roi_pixels.length, roi_pixels[0].length)
                double dist = calculateDistanceBetweenPoints(roi_pixels.length, roi_pixels[0].length, i, j);
                double radial_freq = dist * freq_inc;
                if (radial_freq > high_freq_cut_off) fft[i][j] = complexZero;
            }
        }


        // Now shift the quadrants so that zero is at the centre
        Complex[][] shifted_fft = new Complex[fft.length][fft[0].length];
        System.arraycopy(fft, (int)(fft.length/2.0), shifted_fft, 0, (int)(fft.length/2.0));
        System.arraycopy(fft, 0, shifted_fft, (int)(fft.length/2.0), (int)(fft.length/2.0));

        fft = transpose(shifted_fft);
        shifted_fft = transpose(shifted_fft);
        System.arraycopy(fft, (int)(fft.length/2.0), shifted_fft, 0, (int)(fft.length/2.0));
        System.arraycopy(fft, 0, shifted_fft, (int)(fft.length/2.0), (int)(fft.length/2.0));

        shifted_fft = transpose(shifted_fft);

/*
        // ---------------------------------------------------------------------
        temp = NewImage.createFloatImage("Roi shifted FFT", shifted_fft.length, shifted_fft[0].length, 1, NewImage.FILL_BLACK);
        temp.getProcessor().setFloatArray(convertComplexToFloat(shifted_fft));
        temp.show();
        temp.updateAndDraw();
        IJ.run(temp, "Enhance Contrast", "saturated=0.35");
        IJ.run("Set... ", "zoom=800");
        // ---------------------------------------------------------------------
*/

        // Now inverse Fourier transform the fft to obtain image space again
        fft = twod_fourier_transformer.transform(shifted_fft, DftNormalization.STANDARD, TransformType.INVERSE);

/*
        // ---------------------------------------------------------------------
        temp = NewImage.createFloatImage("Inverse FFT", fft.length, fft[0].length, 1, NewImage.FILL_BLACK);
        temp.getProcessor().setFloatArray(convertComplexToFloat(fft));
        temp.show();
        temp.updateAndDraw();
        IJ.run(temp, "Enhance Contrast", "saturated=0.35");
        IJ.run("Set... ", "zoom=800");
        // ---------------------------------------------------------------------
*/

        // Now find the centre of mass of the filtered image
        double x_wt=0;
        double y_wt=0;
        double total_wt=0;
        double roi_x_base = roi.getXBase();
        double roi_y_base = roi.getYBase();
        for (int x=0; x<fft.length; x++) {
            for (int y=0; y<fft.length; y++) {
                double pixel_val = fft[x][y].abs();
                pixel_val += pv_offset;
                x_wt += (x+roi_x_base) * pixel_val;
                y_wt += (y+roi_y_base) * pixel_val;
                total_wt += pixel_val;
            }
        }
        double x_com = cal.getX(x_wt / total_wt);
        double y_com = cal.getY(y_wt / total_wt);

        return new double[] {x_com, y_com};
    }


    /**
     *
     * @param imp
     * @param obj_roi
     * @param bgd_roi
     * @return
     */
    private CNR_result contrastToNoiseRatio(ImagePlus imp, Roi obj_roi, Roi bgd_roi) {
        // Calculate the contrast to noise ratio between the two regions of interest.
        // CNR calculated as (obj_mean - bgd_mean) / bgd_stdDev as show in equation
        // 1 of Christianson et al, https://doi.org/10.1148/radiol.15132091
        imp.setRoi(obj_roi);
        ImageStatistics obj_stats = imp.getAllStatistics();

        imp.setRoi(bgd_roi);
        ImageStatistics bgd_stats = imp.getAllStatistics();

        CNR_result result = new CNR_result();
        result.contrast = obj_stats.mean - bgd_stats.mean;
        result.noise = bgd_stats.stdDev;
        result.cnr = Math.abs(result.contrast) / result.noise;

        return result;
    }


    /**
     *
     * @param lsf
     * @param lsf_inc_mm
     * @param px_size_mm
     * @return
     */
    private TTF_result TTF(double[] lsf, double lsf_inc_mm, double px_size_mm) {
        // This uses the Fourier transform routine within the Apache Commons
        // Math 3.6 API. See:
        // http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/transform/FastFourierTransformer.html

        // Apply a Hann window before the FFT to match ImaQuest software - see
        // paper by Chen et al: http://dx.doi.org/10.1118/1.4881519
        // Also see http://download.ni.com/evaluation/pxi/Understanding%20FFTs%20and%20Windowing.pdf
        double[] hann = hannWindow(lsf.length);
        lsf = MathArrays.ebeMultiply(lsf, hann);

        FastFourierTransformer fourier_transformer = new FastFourierTransformer(DftNormalization.STANDARD);

        // FastFourierTransformer requires that the number of data elements is
        // a power of 2, so pad the LSF to the next highest power of 2 if
        // necessary.
        int N = lsf.length;
        int M = (int) (Math.log(N) / Math.log(2.0));
        int N_power_two = (int) Math.pow(2, M);

        if (N_power_two != N) { // Pad the ttf_val array with zeros to be a power of 2
            N_power_two = (int) Math.pow(2, M+1);

            // Create the padding part
            double[] zero_padding = new double[N_power_two - N];
            for (int i=0; i<zero_padding.length; i++) zero_padding[i] = 0;

            // Concatenate the padding on to the original LSF
            lsf = MathArrays.concatenate(lsf, zero_padding);
        }

        // Fourier transform the padded LSF
        Complex[] complex_ttf_results = fourier_transformer.transform(lsf, TransformType.FORWARD);

        // Work out the Nyquist frequency and how many TTF elements to plot in
        // order to get to 2 x Nyquist.
        double nyquist_freq = 1.0 / (2.0 * px_size_mm);
        double freq_inc = 1.0 / (lsf.length * lsf_inc_mm);
        int num_ttf_elements_to_inc = (int) Math.floor( (2.0 * nyquist_freq) / freq_inc);

        // Calculate the normalised TTF and the frequency scale up to 2 x Nyquist
        double[] nttf_results = new double[num_ttf_elements_to_inc];
        double[] ttf_freq_scale = new double[num_ttf_elements_to_inc];

        for (int i=0; i<num_ttf_elements_to_inc; i++) {
            nttf_results[i] = complex_ttf_results[i].abs() / complex_ttf_results[0].abs();
            ttf_freq_scale[i] = i * freq_inc;
        }

        TTF_result results = new TTF_result();
        results.freq = ttf_freq_scale;
        results.val = nttf_results;

        return results;
    }


    /**
     *
     * @param locations
     * @param diff_func
     * @return
     */
    private double[] differentiate(double[] locations, UnivariateDifferentiableFunction diff_func) {
        double[] deriv_val = new double[locations.length];
        for (int i=0; i<locations.length; i++) {
            try {
                DerivativeStructure xDS = new DerivativeStructure(1, 1, 0, locations[i]);
                DerivativeStructure yDS = diff_func.value(xDS);
                deriv_val[i] = yDS.getPartialDerivative(1);
            }
            catch (OutOfRangeException e) {
                // This is thrown for the first and last locations if a linear
                // interpolator is used because the interpolator doesn't have a
                // data point at either side of the required position to use
                // for its calculation.
            }
        }
        return deriv_val;
    }


    /**
     *
     * @param ttf
     * @param freq_inc
     * @return
     */
    private TTF_result resampleTTF(TTF_result ttf, double freq_inc) {
        UnivariateInterpolator linear_interpolator = new LinearInterpolator();
        UnivariateFunction fn_ttf_lin = linear_interpolator.interpolate(ttf.freq, ttf.val);

        ArrayList<Double> ttf_freq = new ArrayList<Double>();
        ArrayList<Double> ttf_val = new ArrayList<Double>();

        double freq = 0.0;
        while (freq < StatUtils.max(ttf.freq)) {
            try {
                ttf_val.add(fn_ttf_lin.value(freq));
                ttf_freq.add(freq);
            }
            catch (OutOfRangeException e) {
                // If the interpolator doesn't have data on either side of the
                // required position then it throws an error.
            }
            freq += freq_inc;
        }

        double[] freq_array = new double[ttf_freq.size()];
        double[] val_array = new double[ttf_val.size()];
        for (int i=0; i<freq_array.length; i++) {
            freq_array[i] = ttf_freq.get(i);
            val_array[i] = ttf_val.get(i);
        }

        TTF_result result = new TTF_result();
        result.freq = freq_array;
        result.val = val_array;

        return result;
    }


    /**
     *
     * @param ttf
     * @param ttf_label
     */
    private void logTTF(TTF_result ttf, String ttf_label) {
        IJ.log("Freq (mm^-1)," + ttf_label);

        for (int i=0; i<ttf.freq.length; i++) {
            IJ.log("" + this.two_dp.format(ttf.freq[i]) + ", " + this.three_dp.format(ttf.val[i]));
        }
    }


    /**
     *
     * @param ttf
     * @param req_ttf_val
     * @return
     */
    private double freqAtSpecificTTF(TTF_result ttf, double req_ttf_val) {
        UnivariateInterpolator linear_interpolator = new LinearInterpolator();
        UnivariateFunction fn_ttf_lin = linear_interpolator.interpolate(ttf.freq, ttf.val);

        double current_ttf = 100.0;
        double prev_ttf = 100.0;

        double freq_step = 0.01;

        double current_freq = ttf.freq[0];
        double prev_freq = ttf.freq[0];

        while (current_ttf > req_ttf_val) {
            prev_ttf = current_ttf;
            prev_freq = current_freq - freq_step;

            try {
                current_ttf = fn_ttf_lin.value(current_freq);
            }
            catch (OutOfRangeException e) {
                // If the interpolator doesn't have data on either side of the
                // required position then it throws an error. This may also
                // be thrown if there is a problem with the TTF data that means
                // this while loop never reaches the required ttf value.
                return 999.0;
            }
            current_freq += freq_step;
        }
        double[] ttf_vals = {prev_ttf, current_ttf};
        double[] ttf_freqs = {prev_freq, current_freq-freq_step};

        MathArrays.sortInPlace(ttf_vals, ttf_freqs);
        fn_ttf_lin = linear_interpolator.interpolate(ttf_vals, ttf_freqs);

        return fn_ttf_lin.value(req_ttf_val);
    }


    /**
     *
     * @param pos
     * @param val
     * @param pos_inc
     * @param pos_bin
     * @return
     */
    private ESF_result rebinESF(double[] pos, double[] val, double pos_inc, double pos_bin) {

        double max_pos = StatUtils.max(pos);
        double min_pos = StatUtils.min(pos);
        int num_intervals = (int) Math.round( (max_pos-min_pos) / pos_inc) + 1;

        double[] new_val = new double[num_intervals];
        double[] new_pos = new double[num_intervals];
        for (int i=0; i<num_intervals; i++) {
            new_pos[i] = min_pos + (i * pos_inc);
        }

        for (int i=0; i<num_intervals; i++) {
            double sum = 0;
            double count = 0;
            double lower_edge = new_pos[i] - pos_bin;
            double upper_edge = new_pos[i] + pos_bin;
            for (int j=0; j<pos.length; j++) {
                if (pos[j] >= lower_edge && pos[j] <= upper_edge) {
                    count++;
                    sum += val[j];
                }
            }
            new_val[i] = sum / count;
        }

        ESF_result results = new ESF_result();
        results.pos = new_pos;
        results.val = new_val;

        return results;
    }


    /**
     *
     * @param freq
     * @param val
     * @param f_inc
     * @param f_bin
     * @return
     */
    private TTF_result rebinTTF(double[] freq, double[] val, double f_inc, double f_bin) {

        double max_freq = StatUtils.max(freq);
        int num_intervals = (int) Math.round( max_freq / f_inc) + 1;

        double[] new_val = new double[num_intervals];
        double[] new_freq = new double[num_intervals];
        for (int i=0; i<num_intervals; i++) {
            new_freq[i] = i * f_inc;
        }

        for (int i=0; i<num_intervals; i++) {
            double sum = 0;
            double count = 0;
            double lower_edge = new_freq[i] - f_bin;
            double upper_edge = new_freq[i] + f_bin;
            for (int j=0; j<freq.length; j++) {
                if (freq[j] >= lower_edge && freq[j] <= upper_edge) {
                    count++;
                    sum += val[j];
                }
            }
            new_val[i] = sum / count;
        }

        TTF_result results = new TTF_result();
        results.freq = new_freq;
        results.val = new_val;

        return results;
    }


    /**
     *
     * @param N
     * @return
     */
    private double[] hannWindow(int N) {
        // See https://en.wikipedia.org/wiki/Hann_function
        double[] hann_window = new double[N];

        for (int n=0; n<N; n++) {
            hann_window[n] = 0.5 * (1.0 - Math.cos((2.0 * Math.PI * n) / (N - 1)));
        }

        return hann_window;
    }


    /**
     *
     */
    static class TTF_result {
        double[] freq;
        double[] val;
    }


    /**
     *
     */
    static class ESF_result {
        double[] pos;
        double[] val;
    }


    /**
     *
     */
    static class CNR_result {
        double cnr;
        double contrast;
        double noise;
    }


    /** Method to obtain the index of element val in array arr
     *
     * @param arr
     * @param val
     * @param <T>
     * @return
     */
    public static <T> int indexOf(T[] arr, T val) {
        return Arrays.asList(arr).indexOf(val);
    }


    /**
     *
     * @param values Two-dimensional array of Float values
     * @return Two-dimensional array of Double values
     */
    private static double[][] convertFloatToDouble(final float[][] values) {
        int width = values[0].length;
        int height = values.length;

        final double[][] result = new double[width][height];

        for (int i=0; i<width; i++) {
            for (int j=0; j<width; j++) {
                result[i][j] = (double) values[i][j];
            }
        }
        return result;
    }


    /** Method to transpose a 2d array of complex numbers: the row and column indices are swapped.
     *   1   2   3   4                 1   5   9  13
     *   5   6   7   8    Transpose    2   6  10  14
     *   9  10  11  12  <----------->  3   7  11  15
     *  13  14  15  16                 4   8  12  16
     *
     * @param input A two-dimensional array of complex numbers
     * @return The transpose of the input, also a two-dimensional array of complex numbers
     */
    private Complex[][] transpose(Complex[][] input) {
        Complex[][] output = new Complex[input[0].length][input.length];

        for (int i=0; i<input.length; i++) {
            for (int j=0; j<input[0].length; j++) {
                output[j][i] = input[i][j];
            }
        }

        return output;
    }


    /**
     *
     * @param values Two-dimensional array of Double values
     * @return Two-dimensional array of Float values
     */
    private static float[][] convertDoubleToFloat(final double[][] values) {
        int width = values[0].length;
        int height = values.length;

        final float[][] result = new float[width][height];

        for (int i=0; i<width; i++) {
            for (int j=0; j<width; j++) {
                result[i][j] = (float) values[i][j];
            }
        }
        return result;
    }

    /**
     *
     * @param values Two-dimensional array of Double values
     * @return Two-dimensional array of Complex values
     */
    private static float[][] convertComplexToFloat(final Complex[][] values) {
        int width = values[0].length;
        int height = values.length;

        final float[][] result = new float[width][height];

        for (int i=0; i<width; i++) {
            for (int j=0; j<width; j++) {
                result[i][j] = (float) values[i][j].abs();
            }
        }
        return result;
    }

    /**
     *
     * @param ip The ImageProcessor from which to get the values from
     * @param roi The Roi to use to obtain the values - assumed to be rectangular
     * @return Two-dimensional array containing the calibrated pixel values of the Roi
     */
    private static double[][] getCalibratedPixelValues(ImageProcessor ip, Roi roi) {
        int width = (int) roi.getFloatWidth();
        int height = (int) roi.getFloatHeight();

        final double[][] result = new double[width][height];

        Point[] p = roi.getContainedPoints();

        int current_point = 0;
        for (int x=0; x<width; x++) {
            for (int y=0; y<height; y++) {
                result[x][y] = ip.getPixelValue(p[current_point].x, p[current_point].y);
                current_point++;
            }
        }

        return result;
    }
}
