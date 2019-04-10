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

    private boolean find_com = true;
    private double[] com = new double[2];
    private boolean use_low_pass = false;
    private double low_pass_factor = 0.9;

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
        return this.ttf_result;
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
        return this.cnr_result;
    }

    // Defined to enable retrieval of centre of mass
    double[] getCentreOfMass() {
        return this.com;
    }


    /** Method to get the current frequency binning increment in 1/mm
     *
     * @return The current frequency binning increment
     */
    double getFreqBinInc() {
        return rebin_f_inc;
    }


    /** Method to set up the TTF
     *
     * @param arg Not currently used
     * @param imp The ImagePlus object
     * @return The return value - 0 if unsuccessful
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

        GenericDialog gd = new GenericDialog("Task transfer function options");
        gd.addNumericField("Object diameter (mm) [Mercury 25, Catphan 12]", this.obj_dia_mm, 0);
        gd.addNumericField("Pixel reduction factor", this.pixel_reduction_factor, 0);
        gd.addNumericField("Frequency rebin increment (mm^-1)", this.rebin_f_inc, 2);
        gd.addCheckbox("Find centre of mass", this.find_com);
        gd.addCheckbox("Use low-pass filter for CoM", this.use_low_pass);
        gd.addNumericField("Low-pass filter factor", this.low_pass_factor, 2);
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
        this.use_low_pass = gd.getNextBoolean();
        this.low_pass_factor = gd.getNextNumber();
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
            this.imp.deleteRoi();
            this.imp = z_projector.getProjection();
            this.imp.setCalibration(this.cal);
            this.imp.show();
            this.imp.setRoi(this.main_roi);
        }

        // If using low pass filter for centre of mass then ROI must be a power of two
        if (this.use_low_pass) {
            double N = this.main_roi.getFloatWidth();
            int M = (int) (Math.log(N) / Math.log(2.0));
            int N_power_two = (int) Math.pow(2, M);
            if (N_power_two != N) { // Make the Roi larger using the next power of 2 bigger than the current width
                N_power_two = (int) Math.pow(2, M + 1);
                this.main_roi = new Roi(this.main_roi.getXBase() - (N_power_two - N) / 2.0, this.main_roi.getYBase() - (N_power_two - N) / 2.0, N_power_two, N_power_two);
                this.imp.setRoi(this.main_roi);
                this.imp.updateAndDraw();
            }
        }

        return DOES_ALL;
    }

    /**
     *
     * @param ip The ImageProcessor
     */
    public void run(ImageProcessor ip) {
        if (this.use_average_stack) {
            ip = this.imp.getProcessor();
        }


        //---------------------------------------------------------------------
        // Obtain the pixel size of the image. This will return 1.0 if ImageJ
        // isn't able to read the pixel size from the DICOM tags.
        double pixel_size_in_mm = this.cal.pixelWidth;
        double obj_dia_pix = this.obj_dia_mm/pixel_size_in_mm;
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Determine the centre of mass (Com). Initialise the CoM to the centre
        // of the main_roi.
        this.com[0] = this.cal.getX(this.main_roi.getXBase() + (this.main_roi.getFloatWidth()/2.0));
        this.com[1] = this.cal.getY(this.main_roi.getYBase() + (this.main_roi.getFloatHeight()/2.0));

        double inner_dia = obj_dia_pix * 1.5;//+ (10.0 / pixel_size_in_mm);
        double outer_dia = obj_dia_pix * 2.0;
        Roi inner_roi = new OvalRoi(this.cal.getRawX(this.com[0]) - inner_dia/2.0, this.cal.getRawY(this.com[1]) - inner_dia/2.0, inner_dia, inner_dia);
        Roi outer_roi = new OvalRoi(this.cal.getRawX(this.com[0]) - outer_dia/2.0, this.cal.getRawY(this.com[1]) - outer_dia/2.0, outer_dia, outer_dia);
        Roi annular_roi = new ShapeRoi(outer_roi).xor(new ShapeRoi(inner_roi));

        if (this.find_com) {
            //---------------------------------------------------------------------
            // Obtain the centre of mass of the ROI set by the user.
            this.imp.setRoi(this.main_roi);
            this.com = centreOfMass(ip, this.main_roi, this.cal, 0.0);
            if (this.log_results) IJ.log("Initial CoM:," + this.three_dp.format(this.com[0]) + "," + this.three_dp.format(this.com[1]));
            //---------------------------------------------------------------------


            //---------------------------------------------------------------------
            // Recalculate the CoM after subtracting the background using an
            // annular region of interest centred on the centre of mass that has an
            // internal diameter of object diameter + 10 mm and an external
            // diameter of 2 x object diameter.
            annular_roi.setLocation(this.cal.getRawX(this.com[0]) - annular_roi.getFloatWidth()/2.0, this.cal.getRawY(this.com[1]) - annular_roi.getFloatHeight()/2.0);

            // Obtain the mean pixel value of the annular ROI
            this.imp.setRoi(annular_roi);
            ImageStatistics annular_stats = this.imp.getAllStatistics();
            double bgd_mean = annular_stats.mean;

            // Recalculate the CoM with background subtracted
            outer_roi.setLocation(this.cal.getRawX(this.com[0]) - outer_roi.getFloatWidth()/2.0, this.cal.getRawY(this.com[1]) - outer_roi.getFloatHeight()/2.0);
            this.imp.setRoi(outer_roi);
            this.com = centreOfMass(ip, outer_roi, this.cal, -bgd_mean);
            if (this.log_results) IJ.log("Adjusted CoM:," + this.three_dp.format(this.com[0]) + "," + this.three_dp.format(this.com[1]));

            // Recentre the annualar roi on the new CoM and recalculate CoM  for a final time
            //annular_roi.setLocation(this.cal.getRawX(com[0]) - obj_dia_pix, this.cal.getRawY(com[1]) - obj_dia_pix);
            annular_roi.setLocation(this.cal.getRawX(this.com[0]) - annular_roi.getFloatWidth()/2.0, this.cal.getRawY(this.com[1]) - annular_roi.getFloatHeight()/2.0);
            this.imp.setRoi(annular_roi);
            annular_stats = this.imp.getAllStatistics();
            bgd_mean = annular_stats.mean;
            this.imp.setRoi(outer_roi);
            this.com = centreOfMass(ip, outer_roi, cal, -bgd_mean);

            // Use the low-pass filter method if selected
            if (this.use_low_pass) {
                this.main_roi.setLocation(this.cal.getRawX(this.com[0]) - this.main_roi.getFloatWidth()/2.0, this.cal.getRawY(this.com[1]) - this.main_roi.getFloatHeight()/2.0);
                double[] com_low_pass = centreOfMassLowFreqFilter(ip, this.main_roi, this.cal, this.obj_dia_mm, -bgd_mean);
                if (com_low_pass != null) {
                    this.com = com_low_pass;
                    if (this.log_results)
                        IJ.log("CoM after low-pass filter:," + this.three_dp.format(com[0]) + "," + this.three_dp.format(com[1]));
                } else {
                    IJ.error("Could not work out centre of mass using low-pass filter.\nUsing existing centre of mass.");
                }
            }
        }

        if (this.log_results) IJ.log("Final CoM:," + this.three_dp.format(this.com[0]) + "," + this.three_dp.format(this.com[1]));

        // Add a marker to the image showing the location of the CoM
        Roi com_marker = new PointRoi(this.cal.getRawX(this.com[0]), this.cal.getRawY(this.com[1]));
        Overlay com_overlay = new Overlay(com_marker);
        com_overlay.setStrokeColor(Color.red);
        this.imp.setOverlay(com_overlay);

        // Recentre the outer_roi and annular_roi using the centre of mass.
        outer_roi.setLocation(this.cal.getRawX(this.com[0]) - outer_roi.getFloatWidth()/2.0, this.cal.getRawY(this.com[1]) - outer_roi.getFloatHeight()/2.0);
        annular_roi.setLocation(this.cal.getRawX(this.com[0]) - annular_roi.getFloatWidth()/2.0, this.cal.getRawY(this.com[1]) - annular_roi.getFloatHeight()/2.0);

        // Create an overlay to display the annular ROI and CoM on the image.
        Overlay overlay = new Overlay();
        this.imp.setOverlay(overlay);
        overlay.add(this.main_roi);
        overlay.add(annular_roi);
        overlay.add(com_marker);
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Calculate the CNR of the object
        double obj_roi_dia = obj_dia_pix * 0.5; // - (3.0 / pixel_size_in_mm);
        Roi obj_roi = new OvalRoi(this.cal.getRawX(this.com[0])-obj_roi_dia/2.0, this.cal.getRawY(this.com[1])-obj_roi_dia/2.0, obj_roi_dia, obj_roi_dia);
        CNR_result cnr_results = contrastToNoiseRatio(this.imp, obj_roi, annular_roi);

        if (this.log_results) IJ.log("CNR is: " + this.three_dp.format(cnr_results.cnr));
        if (this.log_results) IJ.log("Contrast is: " + this.three_dp.format(cnr_results.contrast));
        if (this.log_results) IJ.log("Noise is: " + this.three_dp.format(cnr_results.noise));
        if (this.log_results) IJ.log("Background is: " + this.three_dp.format(cnr_results.background));
        if (this.log_results) IJ.log("Signal is: " + this.three_dp.format(cnr_results.signal));

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
            raw_esf_pos[i] = calculateDistanceBetweenPoints(this.com[0], this.com[1], this.cal.getX(p.x + 0.5), this.cal.getY(p.y + 0.5));  // The value is associated with the centre of the pixel, so add 0.5
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
            double[] smoothed_esf;

            try {
                smoothed_esf = new LoessInterpolator(this.loess_bandwidth, this.loess_robustness).smooth(raw_esf_pos, raw_esf_val);
                function = interpolator.interpolate(raw_esf_pos, smoothed_esf);
            }
            catch (Exception e) {
                IJ.error("There was an error when interpolating the raw ESF:\n" + e);
                return;
            }

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
        FiniteDifferencesDifferentiator differentiator = new FiniteDifferencesDifferentiator(2, rebinned_sample_inc);
        try {
            function = interpolator.interpolate(esf_rebinned_pos, esf_rebinned_val);
        }
        catch (Exception e) {
            IJ.error("There was an error when interpolating the rebinned ESF:\n" + e);
            return;
        }
        UnivariateDifferentiableFunction complete_f = differentiator.differentiate(function);
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


    /** Method to calculate the distance between two points in cartesian space
     *
     * @param x1 The x coordinate of point 1
     * @param y1 The y coordinate of point 1
     * @param x2 The x coordinate of point 2
     * @param y2 The y coordinate of point 2
     * @return The distance between points 1 and 2
     */
    private double calculateDistanceBetweenPoints(double x1, double y1, double x2,double y2) {
        // Calculate the distance between two points in cartesian space.
        return Math.sqrt( (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }


    /** Method to return the centre of mass of the supplied Roi on the ImageProcessor
     *
     * @param ip The ImageProcessor to use
     * @param roi The Roi to place on the ImageProcessor
     * @param cal The calibration data for the ImageProcessor
     * @param pv_offset The offset to add to each pixel value before using it in the calculations
     * @return The centre of mass as a two-element array
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
            x_wt += cal.getX(p.x + 0.5) * pixel_val; // The value is associated with the centre of the pixel, so add 0.5
            y_wt += cal.getY(p.y + 0.5) * pixel_val; // The value is associated with the centre of the pixel, so add 0.5
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
        Complex[][] fft;
        try {
            fft = twod_fourier_transformer.transform(roi_pixels, DftNormalization.STANDARD);
        }
        catch (Exception e) {
            IJ.error("There was an error when Fourier transforming the region:\n" + e);
            return(null);
        }

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
        double high_freq_cut_off = this.low_pass_factor / obj_diam;//1.0 / obj_diam;
        //double low_freq_cut_off = this.high_pass_factor / obj_diam;

        Complex complexZero = new Complex(0.0, 0.0);

        // First the top-left quadrant, +ve i and +ve j
        for(int i=0; i<roi_pixels.length/2.0; i++) {
            for(int j=0; j<roi_pixels[0].length/2.0; j++) {
                // Calculate distance from the top-left hand corner (0,0)
                double dist = calculateDistanceBetweenPoints(0,0, i, j);
                double radial_freq = dist * freq_inc;
                if (radial_freq > high_freq_cut_off) fft[i][j] = complexZero;
                //if (radial_freq > high_freq_cut_off || radial_freq < low_freq_cut_off) fft[i][j] = complexZero;
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
                //if (radial_freq > high_freq_cut_off || radial_freq < low_freq_cut_off) fft[i][j] = complexZero;
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
                //if (radial_freq > high_freq_cut_off || radial_freq < low_freq_cut_off) fft[i][j] = complexZero;
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
                //if (radial_freq > high_freq_cut_off || radial_freq < low_freq_cut_off) fft[i][j] = complexZero;
            }
        }


        // Now shift the quadrants so that zero is at the centre
        fft = twod_fourier_transformer.swapQuadrants(fft);

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
        fft = twod_fourier_transformer.transform(fft, DftNormalization.STANDARD, TransformType.INVERSE);


        // ---------------------------------------------------------------------
        ImagePlus temp = NewImage.createFloatImage("Inverse FFT", fft.length, fft[0].length, 1, NewImage.FILL_BLACK);
        temp.getProcessor().setFloatArray(convertComplexToFloat(fft));
        temp.show();
        temp.updateAndDraw();
        IJ.run(temp, "Enhance Contrast", "saturated=0.35");
        IJ.run("Set... ", "zoom=800");
        // ---------------------------------------------------------------------


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
                x_wt += (x+roi_x_base+0.5) * pixel_val; // The value is associated with the centre of the pixel, so add 0.5
                y_wt += (y+roi_y_base+0.5) * pixel_val; // The value is associated with the centre of the pixel, so add 0.5
                total_wt += pixel_val;
            }
        }
        double x_com = cal.getX(x_wt / total_wt);
        double y_com = cal.getY(y_wt / total_wt);

        return new double[] {x_com, y_com};
    }


    /** Method to calculate the contrast to noise ratio
     *
     * @param imp The ImagePlus to use
     * @param obj_roi The Roi of the object (the "signal")
     * @param bgd_roi The Roi of the background (the "background")
     * @return An object containing the CNR, signal value, background value, contrast and noise
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
        result.signal = obj_stats.mean;
        result.background = bgd_stats.mean;
        result.contrast = result.signal - result.background;
        result.noise = bgd_stats.stdDev;
        result.cnr = Math.abs(result.contrast) / result.noise;

        return result;
    }


    /** Method to determine the task transfer function from the line spread function
     *
     * @param lsf The line spread function values
     * @param lsf_inc_mm The spacing of the supplied line spread function in mm
     * @param px_size_mm The pixel size of the image in mm
     * @return The calculated TTF
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

        // Plot the Hann-windowed ESF
        Plot plot_lsf_hann = new Plot("Line spread function (Hann applied)", "Uncalibrated", "Pixel value");
        plot_lsf_hann.add("line", lsf);
        plot_lsf_hann.show();

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


    /** Method to differentiate the supplied values
     *
     * @param locations An array of locations where the differential is required
     * @param diff_func The differentiable function
     * @return The differentials at each of the supplied locations
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


    /** Method to resample the TTF at a specified frequency increment using linear interpolation
     *
     * @param ttf The starting TTF
     * @param freq_inc The required frequency increment
     * @return The resampled TTF
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


    /** Method to log the TTF results to the ImageJ log
     *
     * @param ttf The TTF results
     * @param ttf_label The label to use for the TTF value column
     */
    private void logTTF(TTF_result ttf, String ttf_label) {
        IJ.log("Freq (mm^-1)," + ttf_label);

        for (int i=0; i<ttf.freq.length; i++) {
            IJ.log("" + this.two_dp.format(ttf.freq[i]) + ", " + this.three_dp.format(ttf.val[i]));
        }
        IJ.log("\n");
    }


    /** Method to calculate the frequency at a specified TTF value
     *
     * @param ttf The TTF
     * @param req_ttf_val The TTF value for which the frequency is required
     * @return The frequency of the TTF at the supplied TTF value
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

        try {
            return fn_ttf_lin.value(req_ttf_val);
        }
        catch (Exception e) {
            IJ.error("There was an error when obtaining the TTF frequency at " + req_ttf_val + ":\n" + e);
            return 0.0;
        }
    }


    /** Method to rebin the ESF to a regularly-spaced array
     *
     * @param pos The raw positions of the ESF
     * @param val The value of the ESF at each raw position
     * @param pos_inc The required position increment of the rebinned ESF
     * @param pos_bin The bin width to use when calculating the ESF value at each position
     * @return The rebinned ESF
     */
    private ESF_result rebinESF(double[] pos, double[] val, double pos_inc, double pos_bin) {

        double[][] rebinned_data = rebin(pos, val, pos_inc, pos_bin, true);
        ESF_result results = new ESF_result();
        results.pos = rebinned_data[0];
        results.val = rebinned_data[1];

        return results;
    }


    /** Method to rebin the TTF
     *
     * @param freq The frequencies of the starting TTF
     * @param val The TTF value at each starting frequency
     * @param f_inc The frequency increment for the rebinned TTF
     * @param f_bin The bin width to use when calculating the TTF value at each frequency
     * @return The rebinned TTF
     */
    private TTF_result rebinTTF(double[] freq, double[] val, double f_inc, double f_bin) {

        double[][] rebinned_data = rebin(freq, val, f_inc, f_bin, false);
        TTF_result results = new TTF_result();
        results.freq = rebinned_data[0];
        results.val = rebinned_data[1];

        return results;
    }


    /** Method to rebin a series of x, y values to regular x spacing of x_inc
     *
     * @param x The original x values
     * @param y The value of y at each original x
     * @param x_inc The required x increment
     * @param x_bin The width of the bin to use to calculate the new value of y at each new value of x
     * @param calc_min_x Use min(x) or 0 as the minimum x value
     * @return The rebinned data
     */
    private double[][] rebin(double[] x, double[] y, double x_inc, double x_bin, boolean calc_min_x) {
        double max_x = StatUtils.max(x);
        double min_x = 0.0;
        if (calc_min_x) min_x = StatUtils.min(x);
        int num_intervals = (int) Math.round((max_x - min_x) / x_inc) + 1;

        double [][] rebinned = new double[2][num_intervals];

        // Calculate the new x values
        for (int i = 0; i < num_intervals; i++) {
            rebinned[0][i] = min_x + (i * x_inc);
        }

        // Calculate the new y value at each new x position using the mean of all original y values that fall
        // within x+x_bin >= x >= x-x_bin
        for (int i = 0; i < num_intervals; i++) {
            double sum = 0;
            double count = 0;
            double lower_edge = rebinned[0][i] - x_bin;
            double upper_edge = rebinned[0][i] + x_bin;
            for (int j = 0; j < x.length; j++) {
                if (x[j] >= lower_edge && x[j] <= upper_edge) {
                    count++;
                    sum += y[j];
                }
            }
            rebinned[1][i] = sum / count;
        }

        return rebinned;
    }


    /** Method to calculate a Hann window
     *
     * @param N The number of data points in the window
     * @return The values of the Hann window
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
        double signal;
        double background;
    }


    /** Method to obtain the index of element val in array arr
     *
     * @param arr The array being searched
     * @param val The value being looked for
     * @param <T> The data type of the array and value being looked for
     * @return The array index which contains the supplied value
     */
    private static <T> int indexOf(T[] arr, T val) {
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
