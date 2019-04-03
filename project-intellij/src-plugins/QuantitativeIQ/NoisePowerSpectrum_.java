package QuantitativeIQ;
// ImageJ plugin to calculate the noise power spectrum
//
// The routine will run on a user-selected oval region of interest (ROI). If no
// ROI is selected then the routine uses the whole image.
//
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.*;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.util.MathArrays;

import java.awt.*;
import java.text.DecimalFormat;

/**
 *
 */
public class NoisePowerSpectrum_ implements PlugInFilter {
    private ImagePlus imp;
    private Roi main_roi;

    private int num_records = 16;
    private int record_size = 64;
    private double record_radius = 100.0; // in mm

    private double rebin_f_inc = 0.05; // in mm^-1

    private static NPS_result nps_result;

    private boolean log_results = true;

    private boolean use_average_stack = false;
    private int min_stack_index;
    private int max_stack_index;

    private DecimalFormat two_dp = new DecimalFormat("0.00");
    private DecimalFormat three_dp = new DecimalFormat("0.000");


    /** Method to enable or disable writing to the ImageJ log
     *
     * @param log_contition Boolean value to set logging on or off
     */
    void setLogging(boolean log_contition) {
        this.log_results = log_contition;
    }


    /** Method to obtain the NPS results
     *
     * @return An object of type NPS_result containing the NPS results. The .freq property contains
     * an array of frequency values; the .val property contains the NPS values at the corresponding
     * frequency.
     */
    NPS_result getNPSResult() {
        return nps_result;
    }


    /**
     *
     * @param arg Setup argument; not in current use.
     * @param imp The image
     * @return 0 if successful; 1 if not
     */
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;

        if(this.imp == null) {
            return(0);
        }

        // Remove any overlays that are present
        imp.setOverlay(null);

        // If there is no ROI set then use the whole image
        this.main_roi = imp.getRoi();
        if(this.main_roi == null) {
            imp.setRoi(0, 0, imp.getWidth(), imp.getHeight());
            this.main_roi = imp.getRoi();
        }

        if ((arg != null) && (arg.length() > 0)) {
            if (arg.equals("defaults")) {
                // Use defaults
                if (this.log_results) IJ.log("NoisePowerSpectrum using ROIs " + this.record_radius + " mm from centre");
            }
            else {
                // Assume the argument contains the distance from ROI centre to record in mm.
                this.record_radius = Float.parseFloat(arg);
                if (this.log_results) IJ.log("NoisePowerSpectrum using ROIs " + this.record_radius + " mm from centre");
            }
        }
        else {
            GenericDialog gd = new GenericDialog("Noise power options");
            gd.addNumericField("Distance from centre to records (mm)", this.record_radius, 0);
            gd.addNumericField("Number of records", this.num_records, 0);
            gd.addNumericField("Record size (pixels), must be a power of 2", this.record_size, 0);
            gd.addNumericField("Frequency rebin increment (mm^-1)", this.rebin_f_inc, 2);
            if (this.imp.isStack()) {
                this.use_average_stack = true;
                gd.addCheckbox("Run on the average stack", this.use_average_stack);
                gd.addNumericField("Min slice index", 1, 0);
                gd.addNumericField("to max slice index", this.imp.getImageStackSize(), 0);
            }
            gd.showDialog();
            if(gd.wasCanceled()) {
                return(1);
            }

            this.record_radius = gd.getNextNumber();
            this.num_records = (int) gd.getNextNumber();
            this.record_size = (int) gd.getNextNumber();
            this.rebin_f_inc = gd.getNextNumber();
            if (this.imp.isStack()) {
                this.use_average_stack = gd.getNextBoolean();
                this.min_stack_index = (int) gd.getNextNumber();
                this.max_stack_index = (int) gd.getNextNumber();
            }
        }
        return DOES_ALL;
    }

    /**
     *
     * @param ip The ImageProcessor that is passed to this method
     */
    public void run(ImageProcessor ip) {
        //---------------------------------------------------------------------
        // Obtain the pixel size of the image. This will return 1.0 if ImageJ
        // isn't able to read the pixel size from the DICOM tags.
        Calibration cal = this.imp.getCalibration();
        double pixel_size_in_mm = cal.pixelWidth;
        double record_radius_px = this.record_radius/pixel_size_in_mm;
        //---------------------------------------------------------------------


        double[][] ensemble_nps = new double[this.record_size][this.record_size];


        //---------------------------------------------------------------------
        // Determine the centre of the main_roi and use this as the centre of
        // the region for which we're measuring the NPS.
        double[] centre = new double[2];
        centre[0] = cal.getX(this.main_roi.getXBase() + (this.main_roi.getFloatWidth()/2.0));
        centre[1] = cal.getY(this.main_roi.getYBase() + (this.main_roi.getFloatHeight()/2.0));

        // Add a marker to the image showing the location of the centre
        Roi com_marker = new PointRoi(cal.getRawX(centre[0]), cal.getRawY(centre[1]));
        Overlay com_overlay = new Overlay(com_marker);
        com_overlay.setStrokeColor(Color.red);
        this.imp.setOverlay(com_overlay);

        // Create an overlay to display the record ROIs on the image.
        Overlay overlay = new Overlay();
        this.imp.setOverlay(overlay);
        overlay.add(com_marker);
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Calculate the NPS ensemble
        if (this.use_average_stack) {
            if (this.imp.isStack()) {
                ImageStack stack = this.imp.getStack();
                for (int i=this.min_stack_index; i<=this.max_stack_index; i++) {
                    ImageProcessor im_proc = stack.getProcessor(i);
                    ensemble_nps = add2dArrays(ensemble_nps, ensembleNPS(im_proc, record_radius_px, cal, centre, overlay));
                }
            }
            else {
                ImageProcessor im_proc = this.imp.getProcessor();
                ensemble_nps = ensembleNPS(im_proc, record_radius_px, cal, centre, overlay);
            }
        }
        else {
            ImageProcessor im_proc = this.imp.getProcessor();
            ensemble_nps = ensembleNPS(im_proc, record_radius_px, cal, centre, overlay);
        }
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Multiply the ensemble by dx, dy and divide by M, Nx and Ny
        double multiplication_factor = pixel_size_in_mm * pixel_size_in_mm / this.num_records / this.record_size / this.record_size;

        // Include provision for the number of images in the stack, if this was used
        if (this.use_average_stack) multiplication_factor /= this.imp.getImageStackSize();

        ensemble_nps = multiply2dDoubleArray(ensemble_nps, multiplication_factor);
        // ---------------------------------------------------------------------


        // ---------------------------------------------------------------------
        // Display the resulting 2d NPS
        ImagePlus nps_2d = NewImage.createFloatImage("2d NPS", this.record_size, this.record_size, 1, NewImage.FILL_BLACK);
        nps_2d.getProcessor().setFloatArray(convertDoubleToFloat(ensemble_nps));
        IJ.run(nps_2d, "Swap Quadrants", "");
        nps_2d.show();
        nps_2d.updateAndDraw();
        IJ.run(nps_2d, "Enhance Contrast", "saturated=0.35");
        IJ.run("Set... ", "zoom=800");
        // ---------------------------------------------------------------------


        // ---------------------------------------------------------------------
        // Work out the radial NPS
        //
        // The frequency increment between values in the i or j direction is calculated as:
        //     1/(record_size * pixel_size_in_mm)
        //
        // The top-left quadrant of the nps ensemble contains frequencies:
        //     0 ≤ i ≤ i_max
        //     0 ≤ j ≤ j_max
        //
        // The top-right:
        //     -i_max ≤ i < 0
        //     0 ≤ j ≤ j_max
        //
        // The bottom-left:
        //     0 ≤ i ≤ i_max
        //     -j_max ≥ j > 0
        //
        // The bottom-right:
        //     -i_max ≤ i < 0
        //     -j_max ≥ j > 0
        //
        // In our use-case we only need to use the top half of the FFT as there is symmetry
        // between the top and bottom halves because our image data has a zero imaginary
        // component.
        //
        double[] radial_nps = new double[(int)(this.record_size * this.record_size / 2.0)];

        double[] radial_freq = new double[(int)(this.record_size * this.record_size / 2.0)];
        double freq_inc = 1.0 / (this.record_size * pixel_size_in_mm);
        int current_point = 0;

        // First the +i and +j quadrant
        for(int i=0; i<this.record_size/2.0; i++) {
            for(int j=0; j<this.record_size/2.0; j++) {
                radial_nps[current_point] = ensemble_nps[i][j];
                double dist = calculateDistanceBetweenPoints(0,0, i, j);
                radial_freq[current_point] = dist * freq_inc;
                current_point++;
            }
        }

        // Next the -i and +j quadrant
        int i_start = (int) (this.record_size/2.0);
        for(int i=i_start; i<this.record_size; i++) {
            for(int j=0; j<this.record_size/2.0; j++) {
                radial_nps[current_point] = ensemble_nps[i][j];
                double dist = calculateDistanceBetweenPoints(this.record_size, 0, i, j);
                radial_freq[current_point] = dist * freq_inc;
                current_point++;
            }
        }

        // Sort the radial array data into ascending order of frequency
        MathArrays.sortInPlace(radial_freq, radial_nps);

        // Plot the radial frequency
        Plot plot_nps = new Plot("Noise power spectrum", "Frequency (mm^-1)", "Noise power (HU^2.mm^2)");
        plot_nps.setAxisYLog(true);
        plot_nps.add("dot", radial_freq, radial_nps);
        plot_nps.show();

        double IEC_f_bin = 0.01 / pixel_size_in_mm;
        NPS_result rebinned_NPS = rebinNPS(radial_freq, radial_nps, rebin_f_inc, IEC_f_bin);

        plot_nps.add("line", rebinned_NPS.freq, rebinned_NPS.val);
        // ---------------------------------------------------------------------


        // Write the results to the ImageJ log
        if (this.log_results) logNPS(rebinned_NPS, "NPS (HU^2.mm^2)");

        nps_result = rebinned_NPS;
    }


    private double[][] ensembleNPS(ImageProcessor im_proc, double radius_px, Calibration cal, double[] centre, Overlay overlay) {

        ImagePlus imp = new ImagePlus();
        imp.setProcessor(im_proc);
        imp.setCalibration(cal);

        double[][] ensemble_nps = new double[this.record_size][this.record_size];

        for (int m = 0; m < this.num_records; m++) {
            double current_angle = (2.0 * Math.PI) / this.num_records * m;

            double x_pos = radius_px * Math.cos(current_angle);
            double y_pos = radius_px * Math.sin(current_angle);

            Roi record_roi = new Roi(cal.getRawX(centre[0]) + x_pos - (this.record_size / 2.0), cal.getRawY(centre[1]) + y_pos - (this.record_size / 2.0), this.record_size, this.record_size);

            overlay.add(record_roi);
            overlay.setStrokeColor(Color.green);

            // Set the roi on the main image
            imp.setRoi(record_roi);

            // Carry out a 2d second-order polynomial fit on a new ImagePlus
            // that just contains the current record pixels. The Fit Polynomial
            // command subtracts the fit from the original data.
            ImagePlus record = imp.crop();

            IJ.run(record, "Fit Polynomial", "x=2 y=2 mixed=0");
            record.setRoi(0, 0, record.getWidth(), record.getHeight());
            Roi whole_record_roi = record.getRoi();
            double[][] record_vals = getCalibratedPixelValues(record.getProcessor(), whole_record_roi);

            // Fourier transform the result
            TwoDFFT twod_transformer = new TwoDFFT();
            Complex[][] nps = twod_transformer.transform(record_vals, DftNormalization.STANDARD);

            // Calculate the magnitude and square it
            double[][] nps_abs_sqrd = absSqrdOfComplex2dArray(nps);

            // Add the current nps_abs_squared to the ensemble
            ensemble_nps = add2dArrays(ensemble_nps, nps_abs_sqrd);
        }
        return ensemble_nps;
    }


    /**
     *
     * @param values Two-dimensional array whose elements are to be changed
     * @param mult_val The factor to change the array values by
     * @return Two-dimensional array of the same dimensions as the input
     */
    private static double[][] multiply2dDoubleArray(final double[][] values, final double mult_val) {
        int width = values[0].length;
        int height = values.length;

        final double[][] result = new double[width][height];

        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                result[i][j] = values[i][j] * mult_val;
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
    private static Complex[][] convertDoubleToComplex(final double[][] values) {
        int width = values[0].length;
        int height = values.length;

        final Complex[][] result = new Complex[width][height];

        for (int i=0; i<width; i++) {
            result[i] = ComplexUtils.convertToComplex(values[i]);
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


    /**
     *
     * @param values Two-dimensional Complex array to be worked ou
     * @return Two-dimensional Double array with the resulting values
     */
    private static double[][] absSqrdOfComplex2dArray(final Complex[][] values) {
        int width = values[0].length;
        int height = values.length;

        final double[][] result = new double[width][height];

        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                result[i][j] = values[i][j].abs() * values[i][j].abs();
            }
        }
        return result;
    }


    /**
     *
     * @param arr1 Two-dimensional array
     * @param arr2 Two-dimensional array to be added to arr1. Must be same dimensions as arr1.
     * @return Two-dimensional array which is the sum of corresponding elements in arr1 + arr2.
     */
    private static double[][] add2dArrays(final double[][] arr1, final double[][] arr2) {
        int width = arr1[0].length;
        int height = arr1.length;

        final double[][] result = new double[width][height];

        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                result[i][j] = arr1[i][j] + arr2[i][j];
            }
        }
        return result;
    }


    /**
     *
     * @param x1 X location of point 1
     * @param y1 Y location of point 2
     * @param x2 X location of point 1
     * @param y2 y location of point 2
     * @return The straight-line distance between point 1 and 2
     */
    private double calculateDistanceBetweenPoints(double x1, double y1, double x2,double y2) {
        // Calculate the distance between two points in cartesian space.
        return Math.sqrt( (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }


    /**
     *
     * @param freq The raw frequency data
     * @param val The raw value data
     * @param f_inc The frequency increments of the rebinned NPS
     * @param f_bin The width of the frequency bin to use to find values to average at each binned frequency
     * @return The rebinned array results
     */
    private NPS_result rebinNPS(double[] freq, double[] val, double f_inc, double f_bin) {

        double max_f = StatUtils.max(freq);
        int num_intervals = (int) Math.round(max_f / f_inc) + 1;

        double[] new_vals = new double[num_intervals];
        double[] new_f_scale = new double[num_intervals];
        for (int i=0; i<num_intervals; i++) {
            new_f_scale[i] = i * f_inc;
        }

        for (int i=0; i<num_intervals; i++) {
            double sum = 0;
            double count = 0;
            double lower_edge = new_f_scale[i] - f_bin;
            double upper_edge = new_f_scale[i] + f_bin;
            for (int j=0; j<freq.length; j++) {
                if (freq[j] >= lower_edge && freq[j] <= upper_edge) {
                    count++;
                    sum += val[j];
                }
            }
            new_vals[i] = sum / count;
        }

        NPS_result results = new NPS_result();
        results.freq = new_f_scale;
        results.val = new_vals;

        return results;
    }


    /**
     *
     * @param nps The NPS data to be written to the log
     * @param nps_label The label to use for the NPS value column
     */
    private void logNPS(NPS_result nps, String nps_label) {
        IJ.log("Freq (mm^-1)," + nps_label);

        for (int i=0; i<nps.freq.length; i++) {
            IJ.log("" + this.two_dp.format(nps.freq[i]) + ", " + this.three_dp.format(nps.val[i]));
        }
    }


    /**
     *
     */
    static class NPS_result {
        double[] freq;
        double[] val;
    }
}
