

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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.*;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.util.DicomTools;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.util.MathArrays;

import java.awt.*;
import java.text.DecimalFormat;


/**
 * An ImageJ plugin to calculate the noise power spectrum of a CT image or stack of images.
 */
public class NoisePowerSpectrum_ implements PlugInFilter {
    private ImagePlus imp;
    private Roi main_roi;

    private int num_records = 16;
    private int record_size = 64;
    private double record_radius = 100.0; // in mm
    private double record_stdev_sum = 0.0;

    private double rebin_f_inc = 0.05; // in mm^-1

    private static NPS_result nps_result;

    private boolean log_results = true;
    private boolean show_plots = true;

    private boolean use_substack = false;
    private int min_stack_index;
    private int max_stack_index;

    private DecimalFormat two_dp = new DecimalFormat("0.00");
    private DecimalFormat three_dp = new DecimalFormat("0.000");

    private String recon_filter = "";

    /** Method to enable or disable the display of plots
     *
     * @param show_plots Boolean value to set chart display on or off
     */
    void setPlotting(boolean show_plots) { this.show_plots = show_plots; }


    /** Method to enable or disable writing to the ImageJ log
     *
     * @param log_contition Boolean value to set logging on or off
     */
    void setLogging(boolean log_contition) {
        this.log_results = log_contition;
    }


    /** Method to return the substack stack status
     *
     * @return The substack setting
     */
    boolean getSubstackStatus() { return this.use_substack; }


    double getRecordStdevSum() { return this.record_stdev_sum; }
    int getNumRecords() { return this.num_records; }


    /** Method to return the minimum stack index that was used
     *
     * @return The minimum stack index
     */
    int getMinStackIndex() { return this.min_stack_index; }


    /** Method to return the maximum stack index that was used
     *
     * @return The maximum stack index
     */
    int getMaxStackIndex() { return this.max_stack_index; }


    /** Method to return the reconstruction filter
     *
     * @return The reconstruction filter
     */
    String getReconFilter() { return this.recon_filter; }


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
     * This method is automatically run by ImageJ when the plugin is called.
     * If there is no ROI set then the whole image is selected as an ROI.
     *
     * @param arg If <i>default</i> then analysis is run with default options. If a value then this is used as the
     *            radius of the ROIs. If empty then the user is prompted for the calculation options.
     * @param imp Passed to the routine automatically by ImageJ.
     * @return DONE if unsuccessful; DOES_ALL if complete
     */
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;

        if(this.imp == null) {
            return(DONE);
        }

        this.recon_filter = DicomTools.getTag(this.imp, "0018,1210");

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
                this.use_substack = true;
                gd.addCheckbox("Run on all images in a substack", this.use_substack);
                gd.addNumericField("Min slice index", 1, 0);
                gd.addNumericField("to max slice index", this.imp.getImageStackSize(), 0);
            }
            gd.showDialog();
            if(gd.wasCanceled()) {
                return(DONE);
            }

            this.record_radius = gd.getNextNumber();
            this.num_records = (int) gd.getNextNumber();
            this.record_size = (int) gd.getNextNumber();
            this.rebin_f_inc = gd.getNextNumber();
            if (this.imp.isStack()) {
                this.use_substack = gd.getNextBoolean();
                this.min_stack_index = (int) gd.getNextNumber();
                this.max_stack_index = (int) gd.getNextNumber();
            }
        }
        return DOES_ALL;
    }


    /**
     * This method is called automatically by ImageJ once the <i>setup</i> method has completed.
     * The centre of the main ROI set in the <i>setup</i> method is used as a reference location.
     * The <i>ensembleNPS</i> method is run on either the stack or single image depending on the
     * options the user has chosen. The 2D ensemble NPS is multiplied by:
     * <br><br>
     * d<sub>x</sub> * d<sub>y</sub> / (M * N<sub>x</sub> * N<sub>y</sub>)
     * <br><br>
     * where d<sub>x</sub> and d<sub>y</sub> are the pixel dimensions in the x and y directions (mm),
     * M is the number if records and N<sub>x</sub> and N<sub>y</sub> the x, y pixel dimensions of the record.
     * <br><br>
     * The ensemble NPS is displayed and a one-dimensional radial NPS calculated, plotted and written to the
     * ImageJ log.
     *
     * @param ip The ImageProcessor that is passed to this method automatically by ImageJ
     */
    public void run(ImageProcessor ip) {
        // Used to require at least ImageJ version 1.52b as this was when the isStack()
        // method of the ImagePlus class was introduced:
        // https://imagej.nih.gov/ij/developer/api/ij/ImagePlus.html#isStack--
        //
        // Now requires minimum of 1.52s as this was when a bug that prevented iteration
        // through all points in a rectangular ROI was fixed:
        // https://github.com/imagej/imagej-legacy/issues/222
        if (IJ.versionLessThan("1.52s")) return;


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
        if (this.use_substack) {
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
        this.imp.updateAndDraw(); // Update to show the ring of ROIs
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Multiply the ensemble by dx, dy and divide by M, Nx and Ny
        double multiplication_factor = pixel_size_in_mm * pixel_size_in_mm / this.num_records / this.record_size / this.record_size;

        // Include provision for the number of images the user chose to ue in the stack, if substack was selected
        if (this.use_substack) multiplication_factor /= (this.max_stack_index - this.min_stack_index + 1);

        ensemble_nps = multiply2dDoubleArray(ensemble_nps, multiplication_factor);
        // ---------------------------------------------------------------------



        // ---------------------------------------------------------------------
        if (this.show_plots) {
            // Display the resulting 2d NPS
            ImagePlus nps_2d = NewImage.createFloatImage("2d NPS", this.record_size, this.record_size, 1, NewImage.FILL_BLACK);
            nps_2d.getProcessor().setFloatArray(ArrayConversions.convertDoubleToFloat(ensemble_nps));
            IJ.run(nps_2d, "Swap Quadrants", "");
            nps_2d.show();
            nps_2d.updateAndDraw();
            IJ.run(nps_2d, "Enhance Contrast", "saturated=0.35");
            IJ.run("Set... ", "zoom=800");
        }
        // ---------------------------------------------------------------------


        // ---------------------------------------------------------------------
        // Calculate the radial NPS.
        // The frequency increment between values in the i or j direction is calculated as:
        //     1/(record_size * pixel_size_in_mm)
        //
        // By default the 2d ensemble_nps is arranged in a slightly non-obvious way:
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
        // The quadrants can be swapped so that zero frequency is near the centre.
        // Then loop through every location in the ensemble_nps:
        //      Write the frequency at this location to the nth position in a 1d freq array
        //      Write the value at this location to the nth position in a 1d noise array
        // Then sort the two arrays by ascending frequency and finally rebin into a regular frequency.

        // Swap the quadrants of the ensemble_nps to put zero frequency at the centre
        TwoDFFT twod_transformer = new TwoDFFT();
        ensemble_nps = twod_transformer.swapQuadrants(ensemble_nps);

        // Loop through every element in ensemble_nps
        int zero_u = record_size/2; // u location of zero frequency element
        int zero_v = record_size/2; // v location of zero frequency element
        double[] radial_freq = new double[record_size*record_size];
        double[] radial_noise = new double[record_size*record_size];
        double freq_inc = 1.0 / (this.record_size * pixel_size_in_mm);
        for (int u=0; u<ensemble_nps.length; u++) {
            for (int v=0; v<ensemble_nps[0].length; v++) {
                radial_noise[(u*record_size)+v] = ensemble_nps[u][v];
                radial_freq[(u*record_size)+v] = freq_inc * calculateDistanceBetweenPoints(zero_u, zero_v, u, v);
            }
        }
        // Sort the frequency array in ascending order; sort the values in the same order at the same time
        MathArrays.sortInPlace(radial_freq, radial_noise);

        // Rebin the results
        double IEC_f_bin = 0.01 / pixel_size_in_mm;
        NPS_result rebinned_NPS = rebinNPS(radial_freq, radial_noise, rebin_f_inc, IEC_f_bin);


        if (this.show_plots) {
            // Plot the radial frequency
            Plot plot_nps = new Plot("Noise power spectrum", "Frequency (mm^-1)", "Noise power (HU^2.mm^2)");
            plot_nps.setAxisYLog(true);
            plot_nps.add("dot", radial_freq, radial_noise);
            plot_nps.add("line", rebinned_NPS.freq, rebinned_NPS.val);
            plot_nps.show();
        }
        // ---------------------------------------------------------------------


        // Write the results to the ImageJ log
        if (this.log_results) logNPS(rebinned_NPS, "NPS (HU^2.mm^2)");

        nps_result = rebinned_NPS;
    }


    /**
     * Method to calculate the ensemble 2D NPS from an image using a series of square
     * record ROIs arranged around the perimeter of a circle with specified centre and radius.
     * For each of these records:
     * <ul>
     *     <li>A 2D second-order polynomial fit is carried out and subtracted to remove any signal</li>
     *     <li>The record is Fourier transformed, resulting in a 2D array of complex numbers</li>
     *     <li>The magnitude of each complex value is calculated and squared</li>
     *     <li>This result is added to the ensemble</li>
     * </ul>
     *
     * @param im_proc The image data for analysis
     * @param radius_px The radius to use when positioning the record ROIs
     * @param cal The pixel calibration of the image
     * @param centre The coordinates of the centre of the circle on the perimeter of which the record ROIs are to be placed
     * @param overlay The overlay that will display each record ROI
     * @return The ensemble 2D NPS
     */
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


            this.record_stdev_sum += record.getStatistics().stdDev;
        }
        return ensemble_nps;
    }


    /**
     * Method to multiply every element in a 2D double array by a value.
     *
     * @param values Two-dimensional array whose elements are to be changed
     * @param mult_val The factor to change the array values by
     * @return Two-dimensional array of the same dimensions as the input
     */
    private static double[][] multiply2dDoubleArray(final double[][] values, final double mult_val) {
        int width = values.length;
        int height = values[0].length;

        final double[][] result = new double[width][height];

        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                result[i][j] = values[i][j] * mult_val;
            }
        }
        return result;
    }


    /**
     * Method to obtain the calibrated pixel values contained in an ROI positioned
     * on an ImageProcessor.
     *
     * @param ip The ImageProcessor from which to get the values from
     * @param roi The Roi to use to obtain the values - assumed to be rectangular
     * @return Two-dimensional array containing the calibrated pixel values of the Roi
     */
    private static double[][] getCalibratedPixelValues(ImageProcessor ip, Roi roi) {
        int width = (int) roi.getFloatWidth();
        int height = (int) roi.getFloatHeight();

        final double[][] result = new double[width][height];

        for (Point p : roi) {
            result[p.x][p.y] = ip.getPixelValue(p.x, p.y);
        }

        return result;
    }


    /**
     * Method to calculate the square of the absolute value of a 2D array of complex numbers
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
     * Method to add together a pair of 2D double arrays.
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
     * Method to calculate the distance between two cartesian coordinates
     *
     * @param x1 X location of point 1
     * @param y1 Y location of point 2
     * @param x2 X location of point 1
     * @param y2 y location of point 2
     * @return The straight-line distance between point 1 and 2
     */
    public static double calculateDistanceBetweenPoints(double x1, double y1, double x2,double y2) {
        // Calculate the distance between two points in cartesian space.
        return Math.sqrt( (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }


    /**
     * Method to rebin NPS data to a regular frequency spacing.
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
     * Method to write NPS data out to the ImageJ log.
     *
     * @param nps The NPS data to be written to the log
     * @param nps_label The label to use for the NPS value column
     */
    private void logNPS(NPS_result nps, String nps_label) {
        IJ.log("Freq (mm^-1)," + nps_label);

        for (int i=0; i<nps.freq.length; i++) {
            IJ.log("" + this.two_dp.format(nps.freq[i]) + ", " + this.three_dp.format(nps.val[i]));
        }
        IJ.log("\n");
    }


    /**
     * Definition of a NPS_result class used to store NPS results
     */
    static class NPS_result {
        double[] freq;
        double[] val;
    }
}
