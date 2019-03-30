package QuantitativeIQ;
// ImageJ plugin to calculate the noise power spectrum
//
// The routine will run on a user-selected oval region of interest (ROI). If no
// ROI is selected then the routine uses the whole image.
//
import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.util.MathArrays;

import java.awt.*;
import java.text.DecimalFormat;

/**
 *
 */
public class NoisePowerSpectrum_ implements PlugInFilter {
    private ImagePlus imp;
    private Roi main_roi;

    private int num_records = 8;
    private int record_size = 64;
    private double record_radius = 100.0; // in mm

    private double rebin_f_inc = 0.05; // in mm^-1

    private static NPS_result nps_result;

    private DecimalFormat two_dp = new DecimalFormat("0.00");
    private DecimalFormat three_dp = new DecimalFormat("0.000");

    /**
     *
     * @return
     */
    // Defined to enable other plugins to access the TTF results
    NPS_result getNPSResult() {
        return nps_result;
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
            IJ.error("Please open an image");
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
                IJ.log("NoisePowerSpectrum using ROIs " + this.record_radius + " mm from centre");
            }
            else {
                // Assume the argument contains the distance from ROI centre to record in mm.
                this.record_radius = Float.parseFloat(arg);
                IJ.log("NoisePowerSpectrum using ROIs " + this.record_radius + " mm from centre");
            }
        }
        else {
            GenericDialog gd = new GenericDialog("Processing options");
            gd.addNumericField("Distance from centre to records (mm)", this.record_radius, 0);
            gd.addNumericField("Number of records", this.num_records, 0);
            gd.addNumericField("Record size (pixels), must be a power of 2", this.record_size, 0);
            gd.addNumericField("Frequency rebin increment (mm^-1)", this.rebin_f_inc, 2);

            gd.showDialog();
            if(gd.wasCanceled()) {
                IJ.error("Plugin cancelled");
                return(0);
            }

            this.record_radius = gd.getNextNumber();
            this.num_records = (int) gd.getNextNumber();
            this.record_size = (int) gd.getNextNumber();
            this.rebin_f_inc = gd.getNextNumber();
        }
        return DOES_ALL;
    }

    /**
     *
     * @param ip
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
        // Calculate the NPS of each record
        for (int m=0; m<this.num_records; m++) {
            double current_angle = (2.0 * Math.PI) / this.num_records * m;

            double x_pos = record_radius_px * Math.cos(current_angle);
            double y_pos = record_radius_px * Math.sin(current_angle);

            Roi record_roi = new Roi(cal.getRawX(centre[0]) + x_pos - (this.record_size/2.0), cal.getRawY(centre[1]) + y_pos - (this.record_size/2.0), this.record_size, this.record_size);

            overlay.add(record_roi);
            overlay.setStrokeColor(Color.green);

            // Set the roi on the main image
            this.imp.setRoi(record_roi);

            // Carry out a 2d second-order polynomial fit on a new ImagePlus
            // that just contains the current record pixels. The Fit Polynomial
            // command subtracts the fit from the original data.
            ImagePlus record = this.imp.crop();
            IJ.run(record, "Fit Polynomial", "x=2 y=2 mixed=0");
            record.setRoi(0, 0, record.getWidth(), record.getHeight());
            Roi whole_record_roi = record.getRoi();
            double[][] record_vals = getCalibratedPixelValues(record.getProcessor(), whole_record_roi);

            // Fourier transform the record
            Complex[][] record_vals_complex = convertDoubleToComplex(record_vals);
            FastFourierTransformer fourier_transformer = new FastFourierTransformer(DftNormalization.STANDARD);
            Complex[][] nps = (Complex[][]) fourier_transformer.mdfft(record_vals_complex, TransformType.FORWARD);

            // Calculate magnitude and square it
            double[][] nps_abs_sqrd = absSqrdOfComplex2dArray(nps);

            // Add the current nps_abs_sqrd to the ensemble
            ensemble_nps = add2dArrays(ensemble_nps, nps_abs_sqrd);
        }
        //---------------------------------------------------------------------


        //---------------------------------------------------------------------
        // Multiply the ensemble by dx, dy and divide by M, Nx and Ny
        double multiplication_factor = pixel_size_in_mm * pixel_size_in_mm / this.num_records / this.record_size / this.record_size;
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
        // The results of the Fourier transform are arranged a bit strangely.
        // Element 0 is the zero frequency component.
        // Element 1 is the smallest non-zero positive frequency, 1/(record_size * pixel_size_in_mm)
        // Negative frequencies are stored in the reverse order of positive frequencies, ranging
        // from the highest to lowest negative frequencies.
        // Only need to use the top half of the FFT, as it is symmetrical
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

        MathArrays.sortInPlace(radial_freq, radial_nps);

        Plot plot_nps = new Plot("Noise power spectrum", "Frequency (mm^-1)", "Noise power (HU^2.mm^2)");
        plot_nps.setAxisYLog(true);
        plot_nps.add("dot", radial_freq, radial_nps);
        plot_nps.show();

        double IEC_f_bin = 0.01 / pixel_size_in_mm;
        NPS_result rebinned_NPS = rebinNPS(radial_freq, radial_nps, rebin_f_inc, IEC_f_bin);

        plot_nps.add("line", rebinned_NPS.freq, rebinned_NPS.val);
        // ---------------------------------------------------------------------


        // Write the results to the ImageJ log
        logNPS(rebinned_NPS, "NPS (HU^2.mm^2)");

        nps_result = rebinned_NPS;
    }


    /**
     *
     * @param values
     * @param mult_val
     * @return
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
     * @param values
     * @return
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
     * @param values
     * @return
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
     * @param ip
     * @param roi
     * @return
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
     * @param values
     * @return
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
     * @param arr1
     * @param arr2
     * @return
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
     * @param freq
     * @param val
     * @param f_inc
     * @param f_bin
     * @return
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
     * @param nps
     * @param nps_label
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
