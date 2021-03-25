

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


// ImageJ plugin to calculate the task transfer function of a circular object
// present in an image.
//
// The routine will run on a user-selected oval region of interest (ROI). If no
// ROI is selected then the routine exits.
//
import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.util.MathArrays;

import java.text.DecimalFormat;

/**
 *
 * ImageJ plugin to calculates the radial frequency response of a disc-shaped object of a user-specified diameter and
 * contrast.
 * <br><br>
 * The plugin calculates the absolute value of the Fourier transform of a top-hat function (the profile through the
 * disc), scaled by length / N, where the top-hat function has a width equal to the diameter of the disc, and a height equal to the disc
 * contrast.
 * <br><br>
 * See this stack exchange discussion regarding the scaling: <a href="https://dsp.stackexchange.com/questions/31984/matlab-tt-fft-and-tt-ifft-scaling">https://dsp.stackexchange.com/questions/31984/matlab-tt-fft-and-tt-ifft-scaling</a>
 * <br><br>
 * The Fourier transform of a rectangular function is a sinc function with the height scaled by the height of the
 * rectangular function multiplied by the width. In this case height = task_contrast, and width = 2*task_radius. This
 * plugin can also calculate the sinc function, which should be mathematically equal to the Fourier transform result.
 * <br><br>
 * See:
 * <ul>
 *     <li><a href="https://www.cpp.edu/~zaliyazici/ece307/Fourier%20Transform.pdf">https://www.cpp.edu/~zaliyazici/ece307/Fourier%20Transform.pdf</a></li>
 *     <li><a href="https://www.cpp.edu/~zaliyazici/ece307/Fourier%20Transform.pdf">https://www.math.ubc.ca/~feldman/m267/ft.pdf</a> (example 3)</li>
 *     <li><a href="https://www.cpp.edu/~zaliyazici/ece307/Fourier%20Transform.pdf">http://www.thefouriertransform.com/pairs/box.php</a></li>
 * </ul>
 */
public class TaskFunction_ implements PlugIn {
    private String task_type = "cylinder";

    private double task_radius = 1.5; // Radius of task (mm)
    private double task_contrast = 1.0; // The contrast of the task

    private double cutoff_f = 2.0; // Maximum frequency for which to calculate the task when logging
    private double f_res = 0.01; // Increment for the task function fft (mm^-1)

    private DecimalFormat three_dp = new DecimalFormat("0.000");

    private UnivariateFunction task_function;
    private UnivariateFunction task_function_fft;
    private UnivariateFunction task_function_sinc;

    private boolean log_results = true;
    private boolean show_plots = true;


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


    /** Method to set the frequency cut off for logging
     *
     * @param new_freq The new frequency cut off (per mm)
     */
    void setFreqCutoff(double new_freq) {
        this.cutoff_f = new_freq;
    }


    /** Method to return the task function
     *
     * @return The task function
     */
    // Defined to enable other plugins to access the task function
    UnivariateFunction getTaskFunction() {
        return this.task_function;
    }


    /** Method to return the Fourier transform of the task function
     *
     * @return The Fourier transform of the task function
     */
    // Defined to enable other plugins to access the task function TTF
    UnivariateFunction getTaskFunctionFFT() {
        return this.task_function_fft;
    }


    /** Method to return the sinc of the task function - this mathematically equivalent
     * to the FFT when the function is rectangular
     *
     * @return The task function
     */
    // Defined to enable other plugins to access the sinc of the task function
    UnivariateFunction getTaskFunctionSinc() {
        return this.task_function_sinc;
    }


    /** Method to set the task radius
     *
     * @param new_radius The new task radius in mm
     */
    void setTaskRadius(double new_radius) { this.task_radius = new_radius; }


    /** Method to get the current task radius
     *
     * @return The current task radius in mm
     */
    double getTaskRadius() { return this.task_radius; }


    /** Method to set the task contrast
     *
     * @param new_contrast The task contrast
     */
    void setTaskContrast(double new_contrast) { this.task_contrast = new_contrast; }


    /** Method to get the current task contrast
     *
     * @return The current task contrast
     */
    double getTaskContrast() { return this.task_contrast; }


    /** Method to set the frequency increment used for plotting the task function
     *
     * @param new_f_res The new frequency increment in 1/mm
     */
    void setFreqRes(double new_f_res) { this.f_res = new_f_res; }


    /** Method to get the current frequency increment
     *
     * @return The current frequency increment
     */
    double getFreqRes() { return this.f_res; }


    /**
     *
     * @param arg Not currently in use
     */
    public void run(String arg) {
        // Requires at least ImageJ version 1.49v as this is when the ij.gui.Plot
        // class was introduced:
        // https://imagej.nih.gov/ij/developer/api/ij/gui/Plot.html
        if (IJ.versionLessThan("1.49v")) return;

        // Produce a task function for the user, but don't do anything with it
        double pixel_size = 0.5;

        GenericDialog gd = new GenericDialog("Task function options");
        gd.addNumericField("Task diameter (mm)", this.task_radius * 2.0, 2);
        gd.addNumericField("Task contrast (HU)", this.task_contrast, 1);
        gd.addNumericField("Image pixel size (mm)", pixel_size, 3);

        gd.showDialog();
        if(gd.wasCanceled()) {
            return;
        }

        this.task_radius = gd.getNextNumber() / 2.0;
        this.task_contrast = gd.getNextNumber();
        pixel_size = gd.getNextNumber();

        double nyquist = 1.0 / (2.0*pixel_size);

        disc_function();
        disc_function_fft();
        disc_function_sinc();

        if (this.show_plots) {
            // Plot the task function
            double task_x_inc = 0.01;
            double x_min = -2.0 * this.task_radius;
            double x_max = 2.0 * this.task_radius;
            int n_elements = (int) ((x_max - x_min) / task_x_inc) + 1;

            double[] x_pos = new double[n_elements];
            double[] y_val = new double[n_elements];

            for (int i = 0; i < n_elements; i++) {
                x_pos[i] = x_min + i * task_x_inc;
                y_val[i] = this.task_function.value(x_pos[i]);
            }
            Plot plot_task = new Plot("Imaging task", "Distance (mm)", "Pixel value");
            plot_task.add("line", x_pos, y_val);
            String task_labels = "Imaging task\t";
            plot_task.addLegend(task_labels);
            plot_task.show();


            // Plot the FFT and sinc of the task function
            int n_felements = (int) Math.ceil(nyquist / this.f_res) + 1;

            double[] f_pos = new double[n_felements];
            double[] fft_val = new double[n_felements];
            double[] sinc_val = new double[n_felements];
            for (int i = 0; i < n_felements; i++) {
                f_pos[i] = i * this.f_res;
                fft_val[i] = this.task_function_fft.value(f_pos[i]);
                sinc_val[i] = this.task_function_sinc.value(f_pos[i]);
            }
            Plot plot_task_fft = new Plot("Task function", "Frequency (mm^-1)", "Value");
            plot_task_fft.add("line", f_pos, fft_val);
            String task_fft_labels = "|task function| up to Nyquist\t";
            plot_task_fft.setColor("Red");
            plot_task_fft.add("line", f_pos, sinc_val);
            task_fft_labels += "|sinc| from task size and contrast up to Nyquist\t";
            plot_task_fft.addLegend(task_fft_labels);
            plot_task_fft.show();
        }

        if (this.log_results) logResults("FFT of radius " + this.task_radius + " mm; " + this.task_contrast + " contrast");
    }


    /** Method to produce a function that returns the y value of a rectangular function at position x
     *
     */
    void disc_function() {
        this.task_function = new UnivariateFunction() {
            @Override
            public double value(double x) {
                double y = 0.0;
                if (Math.abs(x) <= task_radius) y = task_contrast;
                return y;
            }
        };
    }


    /** Method to produce a function of the Fourier transform of the current task function
     *
     */
    void disc_function_fft() {
        double pos_range = 1.0 / this.f_res;
        double min_x = -1.0 * pos_range/2.0;

        int n_elements = 4096;

        double x_inc = pos_range / n_elements;

        double[] x_pos = new double[n_elements];
        double[] y_val = new double[n_elements];

        for (int i=0; i<n_elements; i++) {
            x_pos[i] = min_x + i*x_inc;
            y_val[i] = this.task_function.value(x_pos[i]);
        }

        FastFourierTransformer fourier_transformer = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] fft = fourier_transformer.transform(y_val, TransformType.FORWARD);

        double[] freq_scale = new double[n_elements];
        double[] abs_fft = new double[n_elements];
        for (int i=0; i<n_elements; i++) {
            freq_scale[i] = i * this.f_res;
            abs_fft[i] = fft[i].abs() * pos_range / n_elements; // Scaling the result to length / N. See https://dsp.stackexchange.com/questions/31984/matlab-tt-fft-and-tt-ifft-scaling
        }

        UnivariateInterpolator linear_interpoloator = new LinearInterpolator();
        this.task_function_fft = linear_interpoloator.interpolate(freq_scale, abs_fft);
    }


    /** Method to produce a function of the sinc function of the current task function
     *
     */
    private void disc_function_sinc() {
        // The Fourier transform of a rectangular function is a sinc function with
        // the height scaled by the height of the rectangular function multiplied
        // by the width. In this case height = task_contrast, and width = 2*task_radius
        // https://www.cpp.edu/~zaliyazici/ece307/Fourier%20Transform.pdf
        // https://www.math.ubc.ca/~feldman/m267/ft.pdf example 3
        // http://www.thefouriertransform.com/pairs/box.php
        int n_elements = 4096;

        double[] freq_scale = new double[n_elements];
        double[] abs_sinc = new double[n_elements];

        abs_sinc[0] = 1.0; // Manually set zero frequency value to avoid divide by zero error in loop
        for (int i=1; i<n_elements; i++) {
            freq_scale[i] = i * this.f_res;
            double x = freq_scale[i] * Math.PI * 2.0*this.task_radius;
            abs_sinc[i] = Math.abs(Math.sin(x) / x);
        }

        // Scale the height of the sinc function by the height * width of the
        // rectangular function.
        MathArrays.scaleInPlace(2.0*this.task_radius * this.task_contrast, abs_sinc);

        UnivariateInterpolator linear_interpolator = new LinearInterpolator();
        this.task_function_sinc = linear_interpolator.interpolate(freq_scale, abs_sinc);
    }


    /** Method to write the task function FFT results to the ImageJ log
     *
     * @param label The label for the task function FFT values
     */
    private void logResults(String label) {
        IJ.log("Task function using " + this.three_dp.format(this.cutoff_f) + " mm^-1 cutoff  frequency");
        IJ.log("Freq (mm^-1)," + label);

        double f = 0.0;
        while (f <= this.cutoff_f) {
            IJ.log( this.three_dp.format(f) + "," + this.three_dp.format(this.task_function_fft.value(f)));
            f += this.f_res;
        }
        IJ.log("\n");
    }
}