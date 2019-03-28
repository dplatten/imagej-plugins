package QuantitativeIQ;// ImageJ plugin to calculate the task transfer function of a circular object
// present in an image.
//
// The routine will run on a user-selected oval region of interest (ROI). If no
// ROI is selected then the routine exits.
//
//package ij.plugin;

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

public class TaskFunction_ implements PlugIn {
    private String task_type = "cylinder";

    private double task_radius = 1.5; // Radius of task (mm)
    private double task_contrast = 1.0; // The contrast of the task
    private double f_res = 0.01; // Increment for the task function fft (mm^-1)

    private UnivariateFunction task_function;
    private UnivariateFunction task_function_fft;

    // Defined to enable other plugins to access the task function
    UnivariateFunction getTaskFunction() {
        return task_function;
    }

    // Defined to enable other plugins to access the task function TTF
    UnivariateFunction getTaskFunctionTTF() {
        return task_function_fft;
    }

    void setTaskRadius(double new_radius) { task_radius = new_radius; }
    double getTaskRadius() { return task_radius; }

    void setTaskContrast(double new_contrast) { task_contrast = new_contrast; }
    double getTaskContrast() { return task_contrast; }

    void setFreqRes(double new_f_res) { f_res = new_f_res; }
    double getFreqRes() { return f_res; }


    public void run(String arg) {
        // Produce a task function for the user, but don't do anything with it
        double pixel_size = 0.5;

        GenericDialog gd = new GenericDialog("Task function options");
        gd.addNumericField("Task radius (mm)", task_radius, 2);
        gd.addNumericField("Task contrast (HU)", task_contrast, 1);
        gd.addNumericField("Image pixel size (mm)", pixel_size, 3);

        gd.showDialog();
        if(gd.wasCanceled()) {
            IJ.error("Cancelled");
            return;
        }

        task_radius = gd.getNextNumber();
        task_contrast = gd.getNextNumber();
        pixel_size = gd.getNextNumber();

        double nyquist = 1.0 / (2.0*pixel_size);

        disc_function();
        disc_function_fft();

        // Plot the task function
        double task_x_inc = 0.01;
        double x_min = -2.0 * task_radius;
        double x_max = 2.0 * task_radius;
        int n_elements = (int) ((x_max - x_min) / task_x_inc) + 1;

        double[] x_pos = new double[n_elements];
        double[] y_val = new double[n_elements];

        for (int i=0; i<n_elements; i++) {
            x_pos[i] = x_min + i*task_x_inc;
            y_val[i] = task_function.value(x_pos[i]);
        }
        Plot plot_task = new Plot("Task function", "Distance (mm)", "Pixel value");
        plot_task.add("line", x_pos, y_val);
        String task_labels = "Task function\t";
        plot_task.addLegend(task_labels);
        plot_task.show();


        // Plot the FFT of the task function
        int n_felements = (int) Math.ceil(nyquist / f_res) + 1;

        double[] f_pos = new double[n_felements];
        double[] fft_val = new double[n_felements];
        for (int i=0; i<n_felements; i++) {
            f_pos[i] = i*f_res;
            fft_val[i] = task_function_fft.value(f_pos[i]);
        }
        Plot plot_task_fft = new Plot("Task function fft", "Frequency (mm^-1)", "Value");
        plot_task_fft.add("line", f_pos, fft_val);
        String task_fft_labels = "FFT of task function up to Nyquist\t";
        plot_task_fft.addLegend(task_fft_labels);
        plot_task_fft.show();
    }


    public void disc_function() {
        task_function = new UnivariateFunction() {
            @Override
            public double value(double x) {
                double y = 0.0;
                if (Math.abs(x) <= task_radius) y = task_contrast;
                return y;
            }
        };
    }


    public void disc_function_fft() {
        double pos_range = 1.0 / f_res;
        double min_x = -1.0 * pos_range/2.0;

        int n_elements = 4096;

        double x_inc = pos_range / n_elements;

        double[] x_pos = new double[n_elements];
        double[] y_val = new double[n_elements];

        for (int i=0; i<n_elements; i++) {
            x_pos[i] = min_x + i*x_inc;
            y_val[i] = task_function.value(x_pos[i]);
        }

        FastFourierTransformer fourier_transformer = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] fft = fourier_transformer.transform(y_val, TransformType.FORWARD);

        double[] freq_scale = new double[n_elements];
        double[] abs_fft = new double[n_elements];
        for (int i=0; i<n_elements; i++) {
            freq_scale[i] = i * f_res;
            abs_fft[i] = fft[i].abs() * pos_range / n_elements; // Scaling the result to length / N. See https://dsp.stackexchange.com/questions/31984/matlab-tt-fft-and-tt-ifft-scaling
        }

        UnivariateInterpolator linear_interpoloator = new LinearInterpolator();

        task_function_fft = linear_interpoloator.interpolate(freq_scale, abs_fft);
    }
}