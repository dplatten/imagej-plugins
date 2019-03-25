package QuantitativeIQ;// ImageJ plugin to calculate the task transfer function of a circular object
// present in an image.
//
// The routine will run on a user-selected oval region of interest (ROI). If no
// ROI is selected then the routine exits.
//
//package ij.plugin;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

public class TaskFunction_ {
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

        FastFourierTransformer fourier_transformer = new FastFourierTransformer(DftNormalization.UNITARY);
        Complex[] fft = fourier_transformer.transform(y_val, TransformType.FORWARD);

        double[] freq_scale = new double[n_elements];
        double[] abs_fft = new double[n_elements];
        for (int i=0; i<n_elements; i++) {
            freq_scale[i] = i * f_res;
            abs_fft[i] = fft[i].abs();
        }

        UnivariateInterpolator linear_interpoloator = new LinearInterpolator();

        task_function_fft = linear_interpoloator.interpolate(freq_scale, abs_fft);
    }
}