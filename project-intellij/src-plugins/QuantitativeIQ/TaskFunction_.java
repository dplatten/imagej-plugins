package QuantitativeIQ;// ImageJ plugin to calculate the task transfer function of a circular object
// present in an image.
//
// The routine will run on a user-selected oval region of interest (ROI). If no
// ROI is selected then the routine exits.
//
//package ij.plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.measure.Calibration;
import ij.measure.CurveFitter;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.FiniteDifferencesDifferentiator;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.util.MathArrays;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class TaskFunction_ {
    private String task_type = "cylinder";

    private double task_radius = 1.5; // Radius of task (mm)
    private double task_contrast = 1.0; // The contrast of the task
    private double x_inc = 0.01; // Increment for task function calcuiation (mm)

    private UnivariateFunction task_function;
    private UnivariateFunction task_function_fft;

    private DecimalFormat two_dp = new DecimalFormat("0.00");
    private DecimalFormat three_dp = new DecimalFormat("0.000");

    // Defined to enable other plugins to access the task function
    UnivariateFunction getTaskFunction() {
        return task_function;
    }

    // Defined to enable other plugins to access the task function TTF
    UnivariateFunction getTaskFunctionTTF() {
        return task_function_fft;
    }

    void setTaskRadius(double new_radius) { task_radius = new_radius; }
    void setTaskContrast(double new_contrast) { task_contrast = new_contrast; }


    public int disc_function() {
        UnivariateFunction disc_function = new UnivariateFunction() {
            @Override
            public double value(double x) {
                double y = 0.0;
                if (Math.abs(x) <= task_radius) y = task_contrast;
                return y;
            }
        };
        task_function = disc_function;

        return 0;
    }


    public int disc_function_fft() {
        double min_x = -2.0 * task_radius;
        double max_x =  2.0 * task_radius;
        int n_elements = 2048;
        double x_inc = (max_x - min_x) / (n_elements-1);

        double[] x_pos = new double[n_elements];
        double[] y_val = new double[n_elements];

        for (int i=0; i<n_elements; i++) {
            x_pos[i] = min_x + i*x_inc;
            y_val[i] = task_function.value(x_pos[i]);
        }

        FastFourierTransformer fourier_transformer = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] fft = fourier_transformer.transform(y_val, TransformType.FORWARD);

        double freq_inc = 1.0 / (max_x - min_x);
        double[] freq_scale = new double[n_elements];
        double[] abs_fft = new double[n_elements];
        for (int i=0; i<n_elements; i++) {
            abs_fft[i] = fft[i].abs();
            freq_scale[i] = i * freq_inc;
        }

        UnivariateInterpolator linear_interpoloator = new LinearInterpolator();
        UnivariateFunction fft_function = linear_interpoloator.interpolate(freq_scale, abs_fft);

        task_function_fft = fft_function;

        return 0;
    }
}