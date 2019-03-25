package QuantitativeIQ;// ImageJ plugin to calculate the detectability index of a circular object
// present in an image.
//
//package ij.plugin;

import QuantitativeIQ.TaskTransferFunction_.TTF_result;
import QuantitativeIQ.TaskTransferFunction_.CNR_result;
import QuantitativeIQ.NoisePowerSpectrum_.NPS_result;
import QuantitativeIQ.EyeModel_.Eye_result;
import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import org.apache.commons.math3.analysis.FunctionUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import java.text.DecimalFormat;
import ij.gui.Plot;


public class DetectabilityIndex_ implements PlugInFilter {
    private ImagePlus imp;

    private DecimalFormat two_dp = new DecimalFormat("0.00");
    private DecimalFormat three_dp = new DecimalFormat("0.000");

    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;

        if(this.imp == null) {
            IJ.error("Please open an image");
            return(0);
        }

        return DOES_ALL;
    }

    public void run(ImageProcessor ip) {
        // Detectability index needs
        //  Task Transfer Function (TTF)
        //  Noise Power Spectrum (NPS)
        //  Eye model (E)
        //  Task function (TF)

        // First run the TTF analysis using default settings
        TaskTransferFunction_ ttf = new TaskTransferFunction_();
        int result = ttf.setup("defaults", imp);
        if (result == 0) {
            IJ.error("There was a problem setting up the TTF routine");
            return;
        }
        ttf.run(ip);
        TTF_result ttf_result = ttf.getTTFResult();
        CNR_result cnr_result = ttf.getCNRResult(); // We're going to use the cnr_result.contrast when creating the task function


        // Then run the NPS on the same image using default settings
        NoisePowerSpectrum_ nps = new NoisePowerSpectrum_();
        result = nps.setup("defaults", imp);
        if (result == 0) {
            IJ.error("There was a problem setting up the NPS routine");
            return;
        }
        nps.run(ip);
        NPS_result nps_result = nps.getNPSResult();


        // Work out the internal noise of the eye
        // This is considered to be constant with frequency
        // See equation 6 in Richard and Siewerdsen 2008, http://dx.doi.org/10.1118/1.2988161
        double eye_noise = 0.02 * Math.pow(50.0/100.0, 2.0) * nps_result.val[0];


        // Obtain the Nyquist frequency of the image
        Calibration cal = imp.getCalibration();
        double nyquist = 1.0 / (2.0 * cal.pixelWidth);


        // Obtain the eye model data
        EyeModel_ eye_model = new EyeModel_();
        eye_model.setup(two_dp.format(nyquist));
        eye_model.run();
        Eye_result eye_result = eye_model.getEyeResult();


        // Obtain the task function
        double task_radius = 1.5; // mm
        double task_contrast = cnr_result.contrast; // HU
        double task_x_inc = 0.01; // mm

        TaskFunction_ task_function = new TaskFunction_();
        task_function.setTaskRadius(task_radius);
        task_function.setTaskContrast(task_contrast);
        task_function.disc_function();

        UnivariateFunction task_fn = task_function.getTaskFunction();

        // Plot the task function
        double x_min = -2.0 * task_radius;
        double x_max = 2.0 * task_radius;
        int n_elements = (int) ((x_max - x_min) / task_x_inc) + 1;

        double[] x_pos = new double[n_elements];
        double[] y_val = new double[n_elements];

        for (int i=0; i<n_elements; i++) {
            x_pos[i] = x_min + i*task_x_inc;
            y_val[i] = task_fn.value(x_pos[i]);
        }
        Plot plot_task = new Plot("Task function", "Distance (mm)", "Pixel value");
        plot_task.add("line", x_pos, y_val);
        String task_labels = "Task function\t";
        plot_task.addLegend(task_labels);
        plot_task.show();


        // Obtain the FFT of the task function, W
        task_function.disc_function_fft();
        UnivariateFunction task_fn_fft =  task_function.getTaskFunctionTTF();

        // Plot the FFT of the task function
        double f_inc = task_function.getFreqRes();
        int n_felements = (int) Math.ceil(nyquist / f_inc) + 1;

        double[] f_pos = new double[n_felements];
        double[] fft_val = new double[n_felements];
        for (int i=0; i<n_felements; i++) {
            f_pos[i] = i*f_inc;
            fft_val[i] = task_fn_fft.value(f_pos[i]);
        }
        Plot plot_task_fft = new Plot("Task function fft", "Frequency (mm^-1)", "Value");
        plot_task_fft.add("line", f_pos, fft_val);
        String task_fft_labels = "FFT of task function up to Nyquist\t";
        plot_task_fft.addLegend(task_fft_labels);
        plot_task_fft.show();


        // Calculate the detectability index using a modified version of equation 8 in Richard and Siewerdsen 2008, http://dx.doi.org/10.1118/1.2988161
        // We are using radial functions here, so it's only a single integral over the radial frequency, rather than
        // a double integral over u and v.
        UnivariateInterpolator linear_interpolator = new LinearInterpolator();
        UnivariateFunction ttf_function = linear_interpolator.interpolate(ttf_result.freq, ttf_result.val);
        UnivariateFunction eye_function = linear_interpolator.interpolate(eye_result.freq, eye_result.val);
        UnivariateFunction nps_function = linear_interpolator.interpolate(nps_result.freq, nps_result.val);

        // Using equation 8 in Richard and Siewerdsen 2008, http://dx.doi.org/10.1118/1.2988161
        UnivariateFunction numerator = FunctionUtils.multiply(task_fn_fft, task_fn_fft, ttf_function, ttf_function, eye_function, eye_function);
        UnivariateFunction denominator_a = FunctionUtils.multiply(nps_function, ttf_function, ttf_function, task_fn_fft, task_fn_fft, eye_function, eye_function, eye_function, eye_function);
        UnivariateFunction denominator_b = FunctionUtils.multiply(ttf_function, ttf_function, task_fn_fft, task_fn_fft); // Must multiply denominator_b integral by eye_noise later on...


        // Now integrate the numerator and denominator
        UnivariateIntegrator simpson_integrator = new SimpsonIntegrator();
        double numerator_int = simpson_integrator.integrate(20000, numerator, 0, nyquist);
        double denominator_a_int = simpson_integrator.integrate(20000, denominator_a, 0, nyquist);
        double denominator_b_int = simpson_integrator.integrate(20000, denominator_b, 0, nyquist) * eye_noise;

        double detectability = Math.sqrt( Math.pow(numerator_int, 2.0) / (denominator_a_int + denominator_b_int) );

        IJ.log("Detectability is: " + detectability);
        IJ.log("Numerator is: " + numerator_int);
        IJ.log("Denominator_a is: " + denominator_a_int);
        IJ.log("Denominator_b is: " + denominator_b_int);
        IJ.log("Eye noise is: " + eye_noise);

        // Log TTF, NPS, Eye, Task, eye noise
        IJ.log("Freq (mm^-1),TTF,NPS,Eye,Task,Eye Noise");
        double f = 0.0;
        while (f <= nyquist) {
            IJ.log( three_dp.format(f) + "," + ttf_function.value(f) + "," + nps_function.value(f) + "," + eye_function.value(f) + "," + task_fn_fft.value(f) + "," + eye_noise);
            f += 0.01;
        }
    }
}
