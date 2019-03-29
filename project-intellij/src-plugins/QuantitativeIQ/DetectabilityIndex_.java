package QuantitativeIQ;
// ImageJ plugin to calculate the detectability index of a circular object
// present in an image.
//
import QuantitativeIQ.TaskTransferFunction_.TTF_result;
import QuantitativeIQ.TaskTransferFunction_.CNR_result;
import QuantitativeIQ.NoisePowerSpectrum_.NPS_result;
import QuantitativeIQ.EyeModel_.Eye_result;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
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

/**
 *
 */
public class DetectabilityIndex_ implements PlugInFilter {
    private ImagePlus imp;

    private DecimalFormat two_dp = new DecimalFormat("0.00");
    private DecimalFormat three_dp = new DecimalFormat("0.000");

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

        return DOES_ALL;
    }

    /**
     *
     * @param ip
     */
    public void run(ImageProcessor ip) {
        // Calculate the detectability index according to equation 2 in http://dx.doi.org/10.1118/1.4923172,
        // adapted for radial frequency rather than u, v.
        //
        // Detectability index needs:
        //  Task function (W)
        //  Task Transfer Function (TTF)
        //  Noise Power Spectrum (NPS)
        //  Eye model (E)

        // First run the TTF analysis using default settings
        TaskTransferFunction_ TTF = new TaskTransferFunction_();
        int result = TTF.setup("defaults", imp);
        if (result == 0) {
            IJ.error("There was a problem setting up the TTF routine");
            return;
        }
        TTF.run(ip);
        TTF_result ttf_result = TTF.getTTFResult();
        CNR_result cnr_result = TTF.getCNRResult(); // We're going to use the cnr_result.contrast when creating the task function


        // Then run the NPS on the same image using default settings
        NoisePowerSpectrum_ NPS = new NoisePowerSpectrum_();
        result = NPS.setup("defaults", imp);
        if (result == 0) {
            IJ.error("There was a problem setting up the NPS routine");
            return;
        }
        NPS.run(ip);
        NPS_result nps_result = NPS.getNPSResult();


        // Obtain the Nyquist frequency of the image
        Calibration cal = imp.getCalibration();
        double nyquist = 1.0 / (2.0 * cal.pixelWidth);


        // Obtain the eye model data
        EyeModel_ eye_model = new EyeModel_();
        eye_model.setup(two_dp.format(nyquist), 999);
        eye_model.setReconFOV(cal.pixelWidth * imp.getWidth());
        eye_model.setImagePixelWidth(imp.getWidth());
        eye_model.setDisplayZoomFactor(3.0);
        eye_model.setDisplayPixelPitch(0.02);
        eye_model.run("normal");
        Eye_result E = eye_model.getEyeResult();


        // Obtain the task function
        double task_radius = 1.5; // mm
        double task_contrast = cnr_result.contrast; // HU
        double task_x_inc = 0.01; // mm

        // Ask the user what task settings they want to use
        GenericDialog gd = new GenericDialog("Task function settings");
        gd.addNumericField("Task radius (mm)", task_radius, 2);
        gd.addNumericField("Task contrast (HU)", task_contrast, 1);
        gd.showDialog();
        if (! gd.wasCanceled()) { // if the dialogue was not cancelled
            task_radius = gd.getNextNumber();
            task_contrast = gd.getNextNumber();
        }

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


        // Obtain the FFT of the task function, W_fn
        task_function.disc_function_fft();
        UnivariateFunction W_fn =  task_function.getTaskFunctionTTF();

        // Plot the FFT of the task function
        double f_inc = task_function.getFreqRes();
        int n_felements = (int) Math.ceil(nyquist / f_inc) + 1;

        double[] f_pos = new double[n_felements];
        double[] fft_val = new double[n_felements];
        for (int i=0; i<n_felements; i++) {
            f_pos[i] = i*f_inc;
            fft_val[i] = W_fn.value(f_pos[i]);
        }
        Plot plot_task_fft = new Plot("Task function fft", "Frequency (mm^-1)", "Value");
        plot_task_fft.add("line", f_pos, fft_val);
        String task_fft_labels = "FFT of task function up to Nyquist\t";
        plot_task_fft.addLegend(task_fft_labels);
        plot_task_fft.show();


        // Calculate the detectability index using a modified version of equation 2 in http://dx.doi.org/10.1118/1.4923172
        // We are using radial functions here, so it's only a single integral over the radial frequency, rather than
        // a double integral over u and v.
        UnivariateInterpolator linear_interpolator = new LinearInterpolator();
        UnivariateFunction TTF_fn = linear_interpolator.interpolate(ttf_result.freq, ttf_result.val);
        UnivariateFunction E_fn = linear_interpolator.interpolate(E.freq, E.val);
        UnivariateFunction NPS_fn = linear_interpolator.interpolate(nps_result.freq, nps_result.val);

        // Using equation 2 in http://dx.doi.org/10.1118/1.4923172
        UnivariateFunction numerator = FunctionUtils.multiply(W_fn, W_fn, TTF_fn, TTF_fn, E_fn, E_fn);
        UnivariateFunction denominator = FunctionUtils.multiply(W_fn, W_fn, TTF_fn, TTF_fn, NPS_fn, E_fn, E_fn, E_fn, E_fn);

        // Now integrate the numerator and denominator
        UnivariateIntegrator simpson_integrator = new SimpsonIntegrator();
        double numerator_int = simpson_integrator.integrate(100000, numerator, 0, nyquist);
        double denominator_int = simpson_integrator.integrate(100000, denominator, 0, nyquist);

        double detectability = Math.sqrt( Math.pow(numerator_int, 2.0) / denominator_int );

        IJ.log("Detectability is: " + detectability);
        IJ.log("Numerator is: " + numerator_int);
        IJ.log("Denominator is: " + denominator_int);

        // Log TTF, NPS, Eye, Task, eye noise
        IJ.log("Freq (mm^-1),TTF,NPS,Eye,Task");
        double f = 0.0;
        while (f <= nyquist) {
            IJ.log( three_dp.format(f) + "," + TTF_fn.value(f) + "," + NPS_fn.value(f) + "," + E_fn.value(f) + "," + W_fn.value(f));
            f += 0.01;
        }
    }
}
