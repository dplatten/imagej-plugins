package QuantitativeIQ;
// ImageJ plugin to calculate the detectability index of a circular object present in an image.
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

    private double task_radius = 1.5; // mm
    private double task_contrast = 1.0; // In HU

    private DecimalFormat one_dp = new DecimalFormat("0.0");
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
        TTF.setLogging(false); // Switch off writing to the ImageJ log
        int result = TTF.setup("defaults", this.imp);
        if (result == 0) {
            IJ.error("There was a problem setting up the TTF routine");
            return;
        }

        TTF.run(ip);
        TTF_result ttf_result = TTF.getTTFResult();
        CNR_result cnr_result = TTF.getCNRResult();


        // Then run the NPS on the same image using default settings
        this.imp.deleteRoi(); // Remove any currently set ROI
        NoisePowerSpectrum_ NPS = new NoisePowerSpectrum_();
        NPS.setLogging(false); // Switch off writing to the ImageJ log
        result = NPS.setup("defaults", this.imp);
        if (result == 0) {
            IJ.error("There was a problem setting up the NPS routine");
            return;
        }
        NPS.run(ip);
        NPS_result nps_result = NPS.getNPSResult();


        // Obtain the Nyquist frequency of the image
        Calibration cal = this.imp.getCalibration();
        double nyquist = 1.0 / (2.0 * cal.pixelWidth);


        // Obtain the eye model data
        EyeModel_ eye_model = new EyeModel_();
        eye_model.setLogging(false); // Switch off writing to the ImageJ log
        eye_model.setFreqCutoff(nyquist);
        eye_model.setReconFOV(cal.pixelWidth * imp.getWidth());
        eye_model.setImagePixelWidth(imp.getWidth());
        eye_model.setDisplayZoomFactor(3.0);
        eye_model.setDisplayPixelPitch(0.2);
        eye_model.setup(EyeModel_.SOLOMON_ET_AL);
        eye_model.run("normal");
        Eye_result E = eye_model.getEyeResult();


        // Obtain the task function
        this.task_contrast = cnr_result.contrast; // HU
        double task_x_inc = 0.01; // mm


        // Ask the user what task settings they want to use
        GenericDialog gd = new GenericDialog("Task function settings");
        gd.addNumericField("Task radius (mm)", this.task_radius, 2);
        gd.addNumericField("Task contrast (HU)", this.task_contrast, 1);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }
        else {
            this.task_radius = gd.getNextNumber();
            this.task_contrast = gd.getNextNumber();
        }

        TaskFunction_ task_function = new TaskFunction_();
        task_function.setTaskRadius(this.task_radius);
        task_function.setTaskContrast(this.task_contrast);
        task_function.disc_function();

        UnivariateFunction task_fn = task_function.getTaskFunction();


        // Plot the task function
        double x_min = -2.0 * this.task_radius;
        double x_max = 2.0 * this.task_radius;
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


        // Log the results
        IJ.log("Detectability results" +
                "\nDetectability,"    + this.two_dp.format(detectability) +
                "\nTask radius (mm)," + this.two_dp.format(this.task_radius) +
                "\nContrast (HU),"    + this.two_dp.format(cnr_result.contrast) +
                "\nCNR,"              + this.two_dp.format(cnr_result.cnr) +
                "\nNoise (HU),"       + this.two_dp.format(cnr_result.noise) +
                "\nTTF 50 (1/mm),"    + this.two_dp.format(TTF.getTTF50()) +
                "\nTTF 10 (1/mm),"    + this.two_dp.format(TTF.getTTF10())
        );

        IJ.log("\nEye model settings used" +
                "\nViewing distance (mm),"            + this.one_dp.format(eye_model.getViewingDistance()) +
                "\nRecon field of view (mm),"         + this.two_dp.format(eye_model.getReconFOV()) +
                "\nImage pixel width,"                + this.three_dp.format(eye_model.getImagePixelWidth()) +
                "\nDisplay monitor pixel pitch (mm)," + this.three_dp.format(eye_model.getDisplayPixelPitch()) +
                "\nDisplay zoom factor,"              + this.one_dp.format(eye_model.getDisplayZoomFactor()) +
                "\nSize of displayed image (mm),"     + this.two_dp.format(eye_model.getDisplayedImageSize()) +
                "\nRadial (1/mm) to angular frequency (c/deg) factor," + this.three_dp.format(eye_model.getFreqToAngFreqFactor())
        );


        // Log the TTF, NPS, Eye, and Task values
        IJ.log("\nDetectability components" +
                "\nFreq (mm^-1),TTF,NPS,Eye,Task"
        );
        double f = 0.0;
        while (f <= nyquist) {
            IJ.log( this.three_dp.format(f) + "," +
                    this.three_dp.format(TTF_fn.value(f)) + "," +
                    this.three_dp.format(NPS_fn.value(f)) + "," +
                    this.three_dp.format(E_fn.value(f))   + "," +
                    this.three_dp.format(W_fn.value(f))
            );
            f += 0.01;
        }
        IJ.log("\n");
    }
}
