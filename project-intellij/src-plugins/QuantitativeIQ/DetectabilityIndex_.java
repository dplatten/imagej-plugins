package QuantitativeIQ;
// ImageJ plugin to calculate the detectability index of a circular object present in an image.
//
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.frame.RoiManager;
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

    private RoiManager roi_manager = RoiManager.getInstance2();
    private Roi main_roi;

    private boolean show_eye_plots = false;
    private boolean show_task_plots = false;
    private boolean show_nps_plots = false;
    private boolean show_ttf_plots = false;


    /** Method to set up the detectability index plugin
     *
     * @param arg Not currently used
     * @param imp The ImagePlus to use
     * @return 0 if no image is open, DOES_ALL otherwise
     */
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        if (this.imp == null) {
            IJ.error("Please open an image");
            return(0);
        }

        GenericDialog gd = new GenericDialog("Plotting options");
        gd.addCheckbox("Show eye filter plots", this.show_eye_plots);
        gd.addCheckbox("Show task function plots", this.show_task_plots);
        gd.addCheckbox("Show NPS plots", this.show_nps_plots);
        gd.addCheckbox("Show TTF plots", this.show_ttf_plots);
        gd.showDialog();
        if(gd.wasCanceled()) {
            return(1); // Return to the calling method with a non-zero result
        }
        this.show_eye_plots = gd.getNextBoolean();
        this.show_task_plots = gd.getNextBoolean();
        this.show_nps_plots = gd.getNextBoolean();
        this.show_ttf_plots = gd.getNextBoolean();

        this.main_roi = imp.getRoi();
        if (this.roi_manager == null && this.main_roi == null) {
            IJ.error("Please add an ROI around an object or\nadd ROIs to the ROI manager");
            return (0);
        }

        if (this.roi_manager == null) {
            this.roi_manager = new RoiManager();
        }

        if (this.main_roi == null) {
            // There's no ROI drawn, so hopefully there's at least one ROI in the manager
            if (this.roi_manager.getCount() == 0) {
                IJ.error("Please add an ROI around an object or\nadd ROIs to the ROI manager");
                return (0);
            }
        }
        else {
            // If the ROI that is set is not one of the ones in the ROI manager then
            // add it to the ROI manager
            if (this.roi_manager.getSelectedIndexes().length == 0 ) this.roi_manager.addRoi(this.main_roi);
        }

        return DOES_ALL;
    }


    /** Method to run the plugin
     *
     * @param ip The ImageProcessor to use
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


        // Run the NPS on the image
        this.imp.deleteRoi(); // Remove any currently set ROI
        NoisePowerSpectrum_ NPS = new NoisePowerSpectrum_();
        NPS.setPlotting(this.show_nps_plots);
        NPS.setLogging(false); // Switch off writing to the ImageJ log
        int result = NPS.setup("", this.imp);
        if (result == 0) {
            IJ.error("There was a problem setting up the NPS routine");
            return;
        }

        try {
            NPS.run(ip);
        } catch (Exception e) {
            IJ.error("There was an error when running the NPS:\n" + e);
            return;
        }

        this.imp.updateAndDraw();
        NoisePowerSpectrum_.NPS_result nps_result = NPS.getNPSResult();


        // Obtain the Nyquist frequency of the image
        Calibration cal = this.imp.getCalibration();
        double nyquist = 1.0 / (2.0 * cal.pixelWidth);


        // Obtain the eye model data
        EyeModel_ eye_model = new EyeModel_();
        eye_model.setLogging(false); // Switch off writing to the ImageJ log
        eye_model.setPlotting(this.show_eye_plots);
        eye_model.setFreqCutoff(nyquist);
        eye_model.setReconFOV(cal.pixelWidth * imp.getWidth());
        eye_model.setImagePixelWidth(imp.getWidth());
        eye_model.setDisplayZoomFactor(3.0);
        eye_model.setDisplayPixelPitch(0.2);
        eye_model.setup(EyeModel_.SOLOMON_ET_AL);
        eye_model.run("normal");
        EyeModel_.Eye_result E = eye_model.getEyeResult();


        // Ask the user what task size they want to use
        GenericDialog gd = new GenericDialog("Task function settings");
        gd.addNumericField("Task radius (mm)", this.task_radius, 2);
        //gd.addNumericField("Task contrast (HU)", this.task_contrast, 1);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        } else {
            this.task_radius = gd.getNextNumber();
            //this.task_contrast = gd.getNextNumber();
        }


        // Run the TTF routine on all ROIs that are selected in the ROI manager
        int[] selected_rois = this.roi_manager.getIndexes();
        for (int idx : selected_rois) {

            Roi current_roi = this.roi_manager.getRoi(selected_rois[idx]);
            this.imp.setRoi(current_roi);

            // Run the TTF analysis using default settings
            TaskTransferFunction_ TTF = new TaskTransferFunction_();
            TTF.setPlotting(this.show_ttf_plots);
            TTF.setLogging(false); // Switch off writing to the ImageJ log
            result = TTF.setup("", this.imp);
            if (result == 0) {
                IJ.error("There was a problem setting up the TTF routine");
                return;
            }

            TTF.run(ip);
            TaskTransferFunction_.TTF_result ttf_result = TTF.getTTFResult();
            double[] ttf_com = TTF.getCentreOfMass();
            TaskTransferFunction_.CNR_result cnr_result = TTF.getCNRResult();


            // Obtain the task function
            this.task_contrast = cnr_result.contrast; // HU
            double task_x_inc = 0.01; // mm

            TaskFunction_ task_function = new TaskFunction_();
            task_function.setLogging(false);
            task_function.setFreqCutoff(nyquist);
            task_function.setTaskRadius(this.task_radius);
            task_function.setTaskContrast(this.task_contrast);
            task_function.disc_function();

            UnivariateFunction task_fn = task_function.getTaskFunction();

            if (this.show_task_plots) {
                // Plot the task function
                double x_min = -2.0 * this.task_radius;
                double x_max = 2.0 * this.task_radius;
                int n_elements = (int) ((x_max - x_min) / task_x_inc) + 1;

                double[] x_pos = new double[n_elements];
                double[] y_val = new double[n_elements];

                for (int i = 0; i < n_elements; i++) {
                    x_pos[i] = x_min + i * task_x_inc;
                    y_val[i] = task_fn.value(x_pos[i]);
                }
                Plot plot_task = new Plot("Task function", "Distance (mm)", "Pixel value");
                plot_task.add("line", x_pos, y_val);
                String task_labels = "Task function\t";
                plot_task.addLegend(task_labels);
                plot_task.show();
            }


            // Obtain the FFT of the task function, W_fn
            task_function.disc_function_fft();
            UnivariateFunction W_fn = task_function.getTaskFunctionFFT();

            if (this.show_task_plots) {
                // Plot the FFT of the task function
                double f_inc = task_function.getFreqRes();
                int n_felements = (int) Math.ceil(nyquist / f_inc) + 1;

                double[] f_pos = new double[n_felements];
                double[] fft_val = new double[n_felements];
                for (int i = 0; i < n_felements; i++) {
                    f_pos[i] = i * f_inc;
                    fft_val[i] = W_fn.value(f_pos[i]);
                }
                Plot plot_task_fft = new Plot("Task function fft", "Frequency (mm^-1)", "Value");
                plot_task_fft.add("line", f_pos, fft_val);
                String task_fft_labels = "FFT of task function up to Nyquist\t";
                plot_task_fft.addLegend(task_fft_labels);
                plot_task_fft.show();
            }


            // Calculate the detectability index using a modified version of equation 2 in http://dx.doi.org/10.1118/1.4923172
            // We are using radial functions here, so it's only a single integral over the radial frequency, rather than
            // a double integral over u and v.
            UnivariateInterpolator linear_interpolator = new LinearInterpolator();
            UnivariateFunction TTF_fn;
            UnivariateFunction E_fn;
            UnivariateFunction NPS_fn;
            try {
                TTF_fn = linear_interpolator.interpolate(ttf_result.freq, ttf_result.val);
                E_fn = linear_interpolator.interpolate(E.freq, E.val);
                NPS_fn = linear_interpolator.interpolate(nps_result.freq, nps_result.val);
            } catch (Exception e) {
                IJ.error("There was an error when interpolating the results:\n" + e);
                return;
            }

            // Using equation 2 in http://dx.doi.org/10.1118/1.4923172
            UnivariateFunction numerator = FunctionUtils.multiply(W_fn, W_fn, TTF_fn, TTF_fn, E_fn, E_fn);
            UnivariateFunction denominator = FunctionUtils.multiply(W_fn, W_fn, TTF_fn, TTF_fn, NPS_fn, E_fn, E_fn, E_fn, E_fn);


            // Now integrate the numerator and denominator
            UnivariateIntegrator simpson_integrator = new SimpsonIntegrator();
            double numerator_int = simpson_integrator.integrate(100000, numerator, 0, nyquist);
            double denominator_int = simpson_integrator.integrate(100000, denominator, 0, nyquist);
            double detectability = Math.sqrt(Math.pow(numerator_int, 2.0) / denominator_int);


            // Log the results
            IJ.log("Detectability results" +
                    "\nDetectability," + this.two_dp.format(detectability) +
                    "\nTask radius (mm)," + this.two_dp.format(this.task_radius) +
                    "\nContrast (HU)," + this.two_dp.format(cnr_result.contrast) +
                    "\nCNR," + this.two_dp.format(cnr_result.cnr) +
                    "\nNoise (HU)," + this.two_dp.format(cnr_result.noise) +
                    "\nTTF 50 (1/mm)," + this.two_dp.format(TTF.getTTF50()) +
                    "\nTTF 10 (1/mm)," + this.two_dp.format(TTF.getTTF10()) +
                    "\nTTF centre of mass," + this.two_dp.format(ttf_com[0]) + "," + this.two_dp.format(ttf_com[1])
            );

            IJ.log("\nEye model settings used" +
                    "\nViewing distance (mm)," + this.one_dp.format(eye_model.getViewingDistance()) +
                    "\nRecon field of view (mm)," + this.two_dp.format(eye_model.getReconFOV()) +
                    "\nImage pixel width," + this.three_dp.format(eye_model.getImagePixelWidth()) +
                    "\nDisplay monitor pixel pitch (mm)," + this.three_dp.format(eye_model.getDisplayPixelPitch()) +
                    "\nDisplay zoom factor," + this.one_dp.format(eye_model.getDisplayZoomFactor()) +
                    "\nSize of displayed image (mm)," + this.two_dp.format(eye_model.getDisplayedImageSize()) +
                    "\nRadial (1/mm) to angular frequency (c/deg) factor," + this.three_dp.format(eye_model.getFreqToAngFreqFactor())
            );


            // Log the TTF, NPS, Eye, and Task values
            IJ.log("\nDetectability components" +
                    "\nFreq (mm^-1),TTF,NPS,Eye,Task"
            );
            double f = 0.0;
            while (f <= nyquist) {
                IJ.log(this.three_dp.format(f) + "," +
                        this.three_dp.format(TTF_fn.value(f)) + "," +
                        this.three_dp.format(NPS_fn.value(f)) + "," +
                        this.three_dp.format(E_fn.value(f)) + "," +
                        this.three_dp.format(W_fn.value(f))
                );
                f += 0.01;
            }
            IJ.log("\n");
        }
    }
}