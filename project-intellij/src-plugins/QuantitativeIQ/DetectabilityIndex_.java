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
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;

import java.text.DecimalFormat;

public class DetectabilityIndex_ implements PlugInFilter {
    private ImagePlus imp;

    private DecimalFormat two_dp = new DecimalFormat("0.00");

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


        // Obtain the Nyquist frequency of the image
        Calibration cal = imp.getCalibration();
        double nyquist = 1.0 / (2.0 * cal.pixelWidth);


        // Obtain the eye model data
        EyeModel_ eye_model = new EyeModel_();
        eye_model.setup(two_dp.format(nyquist));
        eye_model.run();
        Eye_result eye_result = eye_model.getEyeResult();




        // TEST
        // Play about with integration - note, need to multiply NPS, TTF etc
        // together before integrating when working out detectability index.
        // Integrate the TTF and NPS from 0 to Nyquist
        UnivariateIntegrator simpson_integrator = new SimpsonIntegrator();
        UnivariateInterpolator linear_interpolator = new LinearInterpolator();

        UnivariateFunction ttf_function = linear_interpolator.interpolate(ttf_result.freq, ttf_result.val);
        double ttf_int = simpson_integrator.integrate(2000, ttf_function, 0, nyquist);

        UnivariateFunction nps_function = linear_interpolator.interpolate(nps_result.freq, nps_result.val);
        double nps_int = simpson_integrator.integrate(2000, nps_function, 0, nyquist);
        // TEST end
    }
}
