package QuantitativeIQ;
import ij.IJ;
import ij.gui.*;
import ij.plugin.PlugIn;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.MathArrays;

import java.text.DecimalFormat;

/**
 * ImageJ plugin to calculate eye frequency response using either
 * Richard and Siewerdsen 2008:
 * E(f) = f^n.exp(-c.f^2)
 * where n = 1.3
 *       c = 3.093 for a 50 cm viewing distance for E(f) to peak at 4 cycles/degree
 *
 * or Solomon et al 2015:
 * E(p) = |n.p^(a1).exp(-a2.p^a3)|^2
 * where a1 = 1.5
 *       a2 = 0.98
 *       a3 = 0.68
 *       n is a normalisation constant to make max E(p) = 1.0
 *       p (rho, really) is angular spatial frequency expressed in cycles per degree
 *
 * p can be calculated from radial spatial frequency using:
 * p = (f.FOV.R.pi) / (D.180)
 * where f is radial spatial frequency
 *       FOV is the reconstructed field of view of the image
 *       R is the viewing distance
 *       D is the display size (assumed to be 30.5 cm in the Solomon paper)
 *       Note to self: is D the physical size of the full image as displayed on the screen?
 *                     If so, can calculate this if know the pixel pitch in mm, the pixel
 *                     dimensions of the image and the zoom factor... this is what ImaQuest asks for
 */
public class EyeModel_ implements PlugIn {

    // Richard and Siewerdsen 2008, assuming a 50 cm viewing distance
    private double n = 1.3;
    private double c = 3.093; // So that E(f) peaks at 4 c/deg at 50 cm viewing distance

    // Solomon et al 2015
    private double recon_fov = 400.0; // Reconstructed field of view of image (mm)
    private double viewing_dist = 500.0; // The viewing distance (mm)
    private double display_pixel_pitch = 0.2; // The size of a pixel on the display (mm)
    private double image_pixel_width = 512; // The number of pixels in the image width
    private double display_zoom_factor = 3.0; // The zoom applied to the image on the display
    private double display_size = calc_displayed_image_size(); // The display size (mm)
    private double f_to_rho = calc_f_to_rho_factor(); // Conversion from radial freq to cycles per degree

    private double cutoff_f = 2.0; // Maximum frequency for which to calculate E(f)
    private double f_inc = 0.01; // Frequency increment in mm^-1

    private static Eye_result eye_result;

    private DecimalFormat three_dp = new DecimalFormat("0.000");

    private static int SOLOMON_ET_AL = 0;
    private static int RICHARD_AND_SIEWERDSEN = 1;
    private int model_choice = SOLOMON_ET_AL;


    /** Method to obtain the eye filter results
     *
     * @return An object of type Eye_result containing the eye filter results. The <i>freq</i> and <i>val</i> properties
     * each contain an array of double-precision values.
     */
    Eye_result getEyeResult() {
        return eye_result;
    }


    /** Method to set the image reconstruction field of view
     *
     * @param new_recon_fov The value of the image reconstruction field of view in mm.
     */
    void setReconFOV(double new_recon_fov) {
        this.recon_fov = new_recon_fov;
        this.f_to_rho = calc_f_to_rho_factor();
    }


    /** Method to set the pixel pitch of the display monitor being used to view the image
     *
     * @param new_display_pixel_pitch The physical size of single pixel of the display monitor in mm
     */
    void setDisplayPixelPitch(double new_display_pixel_pitch) {
        this.display_pixel_pitch = new_display_pixel_pitch;
        this.display_size = calc_displayed_image_size();
        this.f_to_rho = calc_f_to_rho_factor();
    }


    /** Method to set the image pixel size
     *
     * @param new_image_pixel_width The size of an image pixel (usually defined in the DICOM tags) n mm
     */
    void setImagePixelWidth(double new_image_pixel_width) {
        this.image_pixel_width = new_image_pixel_width;
        this.display_size = calc_displayed_image_size();
        this.f_to_rho = calc_f_to_rho_factor();
    }


    /** Method to set the zoom factor that is used to display the image on the display monitor
     *
     * @param new_display_zoom_factor The zoom factor of the image on the display monitor
     */
    void setDisplayZoomFactor(double new_display_zoom_factor) {
        this.display_zoom_factor = new_display_zoom_factor;
        this.display_size = calc_displayed_image_size();
        this.f_to_rho = calc_f_to_rho_factor();
    }


    /** Method to set up the eye filter according to the user's requirements
     *
     * @param arg If present, assumed to contain the frequency cut off as a String
     * @param model_type If present used to determine the eye model type to use
     * @return 0 if all is well, 1 if the routine is cancelled
     */
    int setup(String arg, int model_type) {
        if (model_type == SOLOMON_ET_AL) {
            IJ.log("Using Solomon et al model");
            this.model_choice = SOLOMON_ET_AL;
        }
        else if (model_type == RICHARD_AND_SIEWERDSEN) {
            IJ.log("Using Richard and Siewerdsen model");
            this.model_choice = RICHARD_AND_SIEWERDSEN;
        }
        else {
            GenericDialog gd = new GenericDialog("Choose which eye model to use");
            String[] choices = {"Solomon et al", "Richard and Siewerdsen"};
            gd.addChoice("Eye model", choices, choices[0]);
            gd.showDialog();
            if(gd.wasCanceled()) {
                return(1); // Return to the calling method with a non-zero result
            }
            this.model_choice = gd.getNextChoiceIndex();
        }

        int result;
        if (this.model_choice == SOLOMON_ET_AL) {
            result = SolomonEtAlSetup(arg);
        }
        else {
            result = RichardAndSiewerdsenSetup(arg);
        }
        if (result != 0) {
            return(1); // Return to the calling method with a non-zero result
        }
        else return(0);
    }


    /** Method to calculate the selected eye model with the set parameters
     *
     * @param arg If present then create the model; if missing then call the setup method then run
     */
    public void run(String arg) {
        if ((arg != null) && (arg.length() > 0)) {
            if (this.model_choice == SOLOMON_ET_AL) {
                SolomonEtAl();
            } else if (this.model_choice == RICHARD_AND_SIEWERDSEN) {
                RichardAndSiewerdsen();
            } else {
                IJ.error("Unrecognised value passed to run");
            }
        }
        else {
            int result  = setup("2.0", 999);
            if (result == 0) run(String.valueOf(this.model_choice));
        }
    }


    /** Method to setup the  Richard and Siewerdsen model
     *
     * @param arg Not currently used
     * @return 0 if all is well; 1 if the routine is cancelled
     */
    private int RichardAndSiewerdsenSetup(String arg) {
/*
        if ((arg != null) && (arg.length() > 0)) {
            if (arg.equals("defaults")) {
                // Use defaults
                IJ.log("Eye model using " + three_dp.format(cutoff_f) + " mm^-1 cutoff  frequency");
            }
            else {
                // Assume the argument contains the frequency cutoff.
                cutoff_f = Float.parseFloat(arg);
                IJ.log("Eye model using " + three_dp.format(cutoff_f) + " mm^-1 cutoff  frequency");
            }
        }
        else {
*/
        GenericDialog gd = new GenericDialog("Richard and Siewerdsen eye model options");
        gd.addMessage("This assumes a 500 mm viewing distance");
        gd.addNumericField("Frequency cutoff (mm^-1)", this.cutoff_f, 1);
        gd.addNumericField("Frequency increment (mm^-1)", this.f_inc, 2);
        gd.addNumericField("n", n, 3);
        gd.addNumericField("c", c, 3);

        gd.showDialog();
        if(gd.wasCanceled()) {
            return(1); // Return to the calling method with a non-zero result
        }

        this.cutoff_f = gd.getNextNumber();
        this.f_inc = gd.getNextNumber();
        this.n = gd.getNextNumber();
        this.c = gd.getNextNumber();
        //}
        return(0);
    }


    /** Method to calculate the eye model using Richard and Siewerdsen
     *
     */
    private void RichardAndSiewerdsen() {

        int num_elements = (int) Math.ceil(this.cutoff_f / this.f_inc) + 1;

        Eye_result E = new Eye_result();
        E.freq = new double[num_elements];
        E.val = new double[num_elements];

        for (int i=0; i<num_elements; i++) {
            E.freq[i] = this.f_inc * i;
            E.val[i] = Math.pow(E.freq[i], n) * Math.exp(-1.0 * c * E.freq[i] * E.freq[i]);
        }

        Plot plot_eye = new Plot("Richard and Siewerdsen eye response", "Frequency (mm^-1)", "Eye response");
        plot_eye.add("line", E.freq, E.val);
        String eye_labels = "Richard and Siewerdsen eye response\t";
        plot_eye.addLegend(eye_labels);
        plot_eye.show();

        // Write the results to the ImageJ log
        logResults(E, "Richard and Siewerdsen eye response");

        eye_result = E;
    }


    /** Method to setup the Solomon et al eye model
     *
     * @param arg Not currently used
     * @return 0 if all is well; 1 if the routine was cancelled
     */
    private int SolomonEtAlSetup(String arg) {
        // Entirely ignoring arg at the moment and showing the dialogue
        GenericDialog gd = new GenericDialog("Solomon et al eye model options");
        gd.addNumericField("Frequency cutoff (mm^-1)", this.cutoff_f, 2);
        gd.addNumericField("Frequency increment (mm^-1)", this.f_inc, 2);
        gd.addNumericField("Viewing distance (mm)", this.viewing_dist, 0);
        gd.addNumericField("Image reconstructed field of view (mm)", this.recon_fov, 0);
        gd.addMessage("Physical image size on the display monitor:");
        gd.addNumericField("Display pixel pitch (mm)", this.display_pixel_pitch, 2);
        gd.addNumericField("Display zoom factor", this.display_zoom_factor, 1);
        gd.addNumericField("Image pixel width", this.image_pixel_width, 0);

        gd.showDialog();
        if(gd.wasCanceled()) {
            return(1);  // Return to the calling method with a non-zero result
        }

        this.cutoff_f = gd.getNextNumber();
        this.f_inc = gd.getNextNumber();
        this.viewing_dist = gd.getNextNumber();
        this.display_pixel_pitch = gd.getNextNumber();
        this.display_zoom_factor = gd.getNextNumber();
        this.image_pixel_width = gd.getNextNumber();
        this.recon_fov = gd.getNextNumber();

        // Recalculate the conversion from radial frequency to angular frequency
        this.display_size = calc_displayed_image_size();
        this.f_to_rho = calc_f_to_rho_factor();

        return(0);
    }


    /** Method to calculate the eye model using Solomon et al
     *
     */
    private void SolomonEtAl() {
        double a1 = 1.5;
        double a2 = 0.98;
        double a3 = 0.68;

        int num_elements = (int) Math.ceil(this.cutoff_f / this.f_inc) + 1;

        Eye_result E = new Eye_result();
        E.freq = new double[num_elements];
        E.val = new double[num_elements];

        double[] ang_freq = new double[num_elements];

        for (int i=0; i<num_elements; i++) {
            E.freq[i] = this.f_inc * i; // The radial frequency
            ang_freq[i] = E.freq[i] * this.f_to_rho; // The angular frequency
            // E, pre-normalisation and squaring
            E.val[i] = Math.pow(ang_freq[i], a1) * Math.exp(-1.0 * a2 * Math.pow(ang_freq[i], a3));
        }

        // Determine the maximum value of E.val and then use this to normalise the result
        double max_val = StatUtils.max(E.val);
        MathArrays.scaleInPlace(1.0/max_val, E.val);

        // Square the normalised values to give us E
        E.val = MathArrays.ebeMultiply(E.val, E.val);

        Plot plot_eye = new Plot("Solomon et al eye response", "Frequency (mm^-1)", "Eye response");
        plot_eye.add("line", E.freq, E.val);
        String eye_labels = "Solomon et al eye response\t";
        plot_eye.addLegend(eye_labels);
        plot_eye.show();

        // Write the results to the ImageJ log
        logResults(E, "Solomon et al eye response");

        eye_result = E;
    }


    /** Method to write the eye model results to the ImageJ log
     *
     * @param E The eye model
     * @param label The label for the eye model values
     */
    private void logResults(Eye_result E, String label) {
        IJ.log("Eye model using " + this.three_dp.format(this.cutoff_f) + " mm^-1 cutoff  frequency");
        IJ.log("Freq (mm^-1)," + label);

        for (int i=0; i<E.freq.length; i++) {
            IJ.log("" + this.three_dp.format(E.freq[i]) + ", " + this.three_dp.format(E.val[i]));
        }
    }


    /** Method to calculate the conversion from radial frequency (per mm) to angular frequency (cycles per degree)
     *
     * @return The conversion factor
     */
    private double calc_f_to_rho_factor() {
        return (this.recon_fov * this.viewing_dist * Math.PI) / (this.display_size * 180.0);
    }


    /** Method to calculate the displayed image size
     *
     * @return The conversion factor
     */
    private double calc_displayed_image_size() {
        return this.display_pixel_pitch * this.image_pixel_width * this.display_zoom_factor;
    }


    static class Eye_result {
        double[] freq;
        double[] val;
    }
}
