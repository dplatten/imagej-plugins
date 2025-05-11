

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

import ij.IJ;
import ij.gui.*;
import ij.plugin.PlugIn;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.MathArrays;

import java.text.DecimalFormat;

/**
 *
 * ImageJ plugin to calculate eye frequency response using either
 * Richard and Siewerdsen 2008, <a href="http://dx.doi.org/10.1118/1.2988161">http://dx.doi.org/10.1118/1.2988161</a>:
 * <br><br>
 * E(f) = f<sup>n</sup>.exp(-c.f<sup>2</sup>)
 * <pre>
 *
 * where f = radial spatial frequency (mm<sup>-1</sup>)
 *       n = 1.3
 *       c = 3.093 for a 50 cm viewing distance for E(f) to peak at 4 degree<sup>-1</sup>
 * </pre>
 * or Solomon et al 2015, <a href="http://dx.doi.org/10.1118/1.4923172">http://dx.doi.org/10.1118/1.4923172</a>:
 * <br><br>
 * E(&rho;) = |n.&rho;<sup>a<sub>1</sub></sup>.exp(-a<sub>2</sub>.&rho;<sup>a<sub>3</sub></sup>)|<sup>2</sup>
 * <pre>
 *
 * where a<sub>1</sub> = 1.5
 *       a<sub>2</sub> = 0.98
 *       a<sub>3</sub> = 0.68
 *       n is a normalisation constant to make max E(&rho;) = 1.0
 *       &rho; is angular spatial frequency (degree<sup>-1</sup>)
 *
 * &rho; can be calculated from radial spatial frequency using:
 *
 * &rho; = (f.FOV.R.&pi;) / (D.180)
 *
 * where f is radial spatial frequency (mm<sup>-1</sup>)
 *       FOV is the reconstructed field of view of the image (mm)
 *       R is the viewing distance (mm)
 *       D is the display size (mm) (assumed to be 305 mm in the Solomon paper)
 *       Note to self: is D the physical size of the full image as displayed on the screen?
 *                     If so, can calculate this if know the pixel pitch in mm, the pixel
 *                     dimensions of the image and the zoom factor... this is what ImaQuest asks for
 * </pre>
 * <br><br>
 * The <i>Solomon et al 2015</i> option is mathematically the same as that used
 * by the Imquest software <i>Eckstein</i> option. The original Eckstein work
 * did not normalise the function, nor did it square the result
 * (<a href="http://dx.doi.org/10.1364/OE.11.000460">http://dx.doi.org/10.1364/OE.11.000460</a>).
 */

public class EyeModel_ implements PlugIn {

    // Richard and Siewerdsen 2008, http://dx.doi.org/10.1118/1.2988161, assuming a 50 cm viewing distance
    private double n = 1.3;
    private double c = 3.093; // So that E(f) peaks at 4 c/deg at 50 cm viewing distance

    // Solomon et al 2015, http://dx.doi.org/10.1118/1.4923172
    private double recon_fov = 400.0; // Reconstructed field of view of image (mm)
    private double viewing_dist = 500.0; // The viewing distance (mm)
    private double display_pixel_pitch = 0.2; // The size of a pixel on the display (mm)
    private double image_pixel_width = 512; // The number of pixels in the image width
    private double display_zoom_factor = 3.0; // The zoom applied to the image on the display

    private double display_size = this.calc_displayed_image_size(); // The size of the image on the display (mm)
    private double f_to_rho = this.calc_f_to_rho_factor(); // Conversion from radial freq to cycles per degree

    private double cutoff_f = 2.0; // Maximum frequency for which to calculate E(f)
    private double f_inc = 0.01; // Frequency increment in mm^-1

    private static Eye_result eye_result;

    private boolean log_results = true;
    private boolean show_plots = true;

    private DecimalFormat zero_dp = new DecimalFormat("0");
    private DecimalFormat one_dp = new DecimalFormat("0.0");
    private DecimalFormat two_dp = new DecimalFormat("0.00");
    private DecimalFormat three_dp = new DecimalFormat("0.000");

    static final int SOLOMON_ET_AL = 0; // Don't make private as external methods may wish to access these values
    static final int RICHARD_AND_SIEWERDSEN = 1; // Don't make private as external methods may wish to access these values
    private int model_choice = SOLOMON_ET_AL;
    private static String[] model_descriptions = {
            "Solomon et al 2015,http://dx.doi.org/10.1118/1.4923172",
            "Richard and Siewerdsen 2008,http://dx.doi.org/10.1118/1.2988161"};
    private static String[] model_labels = {
            "Solomon et al 2015",
            "Richard and Siewerdsen 2008"};


    /** Method to obtain a text description of the eye model used
     *
     * @return A brief text string with the eye model used
     */
    String getModelDescription() { return model_descriptions[this.model_choice]; }


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


    /** Method to obtain the eye filter results
     *
     * @return An object of type Eye_result containing the eye filter results. The <i>freq</i> and <i>val</i> properties
     * each contain an array of double-precision values.
     */
    Eye_result getEyeResult() {
        return eye_result;
    }


    /** Method to obtain the viewing distance
     *
     * @return The viewing distance.
     */
    double getViewingDistance() {
        return this.viewing_dist;
    }


    /** Method to set the frequency cut off
     *
     * @param new_freq The new frequency cut off (mm<sup>-1</sup>)
     */
    void setFreqCutoff(double new_freq) {
        this.cutoff_f = new_freq;
    }


    /** Method to set the image reconstruction field of view
     *
     * @param new_recon_fov The value of the image reconstruction field of view (mm).
     */
    void setReconFOV(double new_recon_fov) {
        this.recon_fov = new_recon_fov;
        this.f_to_rho = this.calc_f_to_rho_factor();
    }


    /** Method to get the image reconstruction field of view
     *
     * @return Returns the image reconstruction field of view.
     */
    double getReconFOV() {
        return this.recon_fov;
    }


    /** Method to set the pixel pitch of the display monitor being used to view the image
     *
     * @param new_display_pixel_pitch The physical size of single pixel of the display monitor (mm)
     */
    void setDisplayPixelPitch(double new_display_pixel_pitch) {
        this.display_pixel_pitch = new_display_pixel_pitch;
        this.display_size = this.calc_displayed_image_size();
        this.f_to_rho = this.calc_f_to_rho_factor();
    }


    /** Method to get the pixel pitch of the display monitor being used to view the image
     *
     * @return Returns the physical size of single pixel of the display monitor (mm)
     */
    double getDisplayPixelPitch() {
        return this.display_pixel_pitch;
    }


    /** Method to get the physical size of the image on the display screen
     *
     * @return Returns the physical size of the image on the display screen (mm)
     */
    double getDisplayedImageSize() {
        return this.display_size;
    }


    /** Method to get the conversion factor from radial frequency to angular frequency.
     *
     * @return Returns the conversion factor to go from radial (mm<sup>-1</sup>) to angular frequency (degree<sup>-1</sup>)
     */
    double getFreqToAngFreqFactor() {
        return this.f_to_rho;
    }


    /** Method to set the number of pixels in the image width
     *
     * @param new_image_pixel_width The number of pixels in the image width
     */
    void setImagePixelWidth(double new_image_pixel_width) {
        this.image_pixel_width = new_image_pixel_width;
        this.display_size = this.calc_displayed_image_size();
        this.f_to_rho = this.calc_f_to_rho_factor();
    }


    /** Method to set the number of pixels in the image width
     *
     * @return The number of pixels in the image width
     */
    double getImagePixelWidth() {
        return this.image_pixel_width;
    }


    /** Method to set the zoom factor that is used to display the image on the display monitor
     *
     * @param new_display_zoom_factor The zoom factor of the image on the display monitor
     */
    void setDisplayZoomFactor(double new_display_zoom_factor) {
        this.display_zoom_factor = new_display_zoom_factor;
        this.display_size = this.calc_displayed_image_size();
        this.f_to_rho = this.calc_f_to_rho_factor();
    }


    /** Method to get the zoom factor that is used to display the image on the display monitor
     *
     * @return Returns the zoom factor of the image on the display monitor
     */
    double getDisplayZoomFactor() {
        return this.display_zoom_factor;
    }


    /** Method to set up the eye filter according to the user's requirements
     *
     * @param model_type If present used to determine the eye model type to use
     * @return 0 if all is well, 1 if the routine is cancelled
     */
    int setup(int model_type) {

        if (model_type == SOLOMON_ET_AL) {
            this.model_choice = SOLOMON_ET_AL;
        }
        else if (model_type == RICHARD_AND_SIEWERDSEN) {
            this.model_choice = RICHARD_AND_SIEWERDSEN;
        }
        else {
            GenericDialog gd = new GenericDialog("Choose which eye model to use");
            gd.addChoice("Eye model", model_labels, model_labels[SOLOMON_ET_AL]);
            gd.showDialog();
            if(gd.wasCanceled()) {
                return(1); // Return to the calling method with a non-zero result
            }
            this.model_choice = gd.getNextChoiceIndex();
        }

        int result;
        if (this.model_choice == SOLOMON_ET_AL) {
            result = this.SolomonEtAlSetup();
        }
        else {
            result = this.RichardAndSiewerdsenSetup();
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
        // Requires at least ImageJ version 1.49v as this is when the ij.gui.Plot
        // class was introduced:
        // https://imagej.nih.gov/ij/developer/api/ij/gui/Plot.html
        if (IJ.versionLessThan("1.49v")) return;

        if ((arg != null) && (arg.length() > 0)) {
            if (this.model_choice == SOLOMON_ET_AL) {
                this.SolomonEtAl();
            } else if (this.model_choice == RICHARD_AND_SIEWERDSEN) {
                this.RichardAndSiewerdsen();
            } else {
                IJ.error("Unrecognised value passed to run");
            }
        }
        else {
            int result = this.setup(999);
            if (result == 0) this.run(String.valueOf(this.model_choice));
        }
    }


    /** Method to setup the <a href="http://dx.doi.org/10.1118/1.2988161">Richard and Siewerdsen</a> model
     *
     * @return 0 if all is well; 1 if the routine is cancelled
     */
    private int RichardAndSiewerdsenSetup() {
/*
        if ((arg != null) && (arg.length() > 0)) {
            if (arg.equals("defaults")) {
                // Use defaults
                if (this.log_results) IJ.log("Eye model using " + three_dp.format(cutoff_f) + " mm^-1 cutoff  frequency");
            }
            else {
                // Assume the argument contains the frequency cutoff.
                cutoff_f = Float.parseFloat(arg);
                if (this.log_results) IJ.log("Eye model using " + three_dp.format(cutoff_f) + " mm^-1 cutoff  frequency");
            }
        }
        else {
*/
        GenericDialog gd = new GenericDialog(model_labels[this.model_choice] + " options");
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


    /** Method to calculate the eye model using <a href="http://dx.doi.org/10.1118/1.2988161">Richard and Siewerdsen</a>
     *
     */
    private void RichardAndSiewerdsen() {

        int num_elements = (int) Math.ceil(this.cutoff_f / this.f_inc) + 1;

        Eye_result E = new Eye_result();
        E.freq = new double[num_elements];
        E.val = new double[num_elements];

        for (int i=0; i<num_elements; i++) {
            E.freq[i] = this.f_inc * i;
            E.val[i] = Math.pow(E.freq[i], this.n) * Math.exp(-1.0 * this.c * E.freq[i] * E.freq[i]);
        }

        if (this.show_plots) {
            Plot plot_eye = new Plot(model_labels[this.model_choice], "Radial spatial frequency (mm^-1)", "Eye response");
            plot_eye.add("line", E.freq, E.val);
            String eye_labels = model_labels[this.model_choice] + "\t";
            plot_eye.addLegend(eye_labels);
            plot_eye.show();
        }

        // Write the results to the ImageJ log
        if (this.log_results) logResults(E, model_labels[this.model_choice]);

        eye_result = E;
    }


    /** Method to setup the Solomon et al eye model (<a href="http://dx.doi.org/10.1118/1.4923172">http://dx.doi.org/10.1118/1.4923172</a>)
     *
     * @return 0 if all is well; 1 if the routine was cancelled
     */
    private int SolomonEtAlSetup() {
        // Entirely ignoring arg at the moment and showing the dialogue

        GenericDialog gd = new GenericDialog(model_labels[this.model_choice] + " options");
        gd.addNumericField("Radial spatial frequency cutoff (mm^-1)", this.cutoff_f, 2);
        gd.addNumericField("Radial spatial frequency increment (mm^-1)", this.f_inc, 2);
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
        this.recon_fov = gd.getNextNumber();
        this.display_pixel_pitch = gd.getNextNumber();
        this.display_zoom_factor = gd.getNextNumber();
        this.image_pixel_width = gd.getNextNumber();

        // Recalculate the conversion from radial frequency to angular frequency
        this.display_size = this.calc_displayed_image_size();
        this.f_to_rho = this.calc_f_to_rho_factor();

        return(0);
    }


    /** Method to calculate the eye model using Solomon et al 2015. This
     * method implements equation 3 in the <a href="http://dx.doi.org/10.1118/1.4923172">Solomon et al 2015</a> work.
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
            E.freq[i] = this.f_inc * i;              // The radial frequency
            ang_freq[i] = E.freq[i] * this.f_to_rho; // Convert to angular frequency

            // Calculate the inner part of Solomon et al 2015 equation 3, except for normalisation
            E.val[i] = Math.pow(ang_freq[i], a1) * Math.exp(-1.0 * a2 * Math.pow(ang_freq[i], a3));
        }

        // Determine the maximum value of E.val and then use this to normalise the result. This
        // completes the inner part of Solomon et al 2015 equation 3
        double max_val = StatUtils.max(E.val);
        MathArrays.scaleInPlace(1.0/max_val, E.val);

        // Square the normalised values to give us E. I don't see any point in taking the absolute
        // value first, as the squaring will render any negative values positive anyway.
        E.val = MathArrays.ebeMultiply(E.val, E.val);

        if (this.show_plots) {
            Plot plot_eye = new Plot(model_labels[this.model_choice], "Radial spatial frequency (mm^-1)", "Relative eye response");
            plot_eye.add("line", E.freq, E.val);
            plot_eye.addLabel(0.5, 0.2,  model_labels[this.model_choice]+"\n" +
                    this.zero_dp.format(viewing_dist) + " mm viewing distance\n" +
                    this.zero_dp.format(recon_fov) + " mm recon field of view\n" +
                    this.zero_dp.format(image_pixel_width) + " image pixel width\n" +
                    this.two_dp.format(display_pixel_pitch) + " mm display monitor pixel pitch\n" +
                    this.one_dp.format(display_zoom_factor) + " display zoom factor\n" +
                    this.one_dp.format(display_size) + " mm displayed image size");
            plot_eye.show();
        }

        // Write the results to the ImageJ log
        if (this.log_results) logResults(E, model_labels[this.model_choice]);

        eye_result = E;
    }


    /** Method to write the eye model results to the ImageJ log
     *
     * @param E The eye model
     * @param label The label for the eye model values
     */
    private void logResults(Eye_result E, String label) {
        IJ.log("Eye model," + getModelDescription());

        if (this.model_choice == RICHARD_AND_SIEWERDSEN) {
            IJ.log("Viewing distance (mm),Assumed 500");
            IJ.log("Radial spatial frequency cutoff (mm^-1)," + this.three_dp.format(this.cutoff_f));
            IJ.log("n," + this.three_dp.format(this.n));
            IJ.log("c," + this.three_dp.format(this.c));
            IJ.log("Radial spatial frequency (mm^-1)," + label);
            for (int i = 0; i < E.freq.length; i++) {
                IJ.log("" + this.three_dp.format(E.freq[i]) + ", " + this.three_dp.format(E.val[i]));
            }
            IJ.log("\n");
        }
        else {
            IJ.log("Viewing distance (mm)," + this.viewing_dist);
            IJ.log("Radial spatial frequency cutoff (mm^-1)," + this.three_dp.format(this.cutoff_f));
            IJ.log("Image reconstructed field of view (mm)," + this.recon_fov);
            IJ.log("Display pixel pitch (mm)," + this.display_pixel_pitch);
            IJ.log("Display zoom factor," + this.display_zoom_factor);
            IJ.log("Image pixel width (mm)," + this.image_pixel_width);
            IJ.log("Conversion from radial to angular frequency," + this.f_to_rho);
            IJ.log("Angular frequency (deg^-1),Radial spatial frequency (mm^-1)," + label);
            for (int i = 0; i < E.freq.length; i++) {
                IJ.log("" + this.three_dp.format(E.freq[i]*this.f_to_rho) + "," + this.three_dp.format(E.freq[i]) + "," + this.three_dp.format(E.val[i]));
            }
            IJ.log("\n");
        }
    }


    /** Method to calculate the conversion from radial frequency (per mm) to angular frequency (cycles per degree).
     * This method uses equation 4 in <a href="http://dx.doi.org/10.1118/1.4923172">http://dx.doi.org/10.1118/1.4923172</a> to calculate this factor:
     * <br><br>
     * &rho; = r.FOV.R.&pi; / (D.180)
     *
     * <pre>
     * where:
     *     &rho;   is angular spatial frequency (degree<sup>-1</sup>)
     *     r   is radial spatial frequency in (mm<sup>-1</sup>)
     *     FOV is the reconstructed field of view of the image (mm)
     *     R   is the viewing distance (mm)
     *     D   is the display size (mm) (I'm assuming this is the physical size of the image on the display)
     * </pre>
     * This method calculates FOV.R.&pi; / (D.180), so that values of radial frequency can be multiplied by this
     * factor to obtain values of angular spatial frequency.
     *
     * @return The conversion factor
     */
    private double calc_f_to_rho_factor() {
        return (this.recon_fov * this.viewing_dist * Math.PI) / (this.display_size * 180.0);
    }


    /** Method to calculate the displayed image size from the display pixel pitch, the pixel
     * width of the image and the display zoom factor.
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
