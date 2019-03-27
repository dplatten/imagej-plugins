package QuantitativeIQ;
// ImageJ plugin to calculate eye frequency response using
// Richard and Siewerdsen 2008:
// E(f) = f^n.exp(-c.f^2)
// where n = 1.3
//       c = 3.093 for a 50 cm viewing distance for E(f) to peak at 4 cycles/degree
//
// or Solomon et al 2015:
// E(p) = |n.p^(a1).exp(-a2.p^a3)|^2
// where a1 = 1.5
//       a2 = 0.98
//       a3 = 0.68
//       n is a normalisation constant to make max E(p) = 1.0
//       p (rho, really) is angular spatial frequency expressed in cycles per degree
//
// p can be calculated from radial spatial frequency using:
// p = (f.FOV.R.pi) / (D.180)
// where f is radial spatial frequency
//       FOV is the reconstructed field of view of the image
//       R is the viewing distance
//       D is the display size (assumed to be 30.5 cm in the Solomon paper)
//       Note to self: is D the physical size of the full image as displayed on the screen?
//                     If so, can calculate this if know the pixel pitch in mm, the pixel
//                     dimensions of the image and the zoom factor... this is what ImaQuest asks for
//
import ij.IJ;
import ij.gui.*;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.MathArrays;

import java.text.DecimalFormat;

public class EyeModel_ {

    // Richard and Siewerdsen 2008, assuming a 50 cm viewing distance
    private double n = 1.3;
    private double c = 3.093; // So that E(f) peaks at 4 c/deg at 50 cm viewing distance

    // Solomon et al 2015
    private double recon_fov = 400.0; // Reconstructed field of view of image (mm)
    private double viewing_dist = 500.0; // The viewing distance (mm)
    private double display_size = 350.0; // The display size (mm)
    private double f_to_rho = calc_f_to_rho_factor(); // Conversion from radial freq to cycles per degree

    private double cutoff_f = 2.0; // Maximum frequency for which to calculate E(f)
    private double f_inc = 0.01; // Frequency increment in mm^-1

    private static Eye_result eye_result;

    private DecimalFormat three_dp = new DecimalFormat("0.000");

    static int SOLOMON_ET_AL = 0;
    static int RICHARD_AND_SIEWERDSEN = 1;
    private int model_choice = SOLOMON_ET_AL;

    // Defined to enable other plugins to access the results
    Eye_result getEyeResult() {
        return eye_result;
    }


    public int setup(String arg, int model_type) {
        if (model_type == SOLOMON_ET_AL) {
            IJ.log("Using Solomon et al model");
            model_choice = SOLOMON_ET_AL;
        }
        else if (model_type == RICHARD_AND_SIEWERDSEN) {
            IJ.log("Using Richard and Siewerdsen model");
            model_choice = RICHARD_AND_SIEWERDSEN;
        }
        else {
            GenericDialog gd = new GenericDialog("Choose which eye model to use");
            String[] choices = {"Solomon et al", "Richard and Siewerdsen"};
            gd.addChoice("Eye model", choices, choices[0]);
            gd.showDialog();
            if(gd.wasCanceled()) {
                IJ.error("Cancelled");
                return(1);
            }
            model_choice = gd.getNextChoiceIndex();
        }

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

        if (model_choice == SOLOMON_ET_AL) SolomonEtAlSetup(arg);
        else RichardAndSiewerdsenSetup(arg);

        return(0);
    }


    public int run() {
        if (model_choice == SOLOMON_ET_AL) {
            SolomonEtAl();
        }
        else if (model_choice == RICHARD_AND_SIEWERDSEN) {
            RichardAndSiewerdsen();
        }
        else {
            IJ.error("Unrecognised value passed to run");
            return(1);
        }

        return(0);
    }


    private int RichardAndSiewerdsenSetup(String arg) {
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
            GenericDialog gd = new GenericDialog("Richard and Siewerdsen eye model options");
            gd.addNumericField("Frequency cutoff (mm^-1)", cutoff_f, 1);
            gd.addNumericField("Frequency increment (mm^-1)", f_inc, 2);
            gd.addNumericField("n", n, 3);
            gd.addNumericField("c", c, 3);

            gd.showDialog();
            if(gd.wasCanceled()) {
                IJ.error("Plugin cancelled");
                return(1);
            }

            cutoff_f = gd.getNextNumber();
            f_inc = gd.getNextNumber();
            n = gd.getNextNumber();
            c = gd.getNextNumber();
        }
        return(0);
    }


    private void RichardAndSiewerdsen() {

        int num_elements = (int) Math.ceil(cutoff_f / f_inc) + 1;

        Eye_result E = new Eye_result();
        E.freq = new double[num_elements];
        E.val = new double[num_elements];

        for (int i=0; i<num_elements; i++) {
            E.freq[i] = f_inc * i;
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


    public int SolomonEtAlSetup(String arg) {
        // Entirely ignoring arg at the moment and showing the dialogue
        GenericDialog gd = new GenericDialog("Solomon et al eye model options");
        gd.addNumericField("Frequency cutoff (mm^-1)", cutoff_f, 2);
        gd.addNumericField("Frequency increment (mm^-1)", f_inc, 2);
        gd.addNumericField("Viewing distance (mm)", viewing_dist, 0);
        gd.addNumericField("Displayed image size (mm)", display_size, 0);
        gd.addNumericField("Recon field of view (mm)", recon_fov, 0);

        gd.showDialog();
        if(gd.wasCanceled()) {
            IJ.error("Plugin cancelled");
            return(1);
        }

        cutoff_f = gd.getNextNumber();
        f_inc = gd.getNextNumber();
        viewing_dist = gd.getNextNumber();
        display_size = gd.getNextNumber();
        recon_fov = gd.getNextNumber();

        // Recalculate the conversion from radial frequency to angular frequency
        f_to_rho = calc_f_to_rho_factor();

        return(0);
    }


    public void SolomonEtAl() {
        double a1 = 1.5;
        double a2 = 0.98;
        double a3 = 0.68;

        int num_elements = (int) Math.ceil(cutoff_f / f_inc) + 1;

        Eye_result E = new Eye_result();
        E.freq = new double[num_elements];
        E.val = new double[num_elements];

        double[] ang_freq = new double[num_elements];

        for (int i=0; i<num_elements; i++) {
            E.freq[i] = f_inc * i; // The radial frequency
            ang_freq[i] = E.freq[i] * f_to_rho; // The angular frequency
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


    private void logResults(Eye_result E, String label) {
        IJ.log("Freq (mm^-1)," + label);

        for (int i=0; i<E.freq.length; i++) {
            IJ.log("" + three_dp.format(E.freq[i]) + ", " + three_dp.format(E.val[i]));
        }
    }


    private double calc_f_to_rho_factor() {
        return (recon_fov * viewing_dist * Math.PI) / (display_size * 180.0);
    }


    static class Eye_result {
        double[] freq;
        double[] val;
    }
}
