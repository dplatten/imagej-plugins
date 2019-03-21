package QuantitativeIQ;
// ImageJ plugin to calculate eye frequency response using
// Richard and Siewerdsen 2008
//
// E(f) = f^n * exp(-c.f^2)
//
//package ij.plugin;

import ij.IJ;
import ij.gui.*;
import java.text.DecimalFormat;

public class EyeModel_ {

    private double n = 1.3;
    private double c = 3.093; // So that E(f) peaks at 4 c/deg at 50 cm viewing distance
    private double cutoff_f = 2.0; // Maximum frequency for which to calculate E(f)
    private double f_inc = 0.01; // Frequency increment in mm^-1

    private static Eye_result eye_result;

    private DecimalFormat three_dp = new DecimalFormat("0.000");

    // Defined to enable other plugins to access the results
    Eye_result getEyeResult() {
        return eye_result;
    }

    public int setup(String arg) {
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
            GenericDialog gd = new GenericDialog("Eye model options");
            gd.addNumericField("Frequency cutoff (mm^-1)", cutoff_f, 1);
            gd.addNumericField("Frequency increment (mm^-1)", f_inc, 2);
            gd.addNumericField("n", n, 3);
            gd.addNumericField("c", c, 3);

            gd.showDialog();
            if(gd.wasCanceled()) {
                IJ.error("Plugin cancelled");
                return(0);
            }

            cutoff_f = gd.getNextNumber();
            f_inc = gd.getNextNumber();
            n = gd.getNextNumber();
            c = gd.getNextNumber();
        }
        return(1);
    }

    public void run() {

        int num_elements = (int) Math.ceil(cutoff_f / f_inc) + 1;

        Eye_result E = new Eye_result();
        E.freq = new double[num_elements];
        E.val = new double[num_elements];

        for (int i=0; i<num_elements; i++) {
            E.freq[i] = f_inc * i;
            E.val[i] = Math.pow(E.freq[i], n) * Math.exp(-1.0 * c * E.freq[i] * E.freq[i]);
        }

        // Write the results to the ImageJ log
        logResults(E, "NPS (HU^2.mm^2)");

        eye_result = E;
    }


    private void logResults(Eye_result E, String label) {
        IJ.log("Freq (mm^-1)," + label);

        for (int i=0; i<E.freq.length; i++) {
            IJ.log("" + three_dp.format(E.freq[i]) + ", " + three_dp.format(E.val[i]));
        }
    }


    static class Eye_result {
        double[] freq;
        double[] val;
    }

}
