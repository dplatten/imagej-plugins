// ImageJ plugin to calculate the detectability index of a circular object
// present in an image.
//
//package ij.plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class DetectabilityIndex_ implements PlugInFilter {
    private ImagePlus imp;

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

        // First run the TTF analysis
        TaskTransferFunction_ ttf = new TaskTransferFunction_();
        ttf.setup("defaults", imp);
        ttf.run(ip);

        TaskTransferFunction_.TTF_result ttf_results = ttf.getTTFResult();
        TaskTransferFunction_.CNR_result cnr_results = ttf.getCNRResult();
    }
}
