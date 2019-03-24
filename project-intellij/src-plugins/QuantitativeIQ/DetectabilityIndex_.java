package QuantitativeIQ;// ImageJ plugin to calculate the detectability index of a circular object
// present in an image.
//
//package ij.plugin;

import QuantitativeIQ.TaskTransferFunction_;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import org.apache.commons.math3.analysis.UnivariateFunction;

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

        TaskFunction_ task_function = new TaskFunction_();
        task_function.setTaskRadius(1.5);
        task_function.setTaskContrast(30);

        task_function.disc_function();
        UnivariateFunction task_fn = task_function.getTaskFunction();
        double[] x_pos = new double[1001];
        double[] y_val = new double[1001];
        double x_min = -2.0 * 1.5;
        double x_max = 2.0 * 1.5;
        double x_inc = (x_max - x_min) / 1000.0;
        for (int x=0; x<1001; x++) {
            x_pos[x] = x_min + x*x_inc;
            y_val[x] = task_fn.value(x_pos[x]);
        }
        Plot plot_task = new Plot("Task function", "Distance (mm)", "Pixel value");
        plot_task.add("line", x_pos, y_val);
        String task_labels = "Task function\t";
        plot_task.addLegend(task_labels);
        plot_task.show();


        task_function.disc_function_fft();
        UnivariateFunction task_fn_fft =  task_function.getTaskFunctionTTF();
        double[] f_pos = new double[1001];
        double[] fft_val = new double[1001];
        double f_max = 5.0;
        double f_inc = f_max / 1000.0;
        for (int x=0; x<1001; x++) {
            f_pos[x] = x*f_inc;
            fft_val[x] = task_fn_fft.value(f_pos[x]);
        }
        Plot plot_task_fft = new Plot("Task function fft", "Frequency (mm^-1)", "Value");
        plot_task_fft.add("line", f_pos, fft_val);
        String task_fft_labels = "Task function fft\t";
        plot_task_fft.addLegend(task_fft_labels);
        plot_task_fft.show();
    }
}
