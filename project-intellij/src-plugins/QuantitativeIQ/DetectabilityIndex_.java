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
import ij.measure.Calibration;
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
        double task_radius = 1.5; // mm
        double task_contrast = 30.0; // HU

        Calibration cal = imp.getCalibration();
        double pixel_size_in_mm = cal.pixelWidth;
        double nyquist = 1.0 / (2.0*pixel_size_in_mm);
        double task_inc = 0.01; // mm

        task_function.setTaskRadius(task_radius);
        task_function.setTaskContrast(task_contrast);

        task_function.disc_function();
        UnivariateFunction task_fn = task_function.getTaskFunction();
        double x_min = -2.0 * task_radius;
        double x_max = 2.0 * task_radius;
        int n_elements = (int) ((x_max - x_min) / task_inc) + 1;

        double[] x_pos = new double[n_elements];
        double[] y_val = new double[n_elements];

        for (int i=0; i<n_elements; i++) {
            x_pos[i] = x_min + i*task_inc;
            y_val[i] = task_fn.value(x_pos[i]);
        }
        Plot plot_task = new Plot("Task function", "Distance (mm)", "Pixel value");
        plot_task.add("line", x_pos, y_val);
        String task_labels = "Task function\t";
        plot_task.addLegend(task_labels);
        plot_task.show();


        task_function.disc_function_fft();
        UnivariateFunction task_fn_fft =  task_function.getTaskFunctionTTF();
        double f_inc = task_function.getFreqRes();
        int n_felements = (int) Math.ceil(nyquist / f_inc) + 1;

        double[] f_pos = new double[n_felements];
        double[] fft_val = new double[n_felements];
        for (int i=0; i<n_felements; i++) {
            f_pos[i] = i*f_inc;
            fft_val[i] = task_fn_fft.value(f_pos[i]);
        }
        Plot plot_task_fft = new Plot("Task function fft", "Frequency (mm^-1)", "Value");
        plot_task_fft.add("line", f_pos, fft_val);
        String task_fft_labels = "FFT of task function up to Nyquist\t";
        plot_task_fft.addLegend(task_fft_labels);
        plot_task_fft.show();
    }
}
