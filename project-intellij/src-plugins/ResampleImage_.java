// ImageJ plugin to resample an image to a lower number of counts.
//
// David Platten, 2019
//
// Multi-threading tips from here:
// https://albert.rierol.net/imagej_programming_tutorials.html#Multithreaded%20programming%20for%20ImageJ
import ij.*;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;

class MyRunnable implements Runnable {
    private float p;
    private int n_threads, thread_num, width, height, mask_value;
    private boolean rescale_values;
    private ImagePlus image;
    private ImageProcessor image_proc;
    private Calibration cal;

    public MyRunnable(int thread_num, int n_threads, int width, int height, ImagePlus image, Calibration cal, int mask_value, float p, boolean rescale_values) {
        this.thread_num = thread_num;
        this.n_threads = n_threads;
        this.width = width;
        this.height = height;
        this.image = image;
        this.image_proc = image.getProcessor();
        this.cal = cal;
        this.mask_value = mask_value;
        this.p = p;
    }

    public void run() {
        for (int y = thread_num; y < height; y+=n_threads) {
            for (int x = 0; x < width; x++) {
                int current_val = image_proc.getPixel(x, y);
                if (cal.getCValue(current_val) != mask_value) {
                    int count = 0;
                    for (int n = 0; n < current_val; n++) {
                        if (Math.random() < p) count++;
                    }

                    // Rescale count back to the original value if required
                    if (rescale_values) count /= p;

                    // Set this pixel location to be the value of count
                    image_proc.putPixel(x, y, (int) count);
                }
            }
            // Update the image each time a row is complete
            image.updateAndDraw();

            // Update the progress bar
            //float progress = (float) y / (float) height;
            //IJ.showStatus("Calculating row " + y + " of " + height);
            //IJ.showProgress(progress);
        }
    }
}

public class ResampleImage_ implements PlugInFilter {
    private ImagePlus imp;
    private Calibration cal;
    private int mask_value = -2048;
    private float p = 0.5f;
    private boolean rescale_values = false;

    public int setup(String arg, ImagePlus imp) {
    	this.imp = imp;
    	this.cal = imp.getCalibration();

        GenericDialog gd = new GenericDialog("Processing options");
        gd.addNumericField("Resampling factor (0.0 - 1.0", this.p, 2);
        gd.addNumericField("Enter mask value (pixels with this value will be ignored)", this.mask_value, 0);
        gd.addCheckbox("Re-scale results to maintain mean values?", this.rescale_values);

        gd.showDialog();
        if(gd.wasCanceled()) {
            IJ.error("Plugin cancelled");
            return(0);
        }

        this.p = (float) gd.getNextNumber();
        this.mask_value = (int) gd.getNextNumber();
        this.rescale_values = gd.getNextBoolean();

    	return DOES_ALL;
    }


    public void run(final ImageProcessor ip) {
    	// Obtain the width and height of the image.
    	final int width  = ip.getWidth();
    	final int height = ip.getHeight();

    	final float p = this.p;
    	final Calibration cal = this.cal;
        final int mask_value = this.mask_value;
        final boolean rescale_values = this.rescale_values;

        // Create a duplicate of the image
    	final ImagePlus img_copy = new Duplicator().run(this.imp);
        img_copy.show();

    	// Main loop that steps through each pixel in the duplicate image
        final ImageProcessor img_copy_ip = img_copy.getProcessor();


        final Thread[] threads = newThreadArray();
        for (int ithread=0; ithread<threads.length; ithread++) {
            threads[ithread] = new Thread(new MyRunnable(ithread, threads.length, width, height, img_copy, cal, mask_value, p, rescale_values));
        }
        startAndJoin(threads);
    }


    /** Create a Thread[] array as large as the number of processors available.
     * From Stephan Preibisch's Multithreading.java class. See:
     * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
     */
    private Thread[] newThreadArray() {
        int n_cpus = Runtime.getRuntime().availableProcessors();
        //int n_cpus = 8;
        return new Thread[n_cpus];
    }


    /** Start all given threads and wait on each of them until all are done.
     * From Stephan Preibisch's Multithreading.java class. See:
     * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
     */
    public static void startAndJoin(Thread[] threads)
    {
        for (int ithread = 0; ithread < threads.length; ++ithread) {
            threads[ithread].setPriority(Thread.NORM_PRIORITY);
            threads[ithread].start();
        }

        try {
            for (int ithread = 0; ithread < threads.length; ++ithread)
                threads[ithread].join();
        }
        catch (InterruptedException ie) {
            throw new RuntimeException(ie);
        }
    }
}