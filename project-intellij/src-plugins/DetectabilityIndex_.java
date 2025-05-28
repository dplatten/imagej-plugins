

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

// ImageJ plugin to calculate the detectability index of a circular object present in an image.
//
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import ij.io.DirectoryChooser;
import org.apache.commons.math3.analysis.FunctionUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;

import java.text.DecimalFormat;
import ij.gui.Plot;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.TransformType;

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
    //private Roi whole_phantom_roi;

    private int eye_model = EyeModel_.SOLOMON_ET_AL;
    private boolean show_eye_plots = false;
    private boolean show_task_plots = false;
    private boolean show_nps_plots = false;
    private boolean show_ttf_plots = true; // It's a good idea to review these to ensure the results are OK

    private String results_filename = "results.csv";

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
            return(DONE);
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

        if (this.roi_manager == null) {
            this.roi_manager = new RoiManager();
            IJ.error("Please add at least one ROI to the ROI manager");
            return (DONE);
        }

        if (this.roi_manager.getCount() == 0) {
            IJ.error("Please add at least one ROI to the ROI manager");
            return (DONE);
        }

        return DOES_ALL;
    }


    /** Method to run the plugin
     *
     * @param ip The ImageProcessor to use
     */
    public void run(ImageProcessor ip) {
        // Requires at least ImageJ version 1.52b as this is when the isStack()
        // method of the ImagePlus class was introduced:
        // https://imagej.nih.gov/ij/developer/api/ij/ImagePlus.html#isStack--
        if (IJ.versionLessThan("1.52b")) return;


        // Calculate the detectability index according to equation 2 in http://dx.doi.org/10.1118/1.4923172,
        // adapted for radial frequency rather than u, v.
        //
        // Detectability index needs:
        //  Task function (W)
        //  Task Transfer Function (TTF)
        //  Noise Power Spectrum (NPS)
        //  Eye model (E)


        // Run the NPS on the image
        // If the user has set an ROI on the image then the centre of this is used as the centre
        // of the phantom. If no ROI is present then the centre of the image use used.
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
        //eye_model.setup(EyeModel_.SOLOMON_ET_AL);
        eye_model.setup(999);
        eye_model.run("normal");
        EyeModel_.Eye_result E = eye_model.getEyeResult();


        // Ask the user what task size they want to use
        GenericDialog gd = new GenericDialog("Task function settings");
        gd.addNumericField("Task diameter (mm)", this.task_radius * 2.0, 2);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        } else {
            this.task_radius = gd.getNextNumber() / 2.0;
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
            if (!TTF.getTTFsuccess()) {
                IJ.error("There was a problem calculating the TTF.");
                return;
            }
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
            double numerator_int, denominator_int;
            try {
                numerator_int = simpson_integrator.integrate(100000, numerator, 0, nyquist);
            }
            catch (Exception e) {
                IJ.error("There was an error when integrating the numerator\nfor the detectability index calculation:\n" + e);
                return;
            }
            try {
                denominator_int = simpson_integrator.integrate(100000, denominator, 0, nyquist);
            }
            catch (Exception e) {
                IJ.error("There was an error when integrating the denominator\nfor the detectability index calculation:\n" + e);
                return;
            }

            double detectability = Math.sqrt(Math.pow(numerator_int, 2.0) / denominator_int);


            // Create an image of the task including the effect of the TTF and NPS
            double nps_record_stdev = NPS.getRecordStdevSum() / NPS.getNumRecords();
            if (NPS.getSubstackStatus()) nps_record_stdev /= (NPS.getMaxStackIndex() - NPS.getMinStackIndex() + 1);

            calcTaskImageIncTTFandNPS(this.task_radius*2.0, TTF_fn, NPS_fn, this.imp.getWidth(), cal.pixelWidth, nps_record_stdev, detectability);
            //calcTaskImageIncTTFandNPS(this.task_radius*2.0, TTF_fn, NPS_fn, this.imp.getWidth(), cal.pixelWidth, 12.0, detectability);


            // Log the results
            IJ.log("Detectability results" +
                    "\nReconstruction filter," + TTF.getReconFilter() +
                    "\nDetectability," + this.two_dp.format(detectability) +
                    "\nTask diameter (mm)," + this.two_dp.format(this.task_radius * 2.0) +
                    "\nContrast (HU)," + this.two_dp.format(cnr_result.contrast) +
                    "\nCNR," + this.two_dp.format(cnr_result.cnr) +
                    "\nNoise (HU)," + this.two_dp.format(cnr_result.noise) +
                    "\nTTF 50 (1/mm)," + this.two_dp.format(TTF.getTTF50()) +
                    "\nTTF 10 (1/mm)," + this.two_dp.format(TTF.getTTF10()) +
                    "\nTTF centre of mass (x)," + this.two_dp.format(ttf_com[0]) +
                    "\nTTF centre of mass (y)," + this.two_dp.format(ttf_com[1]) +
                    "\nLower integration limit (1/mm),0" +
                    "\nUpper integration limit (1/mm)," + this.two_dp.format(nyquist) +
                    "\nImage pixel size (mm)," + this.two_dp.format(cal.pixelWidth)
            );

            if(TTF.getAverageStackStatus()) {
                IJ.log("\nTTF average stack used,Yes" +
                        "\nTTF min stack index," + TTF.getMinStackIndex() +
                        "\nTTF max stack index," + TTF.getMaxStackIndex()
                );
            }
            else {
                IJ.log("\nTTF average stack used,No" +
                        "\nTTF min stack index,N/A" +
                        "\nTTF max stack index,N/A"
                );
            }

            if(NPS.getSubstackStatus()) {
                IJ.log("\nNPS average stack used,Yes" +
                        "\nNPS min stack index," + NPS.getMinStackIndex() +
                        "\nNPS max stack index," + NPS.getMaxStackIndex()
                );
            }
            else {
                IJ.log("\nNPS average stack used,No" +
                        "\nNPS min stack index,N/A" +
                        "\nNPS max stack index,N/A"
                );
            }

            IJ.log("\nEye model settings used" +
                    "\nEye model," + eye_model.getModelDescription() +
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

            // See if the results filename needs to be updated.
            if(this.results_filename.equals("results.csv") && !TTF.getReconFilter().equals("")) {
                this.results_filename = this.imp.getTitle() + " " + TTF.getReconFilter() + " results.csv";
                this.results_filename = sanitiseFilename(this.results_filename);
            }
        }

        // Save the log file
        IJ.selectWindow("Log");
        String image_path = "";
        try {
            image_path = this.imp.getOriginalFileInfo().directory;
        }
        catch (NullPointerException e) {
            DirectoryChooser dirChooser = new DirectoryChooser("Select folder to save results");
            image_path = dirChooser.getDirectory();
        }


        try {
            IJ.saveAs("Text", image_path + this.results_filename);
        }
        catch(Exception e) {
            IJ.error("There was an error when saving the results file:\n" + e);
        }
    }


    private void calcTaskImageIncTTFandNPS(double task_diameter, UnivariateFunction ttf_fn, UnivariateFunction nps_fn, int imageSize, double pixelSize, double bgdStdev, double d) {
        double freq_inc = 1.0 / (imageSize*pixelSize);

        // Check that imageSize is a power of two to prevent the 2D Fourier transform from failing
        if (!isPowerOfTwo(imageSize)) {
            // Set imageSize to the next power of two if it wasn't already a power of two
            imageSize = getNextPowerOfTwo(imageSize);
        }

        ImagePlus discImage = createIdealDiscImage(task_diameter, pixelSize, this.imp.getCalibration(), imageSize);
        double[][] fullDiscImageData = ArrayConversions.convertFloatToDouble(discImage.getProcessor().getFloatArray());
/*
        //------------------------------------
        // Create perfect disc task image
        ImagePlus discImage = TaskTransferFunction_.createDiscTemplate(task_diameter, pixelSize, this.imp.getCalibration());

        // Pad out the discImage to the correct image size; assuming a square image
        float[][] discImageData = discImage.getProcessor().getFloatArray();
        double[][] fullDiscImageData = new double[imageSize][imageSize];
        if (discImageData.length < imageSize) {
            int offset = (imageSize - discImageData.length)/2;
            int discImageData_row = 0;
            for (int row=offset; row<offset+discImageData.length; row++) {
                int discImageData_col = 0;
                for (int col=offset; col<offset+discImageData.length; col++) {
                    fullDiscImageData[row][col] = discImageData[discImageData_row][discImageData_col];
                    discImageData_col++;
                }
                discImageData_row++;
            }
        }
        //------------------------------------*/


        //------------------------------------
        // Convolve perfect disc with the task transfer function
        TwoDFFT twod_fourier_transformer = new TwoDFFT();
        Complex[][] discFreqImage = twod_fourier_transformer.transform(fullDiscImageData, DftNormalization.STANDARD, TransformType.FORWARD);
        discFreqImage = twod_fourier_transformer.swapQuadrants(discFreqImage);
        int zero_u = imageSize/2; // u location of zero frequency element
        int zero_v = imageSize/2; // v location of zero frequency element
        for (int u=0; u<imageSize; u++) {
            for (int v=0; v<imageSize; v++) {
                double freq = freq_inc * calculateDistanceBetweenPoints(zero_u, zero_v, u, v);
                try {
                    discFreqImage[u][v] = discFreqImage[u][v].multiply(ttf_fn.value(freq));
                } catch (Exception e) {
                    discFreqImage[u][v] = Complex.ZERO;
                }
            }
        }
        //------------------------------------


        //------------------------------------
        // FT the frequency image back to the spatial domain and take abs of the values; scale by the required contrast
        discFreqImage = twod_fourier_transformer.swapQuadrants(discFreqImage);
        Complex[][] discImageIncTTF = twod_fourier_transformer.transform(discFreqImage, DftNormalization.STANDARD, TransformType.INVERSE);
        for (int x=0; x<imageSize; x++) {
            for (int y=0; y<imageSize; y++) {
                fullDiscImageData[x][y] = discImageIncTTF[x][y].abs() * this.task_contrast;
            }
        }
        //------------------------------------


        //------------------------------------
        // Create a noise image with zero mean, appropriate standard deviation and power spectrum
        CreateNoiseImage_ noiseImage = new CreateNoiseImage_();
        noiseImage.setImageSize(imageSize);
        noiseImage.setPixelSize(pixelSize);
        noiseImage.setMean(0.0);
        noiseImage.setSd(bgdStdev);
        noiseImage.createNoiseImage(CreateNoiseImage_.WHITE);
        noiseImage.setNoisePowerFn(nps_fn);
        noiseImage.applyNoisePowerFn();
        noiseImage.applyStDevAndMean();
        noiseImage.showNoiseImage();

        System.out.print("Standard deviation for noise image is: " + bgdStdev);
        //------------------------------------


        //------------------------------------
        // Add the noise to the task_image
        float[][] taskNoiseImage = noiseImage.getNoiseImage().getProcessor().getFloatArray();
        for (int x=0; x<imageSize; x++) {
            for (int y=0; y<imageSize; y++) {
                taskNoiseImage[x][y] = taskNoiseImage[x][y] + (float)fullDiscImageData[x][y];
            }
        }
        //------------------------------------


        //------------------------------------
        // Display the task image
        ImagePlus finalTaskImage = NewImage.createFloatImage("Task image (d'=" + this.three_dp.format(d) + ")", imageSize, imageSize, 1, NewImage.FILL_BLACK);
        IJ.run(finalTaskImage, "Properties...", "unit=mm pixel_width=" + pixelSize + " pixel_height=" + pixelSize);
        finalTaskImage.getProcessor().setFloatArray(taskNoiseImage);
        finalTaskImage.show();
        finalTaskImage.updateAndDraw();
        IJ.run(finalTaskImage, "Enhance Contrast", "saturated=0.35");
        //------------------------------------
    }


    public ImagePlus createIdealDiscImage(double disc_diameter, double image_pixel_size_mm, Calibration cal, int image_pixel_width) {
        ImagePlus idealDiscImage = TaskTransferFunction_.createDiscTemplate(disc_diameter, image_pixel_size_mm, cal);

        // Pad out the idealDiscImage to the correct image size; assuming a square image
        float[][] discImageData = idealDiscImage.getProcessor().getFloatArray();
        float[][] fullDiscImageData = new float[image_pixel_width][image_pixel_width];
        if (discImageData.length < image_pixel_width) {
            int offset = (image_pixel_width - discImageData.length)/2;
            int discImageData_row = 0;
            for (int row=offset; row<offset+discImageData.length; row++) {
                int discImageData_col = 0;
                for (int col=offset; col<offset+discImageData.length; col++) {
                    fullDiscImageData[row][col] = discImageData[discImageData_row][discImageData_col];
                    discImageData_col++;
                }
                discImageData_row++;
            }

            ImagePlus paddedIdealDiscImage = NewImage.createFloatImage("Ideal disc image", image_pixel_width, image_pixel_width, 1, NewImage.FILL_BLACK);
            paddedIdealDiscImage.getProcessor().setFloatArray(fullDiscImageData);
            paddedIdealDiscImage.setCalibration(cal);
            return paddedIdealDiscImage;
        }

        return idealDiscImage;
    }

    /**
     * Method to calculate the distance between two cartesian coordinates
     *
     * @param x1 X location of point 1
     * @param y1 Y location of point 2
     * @param x2 X location of point 1
     * @param y2 y location of point 2
     * @return The straight-line distance between point 1 and 2
     */
    private static double calculateDistanceBetweenPoints(double x1, double y1, double x2,double y2) {
        // Calculate the distance between two points in cartesian space.
        return Math.sqrt( (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }


    boolean isPowerOfTwo(int n) {
        if (n==0) return false;

        return (int)(Math.ceil((Math.log(n) / Math.log(2)))) == (int)(Math.floor(((Math.log(n) / Math.log(2)))));
    }


    int getNextPowerOfTwo(int n) {
        int highestOneBit = Integer.highestOneBit(n);
        if (n == highestOneBit) return n;

        return highestOneBit << 1;
    }


    /**
     * Method to replace every character that is not a letter, number, underscore or dot with an underscore using regex.
     * Taken from:
     * https://stackoverflow.com/questions/1184176/how-can-i-safely-encode-a-string-in-java-to-use-as-a-filename
     * @param inputName
     * @return The sanitised string
     */
    public String sanitiseFilename(String inputName) {
        return inputName.replaceAll("[^a-zA-Z0-9-_\\.]", "_");
    }
}