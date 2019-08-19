import ij.IJ;
import ij.ImagePlus;
import ij.gui.NewImage;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.TransformType;

public class CreateNoiseImage_  implements PlugIn {

    private int imageSize = 512;
    private double pixelSize = 0.5;
    private double meanValue = 0.0;
    private double sdValue = 50.0;
    private ImagePlus noiseImage;
    private UnivariateFunction noisePowerFn;

    void setImageSize(int newSize) { this.imageSize = newSize; }
    void setPixelSize(double newSize) { this.pixelSize = newSize; }
    //void setMean(double newMean) { this.meanValue = newMean; }
    void setSd(double newSd) { this.sdValue = newSd; }
    void setNoisePowerFn(UnivariateFunction newNoisePowerFn) { this.noisePowerFn = newNoisePowerFn; }

    ImagePlus getNoiseImage() { return this.noiseImage; }

    int setup() {
        return 0;
    }

    public void run(String arg) {
    }


    void createWhiteNoiseImage() {
        this.noiseImage = NewImage.createImage("NoiseImage", this.imageSize, this.imageSize, 1, 32, NewImage.FILL_RANDOM);
        ImageProcessor data = this.noiseImage.getProcessor();

        data.multiply(this.sdValue);
    }


    void showNoiseImage() {
        this.noiseImage.show();
    }


    void applyNoisePowerFn() {
        double[][] noiseValues = ArrayConversions.convertFloatToDouble(this.noiseImage.getProcessor().getFloatArray());
        TwoDFFT twod_fourier_transformer = new TwoDFFT();
        Complex[][] freqImage = twod_fourier_transformer.transform(noiseValues, DftNormalization.STANDARD, TransformType.FORWARD);
        freqImage = twod_fourier_transformer.swapQuadrants(freqImage);

        double freq_inc = 1.0 / (this.imageSize*this.pixelSize);
        int zero_u = this.imageSize/2; // u location of zero frequency element
        int zero_v = this.imageSize/2; // v location of zero frequency element
        double maxNpsValue = 0.0;
        for (int u=0; u<this.imageSize; u++) {
            for (int v=0; v<this.imageSize; v++) {
                double freq = freq_inc * calculateDistanceBetweenPoints(zero_u, zero_v, u, v);
                try {
                    double currentNoisePower = this.noisePowerFn.value(freq);
                    freqImage[u][v] = freqImage[u][v].multiply(currentNoisePower);
                    if (currentNoisePower > maxNpsValue) maxNpsValue = currentNoisePower;
                } catch (Exception e) {
                    freqImage[u][v] = freqImage[u][v].multiply(0.0);
                }
            }
        }
        for (int u=0; u<this.imageSize; u++) {
            for (int v=0; v<this.imageSize; v++) {
                freqImage[u][v] = freqImage[u][v].divide(maxNpsValue);
            }
        }

        freqImage = twod_fourier_transformer.swapQuadrants(freqImage);
        Complex[][] filteredNoiseValues = twod_fourier_transformer.transform(freqImage, DftNormalization.STANDARD, TransformType.INVERSE);

        for (int x=0; x<this.imageSize; x++) {
            for (int y=0; y<this.imageSize; y++) {
                noiseValues[x][y] = filteredNoiseValues[x][y].abs();
            }
        }

        // Replace the noiseImage data with the filtered data.
        this.noiseImage.getProcessor().setFloatArray(ArrayConversions.convertDoubleToFloat(noiseValues));

        // Ensure resulting noise has mean of zero and required SD
        ImageStatistics noiseImageStats = this.noiseImage.getStatistics();
        IJ.run(this.noiseImage, "Subtract...", "value="+noiseImageStats.mean);
        IJ.run(this.noiseImage, "Divide...", "value="+noiseImageStats.stdDev);
        IJ.run(this.noiseImage, "Multiply...", "value="+this.sdValue);
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
    public static double calculateDistanceBetweenPoints(double x1, double y1, double x2,double y2) {
        // Calculate the distance between two points in cartesian space.
        return Math.sqrt( (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }
}
