import ij.IJ;
import ij.ImagePlus;
import ij.gui.NewImage;
import ij.plugin.FFT;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.TransformType;
import java.util.Random;

public class CreateNoiseImage_  implements PlugIn {

    static final int WHITE = 0;
    static final int GAUSSIAN = 1;

    private int imageSize = 512;
    private double pixelSize = 0.5;
    private double meanValue = 0.0;
    private double sdValue = 50.0;
    private ImagePlus noiseImage;
    private UnivariateFunction noisePowerFn;

    void setImageSize(int newSize) { this.imageSize = newSize; }
    void setPixelSize(double newSize) { this.pixelSize = newSize; }
    void setMean(double newMean) { this.meanValue = newMean; }
    void setSd(double newSd) { this.sdValue = newSd; }
    void setNoisePowerFn(UnivariateFunction newNoisePowerFn) { this.noisePowerFn = newNoisePowerFn; }

    ImagePlus getNoiseImage() { return this.noiseImage; }

    int setup() {
        return 0;
    }

    public void run(String arg) {
    }


    void createNoiseImage(int noiseType) {
        Random random = new Random(System.currentTimeMillis());

        float[][] noiseData = new float[this.imageSize][this.imageSize];

        for (int x=0; x<this.imageSize; x++) {
            for (int y=0; y<this.imageSize; y++) {
                if (noiseType == WHITE) {
                    noiseData[x][y] = random.nextFloat(); // Random number between 0.0 and 1.0
                }
                else {
                    // nextGaussian values have mean of 0.0 and standard deviation of 1.0
                    // Multiply by the required standard deviation
                    // https://www.javamex.com/tutorials/random_numbers/gaussian_distribution_2.shtml
                    noiseData[x][y] = (float) (random.nextGaussian());
                }
            }
        }

        this.noiseImage = NewImage.createImage("NoiseImage", this.imageSize, this.imageSize, 1, 32, NewImage.FILL_BLACK);
        this.noiseImage.getProcessor().setFloatArray(noiseData);
    }


    void showNoiseImage() {
        this.noiseImage.show();
    }


    void applyNoisePowerFn() {
        ImagePlus noiseFilter = NewImage.createFloatImage("NoiseFilter", imageSize, imageSize, 1, NewImage.FILL_BLACK);

        float[][] noiseFilterValues = noiseFilter.getProcessor().getFloatArray();

        double freq_inc = 1.0 / (this.imageSize*this.pixelSize);
        int zero_u = this.imageSize/2; // u location of zero frequency element
        int zero_v = this.imageSize/2; // v location of zero frequency element

        for (int u=0; u<this.imageSize; u++) {
            for (int v=0; v<this.imageSize; v++) {
                double freq = freq_inc * calculateDistanceBetweenPoints(zero_u, zero_v, u, v);
                try {
                    noiseFilterValues[u][v] = (float) this.noisePowerFn.value(freq);
                } catch (Exception e) {
                    noiseFilterValues[u][v] = 0;
                }
            }
        }
        noiseFilter.getProcessor().setFloatArray(noiseFilterValues);

        FFT.filter(this.noiseImage, noiseFilter);
    }


    void applyStDevAndMean() {
        ImageProcessor data = this.noiseImage.getProcessor();
        ImageStatistics stats = data.getStatistics();

        data.subtract(stats.mean); // Makes the mean zero
        data.multiply(this.sdValue / stats.stdDev); // Sets the standard deviation to sdValue
        data.add(this.meanValue); // Sets mean to meanValue
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
