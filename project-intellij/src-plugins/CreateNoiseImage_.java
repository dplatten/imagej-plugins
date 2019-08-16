import ij.ImagePlus;
import ij.gui.NewImage;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.TransformType;

public class CreateNoiseImage_  implements PlugIn {

    private int imageSize = 512;
    private double pixelSize = 0.5;
    private double meanValue = 100.0;
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

    void createWhiteNoiseImage() {
        this.noiseImage = NewImage.createImage("NoiseImage", this.imageSize, this.imageSize, 1, 32, NewImage.FILL_RANDOM);
        ImageProcessor data = this.noiseImage.getProcessor();

        data.multiply(this.sdValue);
        data.add(this.meanValue);
    }

    void showNoiseImage() {
        this.noiseImage.show();
    }

    void applyNoisePowerFn() {
        double[][] noiseValues = convertFloatToDouble(this.noiseImage.getProcessor().getFloatArray());
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
        this.noiseImage.getProcessor().setFloatArray(convertDoubleToFloat(noiseValues));
        //ImagePlus filteredNoiseImage = NewImage.createImage("NoiseImage", this.imageSize, this.imageSize, 1, 32, NewImage.FILL_BLACK);
        //filteredNoiseImage.getProcessor().setFloatArray(convertDoubleToFloat(noiseValues));
        //filteredNoiseImage.show();
    }


    /**
     *
     * @param values Two-dimensional array of float values
     * @return Two-dimensional array of double values
     */
    private static double[][] convertFloatToDouble(final float[][] values) {
        int width = values.length;
        int height = values[0].length;

        final double[][] result = new double[width][height];

        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                result[i][j] = values[i][j];
            }
        }
        return result;
    }


    /**
     *
     * @param values Two-dimensional array of Double values
     * @return Two-dimensional array of Float values
     */
    private static float[][] convertDoubleToFloat(final double[][] values) {
        int width = values.length;
        int height = values[0].length;

        final float[][] result = new float[width][height];

        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                result[i][j] = (float) values[i][j];
            }
        }
        return result;
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
