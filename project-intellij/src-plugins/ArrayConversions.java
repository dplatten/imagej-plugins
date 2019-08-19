import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;

public class ArrayConversions {

    /**
     *
     * @param values Two-dimensional array of float values
     * @return Two-dimensional array of double values
     */
    public static double[][] convertFloatToDouble(final float[][] values) {
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
    public static float[][] convertDoubleToFloat(final double[][] values) {
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
     *
     * @param values Two-dimensional array of Double values
     * @return Two-dimensional array of Complex values
     */
    public static float[][] convertComplexToFloat(final Complex[][] values) {
        int width = values.length;
        int height = values[0].length;

        final float[][] result = new float[width][height];

        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                result[i][j] = (float) values[i][j].abs();
            }
        }
        return result;
    }


    /**
     * Method to convert a 2d double array to a 2d complex array
     *
     * @param values Two-dimensional array of Double values
     * @return Two-dimensional array of Complex values
     */
    private static Complex[][] convertDoubleToComplex(final double[][] values) {
        int width = values.length;
        int height = values[0].length;

        final Complex[][] result = new Complex[width][height];

        for (int i=0; i<width; i++) {
            result[i] = ComplexUtils.convertToComplex(values[i]);
        }
        return result;
    }
}
