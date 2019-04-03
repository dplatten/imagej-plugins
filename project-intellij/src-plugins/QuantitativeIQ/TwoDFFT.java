package QuantitativeIQ;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

public class TwoDFFT {

    /** A method to carry out a two-dimensional forward Fourier transform.
     *  The frequency scale of the resulting transform can be calculated using:
     *  freq_inc_i [per mm] = 1.0 / (total width of input array [mm]);
     *  freq_inc_j [per mm] = 1.0 / (total height of input array [mm]);
     *
     * @param input_array A two-dimensional array of double-precision values. The length
     *                    of each array dimension must be a power of two
     * @param norm_type The type of normalisation to carry out. Can be DftNormalization.STANDARD
     *                  or DftNormalization.UNITARY.
     * @return A two-dimensional array of complex numbers containing the result of the forward
     *         Fourier transform. The return array has the same dimensions as the input array.
     */
    Complex[][] transform(double[][] input_array, DftNormalization norm_type) {

        Complex[][] two_d_fft = new Complex[input_array.length][input_array[0].length];

        // Need to Fourier transform each column in the input array,
        // and then Fourier transform each row of the result of that.
        FastFourierTransformer fourier_transformer = new FastFourierTransformer(norm_type);

        // Fourier transform each column in the input array
        for (int col=0; col<input_array.length; col++) {
            two_d_fft[col] = fourier_transformer.transform(input_array[col], TransformType.FORWARD);
        }

        // Fourier transform each row of the initial result
        // Mirror the initial results along the diagonal so rows become columns
        two_d_fft = transpose(two_d_fft);

        // Fourier transform each column of the initial result
        for (int col=0; col<two_d_fft.length; col++) {
            two_d_fft[col] = fourier_transformer.transform(two_d_fft[col], TransformType.FORWARD);
        }

        // Mirror the results back to their original orientation
        two_d_fft = transpose(two_d_fft);

        return two_d_fft;
    }


    Complex[][] transform(Complex[][] input_array, DftNormalization norm_type, TransformType trans_type) {
        Complex[][] two_d_fft = new Complex[input_array.length][input_array[0].length];

        // Need to Fourier transform each column in the input array,
        // and then Fourier transform each row of the result of that.
        FastFourierTransformer fourier_transformer = new FastFourierTransformer(norm_type);

        // Fourier transform each column in the input array
        for (int col=0; col<input_array.length; col++) {
            two_d_fft[col] = fourier_transformer.transform(input_array[col], trans_type);
        }

        // Fourier transform each row of the initial result
        // Mirror the initial results along the diagonal so rows become columns
        two_d_fft = transpose(two_d_fft);

        // Fourier transform each column of the initial result
        for (int col=0; col<two_d_fft.length; col++) {
            two_d_fft[col] = fourier_transformer.transform(two_d_fft[col], trans_type);
        }

        // Mirror the results back to their original orientation
        two_d_fft = transpose(two_d_fft);

        return two_d_fft;
    }


    /** Method to transpose a 2d array of complex numbers: the row and column indices are swapped.
     *   1   2   3   4                 1   5   9  13
     *   5   6   7   8    Transpose    2   6  10  14
     *   9  10  11  12  <----------->  3   7  11  15
     *  13  14  15  16                 4   8  12  16
     *
     * @param input A two-dimensional array of complex numbers
     * @return The transpose of the input, also a two-dimensional array of complex numbers
     */
    private Complex[][] transpose(Complex[][] input) {
        Complex[][] output = new Complex[input[0].length][input.length];

        for (int i=0; i<input.length; i++) {
            for (int j=0; j<input[0].length; j++) {
                output[j][i] = input[i][j];
            }
        }

        return output;
    }
}
