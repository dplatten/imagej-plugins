package QuantitativeIQ;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

public class TwoDFFT {

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
        two_d_fft = mirrorOnDiagonal(two_d_fft);

        // Fourier transform each column of the initial result
        for (int col=0; col<two_d_fft.length; col++) {
            two_d_fft[col] = fourier_transformer.transform(two_d_fft[col], TransformType.FORWARD);
        }

        // Mirror the results back to their original orientation
        two_d_fft = mirrorOnDiagonal(two_d_fft);

        return two_d_fft;
    }


    private Complex[][] mirrorOnDiagonal(Complex[][] input) {
        Complex[][] output = new Complex[input[0].length][input.length];

        for (int i=0; i<input.length; i++) {
            for (int j=0; j<input[0].length; j++) {
                output[j][i] = input[i][j];
            }
        }

        return output;
    }
}
