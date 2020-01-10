

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

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

class TwoDFFT {

    /** A method to carry out a two-dimensional forward Fourier transform.
     *  <br><br>
     *  The frequency scale of the resulting transform can be calculated using:
     *  <br><br>
     *  freq_inc_i [per mm] = 1.0 / (total width of input array [mm]);
     *  <br><br>
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


    /** A method to carry out a two-dimensional forward or reverse Fourier transform.
     *
     * @param input_array A two-dimensional Complex array of values. The length of
     *                    each array dimension must be a power of two
     * @param norm_type The type of normalisation to carry out. Can be DftNormalization.STANDARD
     *                  or DftNormalization.UNITARY.
     * @param trans_type TransformType.INVERSE or TransformType.FORWARD
     * @return A two-dimensional array of complex numbers containing the result of the Fourier transform.
     *         The return array has the same dimensions as the input array.
     */
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


    /** A method to carry out a two-dimensional forward or reverse Fourier transform.
     *
     * @param input_array A two-dimensional double array of values. The length of
     *                    each array dimension must be a power of two
     * @param norm_type The type of normalisation to carry out. Can be DftNormalization.STANDARD
     *                  or DftNormalization.UNITARY.
     * @param trans_type TransformType.INVERSE or TransformType.FORWARD
     * @return A two-dimensional array of complex numbers containing the result of the Fourier transform.
     *         The return array has the same dimensions as the input array.
     */
    Complex[][] transform(double[][] input_array, DftNormalization norm_type, TransformType trans_type) {
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


    /** Method to swap the quadrants of a two dimensional frequency domain array that is the result
     * of a 2d Fourier transform of an image.
     *
     * @param input_array Array of frequency data with zero frequency at 0,0
     * @return Array of frequency data with zero frequency at the centre
     */
    Complex[][] swapQuadrants(Complex[][] input_array) {
        // Shift the quadrants so that zero is at the centre
        Complex[][] shifted_fft = new Complex[input_array.length][input_array[0].length];
        System.arraycopy(input_array, (int) (input_array.length / 2.0), shifted_fft, 0, (int) (input_array.length / 2.0));
        System.arraycopy(input_array, 0, shifted_fft, (int) (input_array.length / 2.0), (int) (input_array.length / 2.0));

        input_array = transpose(shifted_fft);
        shifted_fft = transpose(shifted_fft);
        System.arraycopy(input_array, (int) (input_array.length / 2.0), shifted_fft, 0, (int) (input_array.length / 2.0));
        System.arraycopy(input_array, 0, shifted_fft, (int) (input_array.length / 2.0), (int) (input_array.length / 2.0));

        return transpose(shifted_fft);
    }


    /** Method to swap the quadrants of a two dimensional frequency domain array that is the result
     * of a 2d Fourier transform of an image.
     *
     * @param input_array Array of frequency data with zero frequency at 0,0
     * @return Array of frequency data with zero frequency at the centre
     */
    double[][] swapQuadrants(double[][] input_array) {
        // Shift the quadrants so that zero is at the centre
        double[][] shifted_fft = new double[input_array.length][input_array[0].length];
        System.arraycopy(input_array, (int) (input_array.length / 2.0), shifted_fft, 0, (int) (input_array.length / 2.0));
        System.arraycopy(input_array, 0, shifted_fft, (int) (input_array.length / 2.0), (int) (input_array.length / 2.0));

        input_array = transpose(shifted_fft);
        shifted_fft = transpose(shifted_fft);
        System.arraycopy(input_array, (int) (input_array.length / 2.0), shifted_fft, 0, (int) (input_array.length / 2.0));
        System.arraycopy(input_array, 0, shifted_fft, (int) (input_array.length / 2.0), (int) (input_array.length / 2.0));

        return transpose(shifted_fft);
    }


    /** Method to transpose a 2d array of complex numbers: the row and column indices are swapped.
     *   <pre>
     *    1   2   3   4    Transpose       1   5   9  13
     *    5   6   7   8    moves between   2   6  10  14
     *    9  10  11  12    these two       3   7  11  15
     *   13  14  15  16    arrays          4   8  12  16
     *  </pre>
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


    /** Method to transpose a 2d double array of numbers: the row and column indices are swapped.
     *   <pre>
     *    1   2   3   4    Transpose       1   5   9  13
     *    5   6   7   8    moves between   2   6  10  14
     *    9  10  11  12    these two       3   7  11  15
     *   13  14  15  16    arrays          4   8  12  16
     *  </pre>
     *
     * @param input A two-dimensional double array of numbers
     * @return The transpose of the input, also a two-dimensional array of complex numbers
     */
    private double[][] transpose(double[][] input) {
        double[][] output = new double[input[0].length][input.length];

        for (int i=0; i<input.length; i++) {
            for (int j=0; j<input[0].length; j++) {
                output[j][i] = input[i][j];
            }
        }

        return output;
    }

}
