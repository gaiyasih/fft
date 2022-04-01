using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace FFT {
    public class tool_1 {
        public static Complex[] fft(float[] p1) {
            //fft要計算幾輪
            int p1Dim = p1.Length;
            int ffttimes = (int)Math.Ceiling(Math.Log(p1Dim, 2));
            int resultDim = (int)Math.Pow(2, ffttimes);
            Complex[] result = new Complex[resultDim];

            //準備輸入vector
            float[] temp_p1 = new float[resultDim];
            for (int i = 0; i < p1Dim; i++) {
                temp_p1[i] = p1[i];
            }
            int[] temp_p2 = new int[resultDim];
            for (int i = 0; i < resultDim; i++) {
                int[] a1 = new int[ffttimes];
                int temp_a1 = i;
                for (int j = 0; j < ffttimes; j++) {
                    a1[j] = temp_a1 % 2;
                    temp_a1 = temp_a1 / 2;
                }
                Array.Reverse(a1);
                for (int k = 0; k < ffttimes; k++) {
                    temp_p2[i] += (int)(a1[k] * Math.Pow(2, k));
                }
            }
            for (int i = 0; i < resultDim; i++) {
                result[i] = temp_p1[temp_p2[i]];
            }

            //準備旋轉因子
            Complex[] twiddle_factor = new Complex[resultDim / 2];
            for (int i = 0; i < resultDim / 2; i++) {
                twiddle_factor[i] = Complex.Exp(new Complex(0, -2 * Math.PI * ((float)i / resultDim)));
            }

            //fft主演算法
            for (int i = 0; i < ffttimes; i++) {
                float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                int temp_a2 = (int)Math.Pow(2, i);
                int temp_a3 = (int)Math.Pow(2, i + 1);
                int temp_a4 = resultDim / temp_a3;
                Complex[] temp_result = new Complex[resultDim];
                for (int k = 0; k < resultDim; k++) {
                    temp_result[k] = result[k];
                }
                for (int j = 0; j < resultDim; j++) {
                    if ((j % temp_a3) < temp_a1) {
                        result[j] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                    } else {
                        result[j] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                    }
                }
            }

            return result;
        }
        public static Complex[] fft(Complex[] p1) {
            //fft要計算幾輪
            int p1Dim = p1.Length;
            int ffttimes = (int)Math.Ceiling(Math.Log(p1Dim, 2));
            int resultDim = (int)Math.Pow(2, ffttimes);
            Complex[] result = new Complex[resultDim];

            //準備輸入vector
            Complex[] temp_p1 = new Complex[resultDim];
            for (int i = 0; i < p1Dim; i++) {
                temp_p1[i] = p1[i];
            }
            int[] temp_p2 = new int[resultDim];
            for (int i = 0; i < resultDim; i++) {
                int[] a1 = new int[ffttimes];
                int temp_a1 = i;
                for (int j = 0; j < ffttimes; j++) {
                    a1[j] = temp_a1 % 2;
                    temp_a1 = temp_a1 / 2;
                }
                Array.Reverse(a1);
                for (int k = 0; k < ffttimes; k++) {
                    temp_p2[i] += (int)(a1[k] * Math.Pow(2, k));
                }
            }
            for (int i = 0; i < resultDim; i++) {
                result[i] = temp_p1[temp_p2[i]];
            }

            //準備旋轉因子
            Complex[] twiddle_factor = new Complex[resultDim / 2];
            for (int i = 0; i < resultDim / 2; i++) {
                twiddle_factor[i] = Complex.Exp(new Complex(0, -2 * Math.PI * ((float)i / resultDim)));
            }

            //fft主演算法
            for (int i = 0; i < ffttimes; i++) {
                float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                int temp_a2 = (int)Math.Pow(2, i);
                int temp_a3 = (int)Math.Pow(2, i + 1);
                int temp_a4 = resultDim / temp_a3;
                Complex[] temp_result = new Complex[resultDim];
                for (int k = 0; k < resultDim; k++) {
                    temp_result[k] = result[k];
                }
                for (int j = 0; j < resultDim; j++) {
                    if ((j % temp_a3) < temp_a1) {
                        result[j] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                    } else {
                        result[j] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                    }
                }
            }

            return result;
        }
        public static float[] ifft(Complex[] p1) {
            //ifft要計算幾輪
            int p1Dim = p1.Length;
            int ffttimes = (int)Math.Ceiling(Math.Log(p1Dim, 2));
            int resultDim = (int)Math.Pow(2, ffttimes);
            float[] result = new float[resultDim];
            Complex[] complex_result = new Complex[resultDim];

            //準備輸入vector
            Complex[] temp_p1 = new Complex[resultDim];
            for (int i = 0; i < p1Dim; i++) {
                temp_p1[i] = p1[i];
            }
            int[] temp_p2 = new int[resultDim];
            for (int i = 0; i < resultDim; i++) {
                int[] a1 = new int[ffttimes];
                int temp_a1 = i;
                for (int j = 0; j < ffttimes; j++) {
                    a1[j] = temp_a1 % 2;
                    temp_a1 = temp_a1 / 2;
                }
                Array.Reverse(a1);
                for (int k = 0; k < ffttimes; k++) {
                    temp_p2[i] += (int)(a1[k] * Math.Pow(2, k));
                }
            }
            for (int i = 0; i < resultDim; i++) {
                complex_result[i] = temp_p1[temp_p2[i]];
            }

            //準備旋轉因子
            Complex[] twiddle_factor = new Complex[resultDim / 2];
            for (int i = 0; i < resultDim / 2; i++) {
                twiddle_factor[i] = Complex.Exp(new Complex(0, 2 * Math.PI * ((float)i / resultDim)));
            }

            //ifft主演算法
            for (int i = 0; i < ffttimes; i++) {
                float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                int temp_a2 = (int)Math.Pow(2, i);
                int temp_a3 = (int)Math.Pow(2, i + 1);
                int temp_a4 = resultDim / temp_a3;
                Complex[] temp_result = new Complex[resultDim];
                for (int k = 0; k < resultDim; k++) {
                    temp_result[k] = complex_result[k];
                }
                for (int j = 0; j < resultDim; j++) {
                    if ((j % temp_a3) < temp_a1) {
                        complex_result[j] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                    } else {
                        complex_result[j] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                    }
                }
            }

            //將complex轉成real
            for (int i = 0; i < resultDim; i++) {
                result[i] = (float)(complex_result[i].Real / resultDim);
            }

            return result;
        }
        public static Complex[] temp_ifft(Complex[] p1) {
            //ifft要計算幾輪
            int p1Dim = p1.Length;
            int ffttimes = (int)Math.Ceiling(Math.Log(p1Dim, 2));
            int resultDim = (int)Math.Pow(2, ffttimes);
            Complex[] result = new Complex[resultDim];
            Complex[] complex_result = new Complex[resultDim];

            //準備輸入vector
            Complex[] temp_p1 = new Complex[resultDim];
            for (int i = 0; i < p1Dim; i++) {
                temp_p1[i] = p1[i];
            }
            int[] temp_p2 = new int[resultDim];
            for (int i = 0; i < resultDim; i++) {
                int[] a1 = new int[ffttimes];
                int temp_a1 = i;
                for (int j = 0; j < ffttimes; j++) {
                    a1[j] = temp_a1 % 2;
                    temp_a1 = temp_a1 / 2;
                }
                Array.Reverse(a1);
                for (int k = 0; k < ffttimes; k++) {
                    temp_p2[i] += (int)(a1[k] * Math.Pow(2, k));
                }
            }
            for (int i = 0; i < resultDim; i++) {
                complex_result[i] = temp_p1[temp_p2[i]];
            }

            //準備旋轉因子
            Complex[] twiddle_factor = new Complex[resultDim / 2];
            for (int i = 0; i < resultDim / 2; i++) {
                twiddle_factor[i] = Complex.Exp(new Complex(0, 2 * Math.PI * ((float)i / resultDim)));
            }

            //ifft主演算法
            for (int i = 0; i < ffttimes; i++) {
                float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                int temp_a2 = (int)Math.Pow(2, i);
                int temp_a3 = (int)Math.Pow(2, i + 1);
                int temp_a4 = resultDim / temp_a3;
                Complex[] temp_result = new Complex[resultDim];
                for (int k = 0; k < resultDim; k++) {
                    temp_result[k] = complex_result[k];
                }
                for (int j = 0; j < resultDim; j++) {
                    if ((j % temp_a3) < temp_a1) {
                        complex_result[j] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                    } else {
                        complex_result[j] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                    }
                }
            }

            for (int i = 0; i < resultDim; i++) {
                result[i] = complex_result[i] / resultDim;
            }

            return result;
        }
        public static Complex[,] fft2(float[,] p1) {
            //準備輸入資料
            int p1Col = p1.GetLength(0);
            int p1Row = p1.GetLength(1);
            int resultCol = (int)Math.Pow(2, Math.Ceiling(Math.Log(p1Col, 2)));
            int resultRow = (int)Math.Pow(2, Math.Ceiling(Math.Log(p1Row, 2)));
            Complex[,] result = new Complex[resultCol, resultRow];
            Complex[] temp_resultCol = new Complex[resultCol];
            Complex[] temp_resultRow = new Complex[resultRow];

            //將Col分別FFT
            for (int i = 0; i < p1Row; i++) {
                float[] temp_p1 = new float[resultCol];
                for (int j = 0; j < p1Col; j++) {
                    temp_p1[j] = p1[j, i];
                }
                temp_resultCol = fft(temp_p1);
                for (int j = 0; j < resultCol; j++) {
                    result[j, i] = temp_resultCol[j];
                }
            }

            //將Row分別FFT
            for (int i = 0; i < resultCol; i++) {
                Complex[] temp_p1 = new Complex[resultRow];
                for (int j = 0; j < p1Row; j++) {
                    temp_p1[j] = result[i, j];
                }
                temp_resultRow = fft(temp_p1);
                for (int j = 0; j < resultRow; j++) {
                    result[i, j] = temp_resultRow[j];
                }
            }

            return result;
        }
        public static float[,] ifft2(Complex[,] p1) {
            //準備輸入資料
            int resultCol = (int)Math.Pow(2, Math.Ceiling(Math.Log(p1.GetLength(0), 2)));
            int resultRow = (int)Math.Pow(2, Math.Ceiling(Math.Log(p1.GetLength(1), 2)));
            float[,] result = new float[resultCol, resultRow];
            Complex[,] temp_result = new Complex[resultCol, resultRow];
            Complex[] temp_resultCol = new Complex[resultCol];
            float[] temp_resultRow = new float[resultRow];

            //將Col分別ifft
            for (int i = 0; i < resultRow; i++) {
                Complex[] temp_p1 = new Complex[resultCol];
                for (int j = 0; j < resultCol; j++) {
                    temp_p1[j] = p1[j, i];
                }
                temp_resultCol = temp_ifft(temp_p1);
                for (int j = 0; j < resultCol; j++) {
                    temp_result[j, i] = temp_resultCol[j];
                }
            }

            //將Row分別ifft
            for (int i = 0; i < resultCol; i++) {
                Complex[] temp_p1 = new Complex[resultRow];
                for (int j = 0; j < resultRow; j++) {
                    temp_p1[j] = temp_result[i, j];
                }
                temp_resultRow = ifft(temp_p1);
                for (int j = 0; j < resultRow; j++) {
                    result[i, j] = temp_resultRow[j];
                }
            }

            return result;
        }





















    }
}
