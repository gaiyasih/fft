using System;
using System.Threading.Tasks;
using System.Numerics;

namespace FFT {
    public class tool_2 {
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
                Parallel.For(0, resultDim, j => {
                    if ((j % temp_a3) < temp_a1) {
                        result[j] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                    } else {
                        result[j] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                    }
                });
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
                Parallel.For(0, resultDim, j => {
                    if ((j % temp_a3) < temp_a1) {
                        complex_result[j] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                    } else {
                        complex_result[j] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                    }
                });
            }

            //將complex轉成real
            for (int i = 0; i < resultDim; i++) {
                result[i] = (float)(complex_result[i].Real / resultDim);
            }

            return result;
        }
        public static Complex[,] fft2(float[,] p1) {
            //準備輸入輸出data
            int p1Col = p1.GetLength(0);
            int p1Row = p1.GetLength(1);
            int ffttimesCol = (int)Math.Ceiling(Math.Log(p1Col, 2));
            int ffttimesRow = (int)Math.Ceiling(Math.Log(p1Row, 2));
            int resultCol = (int)Math.Pow(2, ffttimesCol);
            int resultRow = (int)Math.Pow(2, ffttimesRow);
            Complex[,] result = new Complex[resultCol, resultRow];
            
            //計算Col fft結果
            Parallel.For(0, p1Row, numberRow => {
                float[] temp_p1 = new float[resultCol];
                for (int i = 0; i < p1Col; i++) {
                    temp_p1[i] = p1[i, numberRow];
                }
                int[] temp_p2 = new int[resultCol];
                for (int i = 0; i < resultCol; i++) {
                    int[] a1 = new int[ffttimesCol];
                    int temp_a1 = i;
                    for (int j = 0; j < ffttimesCol; j++) {
                        a1[j] = temp_a1 % 2;
                        temp_a1 = temp_a1 / 2;
                    }
                    Array.Reverse(a1);
                    for (int k = 0; k < ffttimesCol; k++) {
                        temp_p2[i] += (int)(a1[k] * Math.Pow(2, k));
                    }
                }
                for (int i = 0; i < resultCol; i++) {
                    result[i, numberRow] = temp_p1[temp_p2[i]];
                }

                Complex[] twiddle_factor = new Complex[resultCol / 2];
                for (int i = 0; i < resultCol / 2; i++) {
                    twiddle_factor[i] = Complex.Exp(new Complex(0, -2 * Math.PI * ((float)i / resultCol)));
                }

                for (int i = 0; i < ffttimesCol; i++) {
                    float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                    int temp_a2 = (int)Math.Pow(2, i);
                    int temp_a3 = (int)Math.Pow(2, i + 1);
                    int temp_a4 = resultCol / temp_a3;
                    Complex[] temp_result = new Complex[resultCol];
                    for (int k = 0; k < resultCol; k++) {
                        temp_result[k] = result[k, numberRow];
                    }
                    Parallel.For(0, resultCol, j => {
                        if ((j % temp_a3) < temp_a1) {
                            result[j, numberRow] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                        } else {
                            result[j, numberRow] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                        }
                    });
                }
            });
            
            //計算Row fft結果
            Parallel.For(0, resultCol, numberCol => {
                Complex[] temp_p1 = new Complex[resultRow];
                for (int i = 0; i < p1Row; i++) {
                    temp_p1[i] = result[numberCol, i];
                }
                int[] temp_p2 = new int[resultRow];
                for (int i = 0; i < resultRow; i++) {
                    int[] a1 = new int[ffttimesRow];
                    int temp_a1 = i;
                    for (int j = 0; j < ffttimesRow; j++) {
                        a1[j] = temp_a1 % 2;
                        temp_a1 = temp_a1 / 2;
                    }
                    Array.Reverse(a1);
                    for (int k = 0; k < ffttimesRow; k++) {
                        temp_p2[i] += (int)(a1[k] * Math.Pow(2, k));
                    }
                }
                for (int i = 0; i < resultRow; i++) {
                    result[numberCol, i] = temp_p1[temp_p2[i]];
                }

                Complex[] twiddle_factor = new Complex[resultRow / 2];
                for (int i = 0; i < resultRow / 2; i++) {
                    twiddle_factor[i] = Complex.Exp(new Complex(0, -2 * Math.PI * ((float)i / resultRow)));
                }

                for (int i = 0; i < ffttimesRow; i++) {
                    float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                    int temp_a2 = (int)Math.Pow(2, i);
                    int temp_a3 = (int)Math.Pow(2, i + 1);
                    int temp_a4 = resultRow / temp_a3;
                    Complex[] temp_result = new Complex[resultRow];
                    for (int k = 0; k < resultRow; k++) {
                        temp_result[k] = result[numberCol, k];
                    }
                    Parallel.For(0, resultRow, j => {
                        if ((j % temp_a3) < temp_a1) {
                            result[numberCol, j] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                        } else {
                            result[numberCol, j] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                        }
                    });
                }
            });
            
            return result;
        }
        public static float[,] ifft2(Complex[,] p1) {
            //準備輸入輸出data
            int p1Col = p1.GetLength(0);
            int p1Row = p1.GetLength(1);
            int ffttimesCol = (int)Math.Ceiling(Math.Log(p1Col, 2));
            int ffttimesRow = (int)Math.Ceiling(Math.Log(p1Row, 2));
            int resultCol = (int)Math.Pow(2, ffttimesCol);
            int resultRow = (int)Math.Pow(2, ffttimesRow);
            Complex[,] complex_result = new Complex[resultCol, resultRow];
            float[,] result = new float[resultCol, resultRow];

            //計算Col fft結果
            Parallel.For(0, p1Row, numberRow => {
                Complex[] temp_p1 = new Complex[resultCol];
                for (int i = 0; i < p1Col; i++) {
                    temp_p1[i] = p1[i, numberRow];
                }
                int[] temp_p2 = new int[resultCol];
                for (int i = 0; i < resultCol; i++) {
                    int[] a1 = new int[ffttimesCol];
                    int temp_a1 = i;
                    for (int j = 0; j < ffttimesCol; j++) {
                        a1[j] = temp_a1 % 2;
                        temp_a1 = temp_a1 / 2;
                    }
                    Array.Reverse(a1);
                    for (int k = 0; k < ffttimesCol; k++) {
                        temp_p2[i] += (int)(a1[k] * Math.Pow(2, k));
                    }
                }
                for (int i = 0; i < resultCol; i++) {
                    complex_result[i, numberRow] = temp_p1[temp_p2[i]];
                }

                Complex[] twiddle_factor = new Complex[resultCol / 2];
                for (int i = 0; i < resultCol / 2; i++) {
                    twiddle_factor[i] = Complex.Exp(new Complex(0, 2 * Math.PI * ((float)i / resultCol)));
                }

                for (int i = 0; i < ffttimesCol; i++) {
                    float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                    int temp_a2 = (int)Math.Pow(2, i);
                    int temp_a3 = (int)Math.Pow(2, i + 1);
                    int temp_a4 = resultCol / temp_a3;
                    Complex[] temp_result = new Complex[resultCol];
                    for (int k = 0; k < resultCol; k++) {
                        temp_result[k] = complex_result[k, numberRow];
                    }
                    Parallel.For(0, resultCol, j => {
                        if ((j % temp_a3) < temp_a1) {
                            complex_result[j, numberRow] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                        } else {
                            complex_result[j, numberRow] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                        }
                    });
                }
            });

            //計算Row fft結果
            Parallel.For(0, resultCol, numberCol => {
                Complex[] temp_p1 = new Complex[resultRow];
                for (int i = 0; i < p1Row; i++) {
                    temp_p1[i] = complex_result[numberCol, i];
                }
                int[] temp_p2 = new int[resultRow];
                for (int i = 0; i < resultRow; i++) {
                    int[] a1 = new int[ffttimesRow];
                    int temp_a1 = i;
                    for (int j = 0; j < ffttimesRow; j++) {
                        a1[j] = temp_a1 % 2;
                        temp_a1 = temp_a1 / 2;
                    }
                    Array.Reverse(a1);
                    for (int k = 0; k < ffttimesRow; k++) {
                        temp_p2[i] += (int)(a1[k] * Math.Pow(2, k));
                    }
                }
                for (int i = 0; i < resultRow; i++) {
                    complex_result[numberCol, i] = temp_p1[temp_p2[i]];
                }

                Complex[] twiddle_factor = new Complex[resultRow / 2];
                for (int i = 0; i < resultRow / 2; i++) {
                    twiddle_factor[i] = Complex.Exp(new Complex(0, 2 * Math.PI * ((float)i / resultRow)));
                }

                for (int i = 0; i < ffttimesRow; i++) {
                    float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                    int temp_a2 = (int)Math.Pow(2, i);
                    int temp_a3 = (int)Math.Pow(2, i + 1);
                    int temp_a4 = resultRow / temp_a3;
                    Complex[] temp_result = new Complex[resultRow];
                    for (int k = 0; k < resultRow; k++) {
                        temp_result[k] = complex_result[numberCol, k];
                    }
                    Parallel.For(0, resultRow, j => {
                        if ((j % temp_a3) < temp_a1) {
                            complex_result[numberCol, j] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                        } else {
                            complex_result[numberCol, j] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                        }
                    });
                }

                for (int i = 0; i < resultRow; i++) {
                    result[numberCol, i] = (float)(complex_result[numberCol, i].Real / (resultCol * resultRow));
                }
            });

            return result;
        }
    }
}
