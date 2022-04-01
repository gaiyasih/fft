using System;
using System.Threading.Tasks;
using System.Numerics;

namespace FFT {
    public class tool {
        public static Complex[] fft(float[] p1, int p2 = 0) {
            //fft要計算幾輪
            int p1Dim = p1.Length;
            int ffttimes = (int)Math.Ceiling(Math.Log(p1Dim, 2));
            if (p2 != 0) {
                ffttimes = (int)Math.Ceiling(Math.Log(p2, 2));
            }
            int resultDim = (int)Math.Pow(2, ffttimes);
            Complex[] result = new Complex[resultDim];
            //準備輸入vector
            float[] temp_p1 = new float[resultDim];
            for (int i = 0; i < p1Dim; i++) {
                temp_p1[i] = p1[i];
            }
            int[] change_address = new int[resultDim];
            for (int i = 0; i < resultDim; i++) {
                int[] temp_address = new int[ffttimes];
                int temp_i = i;
                for (int j = 0; j < ffttimes; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int j = 0; j < ffttimes; j++) {
                    change_address[i] += (int)(temp_address[j] * Math.Pow(2, j));
                }
            }
            for (int i = 0; i < resultDim; i++) {
                result[i] = temp_p1[change_address[i]];
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
        public static Complex[,] fft2(float[,] p1, int p2 = 0, int p3 = 0) {
            //準備輸入輸出data
            int p1Col = p1.GetLength(0);
            int p1Row = p1.GetLength(1);
            int ffttimesCol = (int)Math.Ceiling(Math.Log(p1Col, 2));
            int ffttimesRow = (int)Math.Ceiling(Math.Log(p1Row, 2));
            if (p2 != 0) {
                ffttimesCol = (int)Math.Ceiling(Math.Log(p2, 2));
            }
            if (p3 != 0) {
                ffttimesRow = (int)Math.Ceiling(Math.Log(p3, 2));
            }
            int resultCol = (int)Math.Pow(2, ffttimesCol);
            int resultRow = (int)Math.Pow(2, ffttimesRow);
            Complex[,] result = new Complex[resultCol, resultRow];
            int[] change_address;
            Complex[] twiddle_factor;
            //計算Col fft結果
            change_address = new int[resultCol];
            for (int i = 0; i < resultCol; i++) {
                int[] temp_address = new int[ffttimesCol];
                int temp_i = i;
                for (int j = 0; j < ffttimesCol; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int j = 0; j < ffttimesCol; j++) {
                    change_address[i] += (int)(temp_address[j] * Math.Pow(2, j));
                }
            }
            twiddle_factor = new Complex[resultCol / 2];
            for (int i = 0; i < resultCol / 2; i++) {
                twiddle_factor[i] = Complex.Exp(new Complex(0, -2 * Math.PI * ((float)i / resultCol)));
            }
            Parallel.For(0, p1Row, numberRow => {
                float[] temp_p1 = new float[resultCol];
                for (int i = 0; i < p1Col; i++) {
                    temp_p1[i] = p1[i, numberRow];
                }
                for (int i = 0; i < resultCol; i++) {
                    result[i, numberRow] = temp_p1[change_address[i]];
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
                    for (int j = 0; j < resultCol; j++) {
                        if ((j % temp_a3) < temp_a1) {
                            result[j, numberRow] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                        } else {
                            result[j, numberRow] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                        }
                    }
                }
            });
            //計算Row fft結果
            change_address = new int[resultRow];
            for (int i = 0; i < resultRow; i++) {
                int[] temp_address = new int[ffttimesRow];
                int temp_i = i;
                for (int j = 0; j < ffttimesRow; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int k = 0; k < ffttimesRow; k++) {
                    change_address[i] += (int)(temp_address[k] * Math.Pow(2, k));
                }
            }
            twiddle_factor = new Complex[resultRow / 2];
            for (int i = 0; i < resultRow / 2; i++) {
                twiddle_factor[i] = Complex.Exp(new Complex(0, -2 * Math.PI * ((float)i / resultRow)));
            }
            Parallel.For(0, resultCol, numberCol => {
                Complex[] temp_p1 = new Complex[resultRow];
                for (int i = 0; i < p1Row; i++) {
                    temp_p1[i] = result[numberCol, i];
                }
                for (int i = 0; i < resultRow; i++) {
                    result[numberCol, i] = temp_p1[change_address[i]];
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
                    for (int j = 0; j < resultRow; j++) {
                        if ((j % temp_a3) < temp_a1) {
                            result[numberCol, j] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                        } else {
                            result[numberCol, j] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                        }
                    }
                }
            });
            return result;
        }
        public static Complex[,,] fft2(float[,,] p1, int p2 = 0, int p3 = 0) {
            //準備輸入輸出data
            int p1Col = p1.GetLength(1);
            int p1Row = p1.GetLength(2);
            int ffttimesCol = (int)Math.Ceiling(Math.Log(p1Col, 2));
            int ffttimesRow = (int)Math.Ceiling(Math.Log(p1Row, 2));
            if (p2 != 0) {
                ffttimesCol = (int)Math.Ceiling(Math.Log(p2, 2));
            }
            if (p3 != 0) {
                ffttimesRow = (int)Math.Ceiling(Math.Log(p3, 2));
            }
            int resultHigh = p1.GetLength(0);
            int resultCol = (int)Math.Pow(2, ffttimesCol);
            int resultRow = (int)Math.Pow(2, ffttimesRow);
            Complex[,,] result = new Complex[resultHigh, resultCol, resultRow];
            int[] change_addressCol = new int[resultCol];
            for (int i = 0; i < resultCol; i++) {
                int[] temp_address = new int[ffttimesCol];
                int temp_i = i;
                for (int j = 0; j < ffttimesCol; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int j = 0; j < ffttimesCol; j++) {
                    change_addressCol[i] += (int)(temp_address[j] * Math.Pow(2, j));
                }
            }
            Complex[] twiddle_factorCol = new Complex[resultCol / 2];
            for (int i = 0; i < resultCol / 2; i++) {
                twiddle_factorCol[i] = Complex.Exp(new Complex(0, -2 * Math.PI * ((float)i / resultCol)));
            }
            int[] change_addressRow = new int[resultRow];
            for (int i = 0; i < resultRow; i++) {
                int[] temp_address = new int[ffttimesRow];
                int temp_i = i;
                for (int j = 0; j < ffttimesRow; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int k = 0; k < ffttimesRow; k++) {
                    change_addressRow[i] += (int)(temp_address[k] * Math.Pow(2, k));
                }
            }
            Complex[] twiddle_factorRow = new Complex[resultRow / 2];
            for (int i = 0; i < resultRow / 2; i++) {
                twiddle_factorRow[i] = Complex.Exp(new Complex(0, -2 * Math.PI * ((float)i / resultRow)));
            }
            Parallel.For(0, resultHigh, numberHigh => {
                //計算Col fft結果
                for (int numberRow = 0; numberRow < p1Row; numberRow++) {
                    float[] temp_p1 = new float[resultCol];
                    for (int i = 0; i < p1Col; i++) {
                        temp_p1[i] = p1[numberHigh, i, numberRow];
                    }
                    for (int i = 0; i < resultCol; i++) {
                        result[numberHigh, i, numberRow] = temp_p1[change_addressCol[i]];
                    }
                    for (int i = 0; i < ffttimesCol; i++) {
                        float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                        int temp_a2 = (int)Math.Pow(2, i);
                        int temp_a3 = (int)Math.Pow(2, i + 1);
                        int temp_a4 = resultCol / temp_a3;
                        Complex[] temp_result = new Complex[resultCol];
                        for (int k = 0; k < resultCol; k++) {
                            temp_result[k] = result[numberHigh, k, numberRow];
                        }
                        for (int j = 0; j < resultCol; j++) {
                            if ((j % temp_a3) < temp_a1) {
                                result[numberHigh, j, numberRow] = temp_result[j] + twiddle_factorCol[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                            } else {
                                result[numberHigh, j, numberRow] = -temp_result[j] * twiddle_factorCol[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                            }
                        }
                    }
                }
                //計算Row fft結果
                for (int numberCol = 0; numberCol < resultCol; numberCol++) {
                    Complex[] temp_p1 = new Complex[resultRow];
                    for (int i = 0; i < p1Row; i++) {
                        temp_p1[i] = result[numberHigh, numberCol, i];
                    }
                    for (int i = 0; i < resultRow; i++) {
                        result[numberHigh, numberCol, i] = temp_p1[change_addressRow[i]];
                    }
                    for (int i = 0; i < ffttimesRow; i++) {
                        float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                        int temp_a2 = (int)Math.Pow(2, i);
                        int temp_a3 = (int)Math.Pow(2, i + 1);
                        int temp_a4 = resultRow / temp_a3;
                        Complex[] temp_result = new Complex[resultRow];
                        for (int k = 0; k < resultRow; k++) {
                            temp_result[k] = result[numberHigh, numberCol, k];
                        }
                        for (int j = 0; j < resultRow; j++) {
                            if ((j % temp_a3) < temp_a1) {
                                result[numberHigh, numberCol, j] = temp_result[j] + twiddle_factorRow[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                            } else {
                                result[numberHigh, numberCol, j] = -temp_result[j] * twiddle_factorRow[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                            }
                        }
                    }
                }
            });
            return result;
        }
        public static Complex[,,,] fft2(float[,,,] p1, int p2 = 0, int p3 = 0) {
            //準備輸入輸出data
            int p1Col = p1.GetLength(2);
            int p1Row = p1.GetLength(3);
            int ffttimesCol = (int)Math.Ceiling(Math.Log(p1Col, 2));
            int ffttimesRow = (int)Math.Ceiling(Math.Log(p1Row, 2));
            if (p2 != 0) {
                ffttimesCol = (int)Math.Ceiling(Math.Log(p2, 2));
            }
            if (p3 != 0) {
                ffttimesRow = (int)Math.Ceiling(Math.Log(p3, 2));
            }
            int resultKer = p1.GetLength(0);
            int resultHigh = p1.GetLength(1);
            int resultCol = (int)Math.Pow(2, ffttimesCol);
            int resultRow = (int)Math.Pow(2, ffttimesRow);
            Complex[,,,] result = new Complex[resultKer, resultHigh, resultCol, resultRow];
            int[] change_addressCol = new int[resultCol];
            for (int i = 0; i < resultCol; i++) {
                int[] temp_address = new int[ffttimesCol];
                int temp_i = i;
                for (int j = 0; j < ffttimesCol; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int j = 0; j < ffttimesCol; j++) {
                    change_addressCol[i] += (int)(temp_address[j] * Math.Pow(2, j));
                }
            }
            Complex[] twiddle_factorCol = new Complex[resultCol / 2];
            for (int i = 0; i < resultCol / 2; i++) {
                twiddle_factorCol[i] = Complex.Exp(new Complex(0, -2 * Math.PI * ((float)i / resultCol)));
            }
            int[] change_addressRow = new int[resultRow];
            for (int i = 0; i < resultRow; i++) {
                int[] temp_address = new int[ffttimesRow];
                int temp_i = i;
                for (int j = 0; j < ffttimesRow; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int k = 0; k < ffttimesRow; k++) {
                    change_addressRow[i] += (int)(temp_address[k] * Math.Pow(2, k));
                }
            }
            Complex[] twiddle_factorRow = new Complex[resultRow / 2];
            for (int i = 0; i < resultRow / 2; i++) {
                twiddle_factorRow[i] = Complex.Exp(new Complex(0, -2 * Math.PI * ((float)i / resultRow)));
            }
            Parallel.For(0, resultKer, numberKer => {
                for (int numberHigh = 0; numberHigh < resultHigh; numberHigh++) {
                    //計算Col fft結果
                    for (int numberRow = 0; numberRow < p1Row; numberRow++) {
                        float[] temp_p1 = new float[resultCol];
                        for (int i = 0; i < p1Col; i++) {
                            temp_p1[i] = p1[numberKer, numberHigh, i, numberRow];
                        }
                        for (int i = 0; i < resultCol; i++) {
                            result[numberKer, numberHigh, i, numberRow] = temp_p1[change_addressCol[i]];
                        }
                        for (int i = 0; i < ffttimesCol; i++) {
                            float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                            int temp_a2 = (int)Math.Pow(2, i);
                            int temp_a3 = (int)Math.Pow(2, i + 1);
                            int temp_a4 = resultCol / temp_a3;
                            Complex[] temp_result = new Complex[resultCol];
                            for (int k = 0; k < resultCol; k++) {
                                temp_result[k] = result[numberKer, numberHigh, k, numberRow];
                            }
                            for (int j = 0; j < resultCol; j++) {
                                if ((j % temp_a3) < temp_a1) {
                                    result[numberKer, numberHigh, j, numberRow] = temp_result[j] + twiddle_factorCol[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                                } else {
                                    result[numberKer, numberHigh, j, numberRow] = -temp_result[j] * twiddle_factorCol[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                                }
                            }
                        }
                    }
                    //計算Row fft結果
                    for (int numberCol = 0; numberCol < resultCol; numberCol++) {
                        Complex[] temp_p1 = new Complex[resultRow];
                        for (int i = 0; i < p1Row; i++) {
                            temp_p1[i] = result[numberKer, numberHigh, numberCol, i];
                        }
                        for (int i = 0; i < resultRow; i++) {
                            result[numberKer, numberHigh, numberCol, i] = temp_p1[change_addressRow[i]];
                        }
                        for (int i = 0; i < ffttimesRow; i++) {
                            float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                            int temp_a2 = (int)Math.Pow(2, i);
                            int temp_a3 = (int)Math.Pow(2, i + 1);
                            int temp_a4 = resultRow / temp_a3;
                            Complex[] temp_result = new Complex[resultRow];
                            for (int k = 0; k < resultRow; k++) {
                                temp_result[k] = result[numberKer, numberHigh, numberCol, k];
                            }
                            for (int j = 0; j < resultRow; j++) {
                                if ((j % temp_a3) < temp_a1) {
                                    result[numberKer, numberHigh, numberCol, j] = temp_result[j] + twiddle_factorRow[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                                } else {
                                    result[numberKer, numberHigh, numberCol, j] = -temp_result[j] * twiddle_factorRow[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                                }
                            }
                        }
                    }
                }
            });
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
            int[] change_address = new int[resultDim];
            for (int i = 0; i < resultDim; i++) {
                int[] temp_address = new int[ffttimes];
                int temp_i = i;
                for (int j = 0; j < ffttimes; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int k = 0; k < ffttimes; k++) {
                    change_address[i] += (int)(temp_address[k] * Math.Pow(2, k));
                }
            }
            for (int i = 0; i < resultDim; i++) {
                complex_result[i] = temp_p1[change_address[i]];
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
            int[] change_address;
            Complex[] twiddle_factor;
            //計算Col fft結果
            change_address = new int[resultCol];
            for (int i = 0; i < resultCol; i++) {
                int[] temp_address = new int[ffttimesCol];
                int temp_i = i;
                for (int j = 0; j < ffttimesCol; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int k = 0; k < ffttimesCol; k++) {
                    change_address[i] += (int)(temp_address[k] * Math.Pow(2, k));
                }
            }
            twiddle_factor = new Complex[resultCol / 2];
            for (int i = 0; i < resultCol / 2; i++) {
                twiddle_factor[i] = Complex.Exp(new Complex(0, 2 * Math.PI * ((float)i / resultCol)));
            }
            Parallel.For(0, p1Row, numberRow => {
                Complex[] temp_p1 = new Complex[resultCol];
                for (int i = 0; i < p1Col; i++) {
                    temp_p1[i] = p1[i, numberRow];
                }
                for (int i = 0; i < resultCol; i++) {
                    complex_result[i, numberRow] = temp_p1[change_address[i]];
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
                    for (int j = 0; j < resultCol; j++) {
                        if ((j % temp_a3) < temp_a1) {
                            complex_result[j, numberRow] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                        } else {
                            complex_result[j, numberRow] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                        }
                    }
                }
            });
            //計算Row fft結果
            change_address = new int[resultRow];
            for (int i = 0; i < resultRow; i++) {
                int[] temp_address = new int[ffttimesRow];
                int temp_i = i;
                for (int j = 0; j < ffttimesRow; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int k = 0; k < ffttimesRow; k++) {
                    change_address[i] += (int)(temp_address[k] * Math.Pow(2, k));
                }
            }
            twiddle_factor = new Complex[resultRow / 2];
            for (int i = 0; i < resultRow / 2; i++) {
                twiddle_factor[i] = Complex.Exp(new Complex(0, 2 * Math.PI * ((float)i / resultRow)));
            }
            Parallel.For(0, resultCol, numberCol => {
                Complex[] temp_p1 = new Complex[resultRow];
                for (int i = 0; i < p1Row; i++) {
                    temp_p1[i] = complex_result[numberCol, i];
                }
                for (int i = 0; i < resultRow; i++) {
                    complex_result[numberCol, i] = temp_p1[change_address[i]];
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
                    for (int j = 0; j < resultRow; j++) {
                        if ((j % temp_a3) < temp_a1) {
                            complex_result[numberCol, j] = temp_result[j] + twiddle_factor[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                        } else {
                            complex_result[numberCol, j] = -temp_result[j] * twiddle_factor[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                        }
                    }
                }
                for (int i = 0; i < resultRow; i++) {
                    result[numberCol, i] = (float)(complex_result[numberCol, i].Real / (resultCol * resultRow));
                }
            });
            return result;
        }
        public static float[,,] ifft2(Complex[,,] p1) {
            //準備輸入輸出data
            int p1Col = p1.GetLength(1);
            int p1Row = p1.GetLength(2);
            int ffttimesCol = (int)Math.Ceiling(Math.Log(p1Col, 2));
            int ffttimesRow = (int)Math.Ceiling(Math.Log(p1Row, 2));
            int resultHigh = p1.GetLength(0);
            int resultCol = (int)Math.Pow(2, ffttimesCol);
            int resultRow = (int)Math.Pow(2, ffttimesRow);
            Complex[,,] complex_result = new Complex[resultHigh, resultCol, resultRow];
            float[,,] result = new float[resultHigh, resultCol, resultRow];
            int[] change_addressCol = new int[resultCol];
            for (int i = 0; i < resultCol; i++) {
                int[] temp_address = new int[ffttimesCol];
                int temp_i = i;
                for (int j = 0; j < ffttimesCol; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int k = 0; k < ffttimesCol; k++) {
                    change_addressCol[i] += (int)(temp_address[k] * Math.Pow(2, k));
                }
            }
            Complex[] twiddle_factorCol = new Complex[resultCol / 2];
            for (int i = 0; i < resultCol / 2; i++) {
                twiddle_factorCol[i] = Complex.Exp(new Complex(0, 2 * Math.PI * ((float)i / resultCol)));
            }
            int[] change_addressRow = new int[resultRow];
            for (int i = 0; i < resultRow; i++) {
                int[] temp_address = new int[ffttimesRow];
                int temp_i = i;
                for (int j = 0; j < ffttimesRow; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int k = 0; k < ffttimesRow; k++) {
                    change_addressRow[i] += (int)(temp_address[k] * Math.Pow(2, k));
                }
            }
            Complex[] twiddle_factorRow = new Complex[resultRow / 2];
            for (int i = 0; i < resultRow / 2; i++) {
                twiddle_factorRow[i] = Complex.Exp(new Complex(0, 2 * Math.PI * ((float)i / resultRow)));
            }
            Parallel.For(0, resultHigh, numberHigh => {
                //計算Col fft結果
                for (int numberRow = 0; numberRow < p1Row; numberRow++) {
                    Complex[] temp_p1 = new Complex[resultCol];
                    for (int i = 0; i < p1Col; i++) {
                        temp_p1[i] = p1[numberHigh, i, numberRow];
                    }
                    for (int i = 0; i < resultCol; i++) {
                        complex_result[numberHigh, i, numberRow] = temp_p1[change_addressCol[i]];
                    }
                    for (int i = 0; i < ffttimesCol; i++) {
                        float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                        int temp_a2 = (int)Math.Pow(2, i);
                        int temp_a3 = (int)Math.Pow(2, i + 1);
                        int temp_a4 = resultCol / temp_a3;
                        Complex[] temp_result = new Complex[resultCol];
                        for (int k = 0; k < resultCol; k++) {
                            temp_result[k] = complex_result[numberHigh, k, numberRow];
                        }
                        for (int j = 0; j < resultCol; j++) {
                            if ((j % temp_a3) < temp_a1) {
                                complex_result[numberHigh, j, numberRow] = temp_result[j] + twiddle_factorCol[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                            } else {
                                complex_result[numberHigh, j, numberRow] = -temp_result[j] * twiddle_factorCol[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                            }
                        }
                    }
                }
                //計算Row fft結果
                for (int numberCol = 0; numberCol < resultCol; numberCol++) {
                    Complex[] temp_p1 = new Complex[resultRow];
                    for (int i = 0; i < p1Row; i++) {
                        temp_p1[i] = complex_result[numberHigh, numberCol, i];
                    }
                    for (int i = 0; i < resultRow; i++) {
                        complex_result[numberHigh, numberCol, i] = temp_p1[change_addressRow[i]];
                    }
                    for (int i = 0; i < ffttimesRow; i++) {
                        float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                        int temp_a2 = (int)Math.Pow(2, i);
                        int temp_a3 = (int)Math.Pow(2, i + 1);
                        int temp_a4 = resultRow / temp_a3;
                        Complex[] temp_result = new Complex[resultRow];
                        for (int k = 0; k < resultRow; k++) {
                            temp_result[k] = complex_result[numberHigh, numberCol, k];
                        }
                        for (int j = 0; j < resultRow; j++) {
                            if ((j % temp_a3) < temp_a1) {
                                complex_result[numberHigh, numberCol, j] = temp_result[j] + twiddle_factorRow[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                            } else {
                                complex_result[numberHigh, numberCol, j] = -temp_result[j] * twiddle_factorRow[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                            }
                        }
                    }
                    for (int i = 0; i < resultRow; i++) {
                        result[numberHigh, numberCol, i] = (float)(complex_result[numberHigh, numberCol, i].Real / (resultCol * resultRow));
                    }
                }
            });
            return result;
        }
        public static float[,,,] ifft2(Complex[,,,] p1) {
            //準備輸入輸出data
            int p1Col = p1.GetLength(2);
            int p1Row = p1.GetLength(3);
            int ffttimesCol = (int)Math.Ceiling(Math.Log(p1Col, 2));
            int ffttimesRow = (int)Math.Ceiling(Math.Log(p1Row, 2));
            int resultKer = p1.GetLength(0);
            int resultHigh = p1.GetLength(1);
            int resultCol = (int)Math.Pow(2, ffttimesCol);
            int resultRow = (int)Math.Pow(2, ffttimesRow);
            Complex[,,,] complex_result = new Complex[resultKer, resultHigh, resultCol, resultRow];
            float[,,,] result = new float[resultKer, resultHigh, resultCol, resultRow];
            int[] change_addressCol = new int[resultCol];
            for (int i = 0; i < resultCol; i++) {
                int[] temp_address = new int[ffttimesCol];
                int temp_i = i;
                for (int j = 0; j < ffttimesCol; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int k = 0; k < ffttimesCol; k++) {
                    change_addressCol[i] += (int)(temp_address[k] * Math.Pow(2, k));
                }
            }
            Complex[] twiddle_factorCol = new Complex[resultCol / 2];
            for (int i = 0; i < resultCol / 2; i++) {
                twiddle_factorCol[i] = Complex.Exp(new Complex(0, 2 * Math.PI * ((float)i / resultCol)));
            }
            int[] change_addressRow = new int[resultRow];
            for (int i = 0; i < resultRow; i++) {
                int[] temp_address = new int[ffttimesRow];
                int temp_i = i;
                for (int j = 0; j < ffttimesRow; j++) {
                    temp_address[j] = temp_i % 2;
                    temp_i = temp_i / 2;
                }
                Array.Reverse(temp_address);
                for (int k = 0; k < ffttimesRow; k++) {
                    change_addressRow[i] += (int)(temp_address[k] * Math.Pow(2, k));
                }
            }
            Complex[] twiddle_factorRow = new Complex[resultRow / 2];
            for (int i = 0; i < resultRow / 2; i++) {
                twiddle_factorRow[i] = Complex.Exp(new Complex(0, 2 * Math.PI * ((float)i / resultRow)));
            }
            Parallel.For(0, resultKer, numberKer => {
                Parallel.For(0, resultHigh, numberHigh => {
                    //計算Col fft結果
                    for (int numberRow = 0; numberRow < p1Row; numberRow++) {
                        Complex[] temp_p1 = new Complex[resultCol];
                        for (int i = 0; i < p1Col; i++) {
                            temp_p1[i] = p1[numberKer, numberHigh, i, numberRow];
                        }
                        for (int i = 0; i < resultCol; i++) {
                            complex_result[numberKer, numberHigh, i, numberRow] = temp_p1[change_addressCol[i]];
                        }
                        for (int i = 0; i < ffttimesCol; i++) {
                            float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                            int temp_a2 = (int)Math.Pow(2, i);
                            int temp_a3 = (int)Math.Pow(2, i + 1);
                            int temp_a4 = resultCol / temp_a3;
                            Complex[] temp_result = new Complex[resultCol];
                            for (int k = 0; k < resultCol; k++) {
                                temp_result[k] = complex_result[numberKer, numberHigh, k, numberRow];
                            }
                            for (int j = 0; j < resultCol; j++) {
                                if ((j % temp_a3) < temp_a1) {
                                    complex_result[numberKer, numberHigh, j, numberRow] = temp_result[j] + twiddle_factorCol[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                                } else {
                                    complex_result[numberKer, numberHigh, j, numberRow] = -temp_result[j] * twiddle_factorCol[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                                }
                            }
                        }
                    }
                    //計算Row fft結果
                    for (int numberCol = 0; numberCol < resultCol; numberCol++) {
                        Complex[] temp_p1 = new Complex[resultRow];
                        for (int i = 0; i < p1Row; i++) {
                            temp_p1[i] = complex_result[numberKer, numberHigh, numberCol, i];
                        }
                        for (int i = 0; i < resultRow; i++) {
                            complex_result[numberKer, numberHigh, numberCol, i] = temp_p1[change_addressRow[i]];
                        }
                        for (int i = 0; i < ffttimesRow; i++) {
                            float temp_a1 = (float)((Math.Pow(2, i + 1) - 1) / 2);
                            int temp_a2 = (int)Math.Pow(2, i);
                            int temp_a3 = (int)Math.Pow(2, i + 1);
                            int temp_a4 = resultRow / temp_a3;
                            Complex[] temp_result = new Complex[resultRow];
                            for (int k = 0; k < resultRow; k++) {
                                temp_result[k] = complex_result[numberKer, numberHigh, numberCol, k];
                            }
                            for (int j = 0; j < resultRow; j++) {
                                if ((j % temp_a3) < temp_a1) {
                                    complex_result[numberKer, numberHigh, numberCol, j] = temp_result[j] + twiddle_factorRow[temp_a4 * (j % temp_a2)] * temp_result[j + temp_a2];
                                } else {
                                    complex_result[numberKer, numberHigh, numberCol, j] = -temp_result[j] * twiddle_factorRow[temp_a4 * (j % temp_a2)] + temp_result[j - temp_a2];
                                }
                            }
                        }
                        for (int i = 0; i < resultRow; i++) {
                            result[numberKer, numberHigh, numberCol, i] = (float)(complex_result[numberKer, numberHigh, numberCol, i].Real / (resultCol * resultRow));
                        }
                    }
                });
            });
            return result;
        }
    }
}
