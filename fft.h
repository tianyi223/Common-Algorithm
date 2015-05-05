// FFT Defines

#ifndef _FFT_H_JXJ_201503041021_
#define _FFT_H_JXJ_201503041021_



//* 定义类型:  FFT数据, 视具体情况 Select one in double or float.
typedef double TFFTData;

//* 定义类型:  复数(Complex)
typedef struct _TFFTComp{
    TFFTData real;
    TFFTData imag;
}TFFTComp;//* 定义类型:  复数

//* FFT:   x[]-原始数据  ; dft-Descrete Fourier Transform  ; n-原始数据长度, 要求为2的整数幂2^m ;
void fft(const TFFTComp x[],TFFTComp dft[],const unsigned long n);



#endif