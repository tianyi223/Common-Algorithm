// FFT Defines

#ifndef _FFT_H_JXJ_201503041021_
#define _FFT_H_JXJ_201503041021_



//* ��������:  FFT����, �Ӿ������ Select one in double or float.
typedef double TFFTData;

//* ��������:  ����(Complex)
typedef struct _TFFTComp{
    TFFTData real;
    TFFTData imag;
}TFFTComp;//* ��������:  ����

//* FFT:   x[]-ԭʼ����  ; dft-Descrete Fourier Transform  ; n-ԭʼ���ݳ���, Ҫ��Ϊ2��������2^m ;
void fft(const TFFTComp x[],TFFTComp dft[],const unsigned long n);



#endif