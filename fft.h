//

#ifndef _FFT_H_JXJ_201503041021_
#define _FFT_H_JXJ_201503041021_



//* ��������:  FFT����, �Ӿ������ ��ѡ��double,float
typedef double TFFTData;

//* ��������:  ����
typedef struct _TFFTComp{
    TFFTData real;
    TFFTData imag;
}TFFTComp;//* ��������:  ����

//* FFT:   x[]-ԭʼ����  ; dft-Descrete Fourier Transform  ; n-ԭʼ���ݳ���, Ҫ��Ϊ2��������2^m ;
void fft(const TFFTComp x[],TFFTComp dft[],const unsigned long n);




#endif