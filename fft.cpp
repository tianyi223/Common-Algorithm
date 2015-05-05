
#include "math.h"
#include "fft.h"

#define PI ((TFFTData)3.14159265358979324) //����Բ����ֵ

//* FFT:   x[]-ԭʼ����  ; dft-Descrete Fourier Transform  ; n-ԭʼ���ݳ���, Ҫ��Ϊ2��������2^m ;
void fft(const TFFTComp x[],TFFTComp dft[],const unsigned long n)
{
    unsigned long ix1 = 0, ix2 = 0;    // x�ľ�������ֵ
    unsigned int  im  = 0, m   = 0; // n == 2^m
    unsigned long k = 0;
    unsigned long ngroup = 0; // ��ǰ�����ε����� = 2^(n-im-1)
    unsigned long igroup = 0; // ��igroup��
    unsigned long nnow = 0;   // ��ǰÿ����� ����Ԫ�ظ���
    TFFTComp wnk = {0}, wnkx2k = {0};
    TFFTData angn = (TFFTData)0;

    // ���㼶�� ��log2(n)��Ϊ����, ������ȡ��
    for (m = 0; (1<<m) < n ; m++) { ; }

    // ׼����λ��
    for (ix1 = 0; ix1 < n ; ix1++)
    {
        for(im = 0,ix2 = 0; im < m ; im++) // ����ix1��λֵ,��: 0101->1010
            { ix2 |= (((ix1>>im)&0X00000001)<<(m-im-1)); }
        
        dft[ix2] = x[ix1];
    }

    // �����㷨 ����FFT
    for (im = 0; im < m ; im++) // imΪ��ǰ�ļ��� һ��log2(n)��
    {
        nnow   = (1<<(im+1));     // ��ǰÿ����� ����Ԫ�ظ���
        ngroup = (1<<(m-im-1)); // ��������
        
        angn = -2*PI/((TFFTData)nnow); // ���� -2pi/N, ��ʾ WNk = exp(-j*2pi/N*k), {ŷ����ʽ: exp(-j*x)=cos(x)+j*sin(x)}
        for (igroup = 0; igroup < ngroup ; igroup++)
        {
            for (k = 0; k < (nnow>>1) ; k++) // kΪ��ǰ�������λ��
            {
                wnk.real = cos(angn*((TFFTData)k)); // ����WNk = exp(-j*2pi/N * k),  {ŷ����ʽ: exp(-j*x)=cos(x)+j*sin(x)}
                wnk.imag = sin(angn*((TFFTData)k));
                
                ix1 = igroup*nnow + k; // X1(k)�ľ�������ֵ
                ix2 = ix1 + (nnow>>1); // X2(k)�ľ�������ֵ
                
                wnkx2k.real = wnk.real*dft[ix2].real - wnk.imag*dft[ix2].imag; // WNk*X2(k) ʵ��
                wnkx2k.imag = wnk.real*dft[ix2].imag + wnk.imag*dft[ix2].real; // WNk*X2(k) �鲿
                
                dft[ix2].real = dft[ix1].real - wnkx2k.real; // ���� X1(k) - WNk*X2(k)
                dft[ix2].imag = dft[ix1].imag - wnkx2k.imag;

                dft[ix1].real += wnkx2k.real; // ���� X1(k) + WNk*X2(k)
                dft[ix1].imag += wnkx2k.imag;
            }
        }
    }
}