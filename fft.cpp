
#include "math.h"
#include "fft.h"

#define PI ((TFFTData)3.14159265358979324) //定义圆周率值

//* FFT:   x[]-原始数据  ; dft-Descrete Fourier Transform  ; n-原始数据长度, 要求为2的整数幂2^m ;
void fft(const TFFTComp x[],TFFTComp dft[],const unsigned long n)
{
    unsigned long ix1 = 0, ix2 = 0;    // x的绝对索引值
    unsigned int  im  = 0, m   = 0; // n == 2^m
    unsigned long k = 0;
    unsigned long ngroup = 0; // 当前级蝶形的组数 = 2^(n-im-1)
    unsigned long igroup = 0; // 第igroup组
    unsigned long nnow = 0;   // 当前每组蝶形 包含元素个数
    TFFTComp wnk = {0}, wnkx2k = {0};
    TFFTData angn = (TFFTData)0;

    // 计算级数 若log2(n)不为整数, 则向上取整
    for (m = 0; (1<<m) < n ; m++) { ; }

    // 准备倒位序
    for (ix1 = 0; ix1 < n ; ix1++)
    {
        for(im = 0,ix2 = 0; im < m ; im++) // 计算ix1倒位值,如: 0101->1010
            { ix2 |= (((ix1>>im)&0X00000001)<<(m-im-1)); }
        
        dft[ix2] = x[ix1];
    }

    // 蝶形算法 计算FFT
    for (im = 0; im < m ; im++) // im为当前的级数 一共log2(n)级
    {
        nnow   = (1<<(im+1));     // 当前每组蝶形 包含元素个数
        ngroup = (1<<(m-im-1)); // 蝶形组数
        
        angn = -2*PI/((TFFTData)nnow); // 计算 -2pi/N, 提示 WNk = exp(-j*2pi/N*k), {欧拉公式: exp(-j*x)=cos(x)+j*sin(x)}
        for (igroup = 0; igroup < ngroup ; igroup++)
        {
            for (k = 0; k < (nnow>>1) ; k++) // k为当前级处理的位置
            {
                wnk.real = cos(angn*((TFFTData)k)); // 计算WNk = exp(-j*2pi/N * k),  {欧拉公式: exp(-j*x)=cos(x)+j*sin(x)}
                wnk.imag = sin(angn*((TFFTData)k));
                
                ix1 = igroup*nnow + k; // X1(k)的绝对索引值
                ix2 = ix1 + (nnow>>1); // X2(k)的绝对索引值
                
                wnkx2k.real = wnk.real*dft[ix2].real - wnk.imag*dft[ix2].imag; // WNk*X2(k) 实部
                wnkx2k.imag = wnk.real*dft[ix2].imag + wnk.imag*dft[ix2].real; // WNk*X2(k) 虚部
                
                dft[ix2].real = dft[ix1].real - wnkx2k.real; // 计算 X1(k) - WNk*X2(k)
                dft[ix2].imag = dft[ix1].imag - wnkx2k.imag;

                dft[ix1].real += wnkx2k.real; // 计算 X1(k) + WNk*X2(k)
                dft[ix1].imag += wnkx2k.imag;
            }
        }
    }
}