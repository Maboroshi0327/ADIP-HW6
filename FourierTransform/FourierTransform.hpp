#ifndef _FourierTransform_
#define _FourierTransform_

#include <cmath>
#include <cstdint>

#include <complex>
#include <vector>

typedef std::complex<double> Complex;
typedef std::vector<Complex> Complex1D;
typedef std::vector<Complex1D> Complex2D;

typedef std::vector<double> Double1D;
typedef std::vector<Double1D> Double2D;

typedef uint8_t Byte;
typedef std::vector<Byte> Byte1D;
typedef std::vector<Byte1D> Byte2D;

template<typename T>
std::vector<std::vector<T>> minMaxNormalization(const Double2D& data, double min_v, double max_v);
template<typename T>
std::vector<std::vector<T>> minMaxNormalization(const Complex2D& data, double min_v, double max_v);
template<typename T>
Complex2D vector2Complex2D(const std::vector<std::vector<T>>& data);
template<typename T>
Double2D vector2Double2D(const std::vector<std::vector<T>>& data);
template<typename T>
std::vector<std::vector<T>> complex2dReal2vector(const Complex2D& data);
template<typename T>
std::vector<std::vector<T>> complex2dAbs2vector(const Complex2D& data);
template<typename T>
std::vector<std::vector<T>> double2d2vector(const Double2D& data);

namespace dft
{
    Complex1D fft1D(const Complex1D& data, bool inverse);
    Complex2D fft2D(const Complex2D& data, bool shift, bool inverse);
    void fft2D(const Complex2D& data, Complex2D& output, bool shift, bool inverse);

    /**
     * @brief Cooley¡VTukey FFT algorithm
     */
    Complex1D radix2FFT1D(const Complex1D& data, bool inverse);

    /**
     * @brief Bluestein's algorithm(Non-Power-of-Two FFT)
     * @brief Reference: https://en.wikipedia.org/wiki/Chirp_Z-transform
     */
    Complex1D BluesteinFFT(const Complex1D& data, bool inverse);

    Complex2D Filtering(const Complex2D& fftArray, const Double2D& filter);
    Complex2D Filtering(const Complex2D& fftArray, const Complex2D& filter);
    Complex2D InverseFiltering(const Complex2D& fftArray, const Complex2D& filter, double radius);
    Complex2D WienerFiltering(const Complex2D& fftArray, const Complex2D& filter, double K);
    Complex2D ConstrainedLeastSquareFilter(const Complex2D& fftArray, const Complex2D& filter, const Complex2D& Laplacian, double gamma);
    Byte2D drawFFT(const Complex2D& data);

    /**
     * @brief Calculate D(u, v)
     * @param u: Frequency domain coordinate(row)
     * @param v: Frequency domain coordinate(column)
     * @param M: Output height(rows)
     * @param N: Output width(columns)
     * @return D(u, v)
     */
    double Duv(double u, double v, double M, double N);

    Double2D makeGaussianLPF(double D0, int M, int N);
    Double2D makeGaussianHPF(double D0, int M, int N);

    /**
     * @brief Butterworth Low-pass filter (Frequency-Domain)
     * @param D0: Cutoff frequency
     * @param n: Orders
     * @param M: Output height(rows)
     * @param N: Output width(columns)
     * @return BLPF(row, col)
     */
    Double2D makeButterworthLPF(double D0, double n, int M, int N);
    Double2D makeButterworthHPF(double D0, double n, int M, int N);

    Double2D makeIdealLPF(double D0, int M, int N);
    Double2D makeIdealHPF(double D0, int M, int N);
    Complex2D makeSobelFilter(int angle, int M, int N);
    void makeSobelFilter(int angle, int M, int N, Complex2D& output);

    /**
     * @brief T/pi(ua+vb) * sin(pi(ua+vb)) * e^(-j*pi(ua+vb))
     * @param T: Period
     * @param a: speed in x direction
     * @param b: speed in y direction
     * @param M: Output height(rows)
     * @param N: Output width(columns)
     * @return ULMF(row, col)
     */
    Complex2D makeUniformLinearMotionFilter(double T, double a, double b, int M, int N);

    Complex2D makeLaplacianFilter(int neighbor, int M, int N);
}

namespace dct
{
    Double1D fct1D(const Double1D& data, bool inverse);
    Double2D fct2D(const Double2D& data, bool inverse);

    /**
     * @brief Reference:
     * @brief M. Narasimha and A. Peterson, 'On the computation of the discrete cosine transform',
     * IEEE Transactions on Communications, vol. 26, no. 6, pp. 934¡V936, 1978.
     */
    Double1D radix2FCT1D(const Double1D& data, bool inverse);
    Double2D radix2FCT2D(const Double2D& data, bool inverse);
    Byte2D drawFCT(const Double2D& data);
}

#include "FourierTransform.ipp"

#endif // _FourierTransform_