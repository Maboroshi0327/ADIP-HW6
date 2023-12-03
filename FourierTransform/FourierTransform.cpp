#include "FourierTransform.hpp"

#include <iostream>

namespace dft
{
    Complex1D fft1D(const Complex1D& data, bool inverse)
    {
        int N = (int)data.size();
        Complex1D result;

        if (1 << (int)std::log2(N) == N)
            result = radix2FFT1D(data, inverse);
        else
            result = BluesteinFFT(data, inverse);

        double invRootN = 1.0 / std::sqrt(data.size());
        for (auto& i : result)
            i *= invRootN;

        return result;
    }

    Complex2D fft2D(const Complex2D& data, bool shift, bool inverse)
    {
        int rows = (int)data.size();
        int cols = (int)data[0].size();
        Complex2D result = data;

        // Shifting in Frequency
        if (inverse == false && shift == true)
        {
            int shiftR = 1;
            for (int r = 0; r < rows; r++)
            {
                int shiftC = 1;
                for (int c = 0; c < cols; c++)
                {
                    result[r][c] *= shiftR * shiftC;
                    shiftC *= -1;
                }
                shiftR *= -1;
            }
        }

        // FFT rows
        bool flag = (1 << (int)std::log2(rows) == rows) ? true : false;
        for (auto& row : result)
        {
            if (flag)
                row = radix2FFT1D(row, inverse);
            else
                row = BluesteinFFT(row, inverse);
        }

        // FFT columns
        flag = (1 << (int)std::log2(cols) == cols) ? true : false;
        for (int c = 0; c < cols; c++)
        {
            Complex1D col(rows);
            for (int r = 0; r < rows; r++)
                col[r] = result[r][c];

            if (flag)
                col = radix2FFT1D(col, inverse);
            else
                col = BluesteinFFT(col, inverse);

            for (int r = 0; r < rows; r++)
                result[r][c] = col[r];
        }

        // Shifting in Frequency
        if (inverse == true && shift == true)
        {
            int shiftR = 1;
            for (int r = 0; r < rows; r++)
            {
                int shiftC = 1;
                for (int c = 0; c < cols; c++)
                {
                    result[r][c] *= shiftR * shiftC;
                    shiftC *= -1;
                }
                shiftR *= -1;
            }
        }

        double invRootMN = 1.0 / std::sqrt(rows * cols);
        for (auto& i : result)
            for (auto& j : i)
                j *= invRootMN;

        return result;
    }

    void fft2D(const Complex2D& data, Complex2D& output, bool shift, bool inverse)
    {
        int rows = (int)data.size();
        int cols = (int)data[0].size();
        Complex2D result = data;

        // Shifting in Frequency
        if (inverse == false && shift == true)
        {
            int shiftR = 1;
            for (int r = 0; r < rows; r++)
            {
                int shiftC = 1;
                for (int c = 0; c < cols; c++)
                {
                    result[r][c] *= shiftR * shiftC;
                    shiftC *= -1;
                }
                shiftR *= -1;
            }
        }

        // FFT rows
        bool flag = (1 << (int)std::log2(rows) == rows) ? true : false;
        for (auto& row : result)
        {
            if (flag)
                row = radix2FFT1D(row, inverse);
            else
                row = BluesteinFFT(row, inverse);
        }

        // FFT columns
        flag = (1 << (int)std::log2(cols) == cols) ? true : false;
        for (int c = 0; c < cols; c++)
        {
            Complex1D col(rows);
            for (int r = 0; r < rows; r++)
                col[r] = result[r][c];

            if (flag)
                col = radix2FFT1D(col, inverse);
            else
                col = BluesteinFFT(col, inverse);

            for (int r = 0; r < rows; r++)
                result[r][c] = col[r];
        }

        // Shifting in Frequency
        if (inverse == true && shift == true)
        {
            int shiftR = 1;
            for (int r = 0; r < rows; r++)
            {
                int shiftC = 1;
                for (int c = 0; c < cols; c++)
                {
                    result[r][c] *= shiftR * shiftC;
                    shiftC *= -1;
                }
                shiftR *= -1;
            }
        }

        double invRootMN = 1.0 / std::sqrt(rows * cols);
        for (auto& i : result)
            for (auto& j : i)
                j *= invRootMN;

        output = result;
    }

    Complex1D radix2FFT1D(const Complex1D& data, bool inverse)
    {
        int N = (int)data.size();
        if (1 == N)
            return data;

        Complex1D even(N / 2);
        Complex1D odd(N / 2);
        for (int i = 0; i < N / 2; i++)
        {
            even[i] = data[i * 2];
            odd[i] = data[i * 2 + 1];
        }

        even = radix2FFT1D(even, inverse);
        odd = radix2FFT1D(odd, inverse);

        Complex Wn;
        if (inverse == false)
            Wn = std::polar(1.0, -2.0 * M_PI / N);
        else
            Wn = std::polar(1.0, 2.0 * M_PI / N);

        Complex W(1, 0), WnOdd;
        Complex1D result(N);
        for (int i = 0; i < N / 2; i++)
        {
            WnOdd = W * odd[i];
            result[i] = even[i] + WnOdd;
            result[i + N / 2] = even[i] - WnOdd;
            W = W * Wn;
        }

        return result;
    }

    Complex1D BluesteinFFT(const Complex1D& data, bool inverse)
    {
        int N = (int)data.size();
        if (1 == N)
            return data;

        double inv = (inverse == false) ? -2.0 : 2.0;
        Complex w = std::polar(1.0, inv * M_PI / N);

        // Pad zeros to the next power of 2
        int N1 = 2 * N - 1;
        int N2 = (int)std::pow(2, std::ceil(std::log2(N1)));

        // Make A and B arrays for convolution
        Complex1D A(N2, Complex(0.0, 0.0)), B(N2, Complex(0.0, 0.0)), chirp(N);

        chirp[0] = Complex(1.0, 0.0);
        A[0] = data[0] * chirp[0];
        B[0] = Complex(1.0, 0.0);
        for (int n = 1; n < N; n++)
        {
            chirp[n] = std::pow(w, n * n / 2.0);
            A[n] = data[n] * chirp[n];
            B[n] = B[N2 - n] = 1.0 / chirp[n];
        }

        // convolution
        A = radix2FFT1D(A, false);
        B = radix2FFT1D(B, false);
        Complex1D conv(N2);
        for (int n = 0; n < N2; n++)
            conv[n] = A[n] * B[n];
        conv = radix2FFT1D(conv, true);

        Complex1D result(N);
        for (int k = 0; k < N; k++)
            result[k] = conv[k] * chirp[k] / (double)N2;

        return result;
    }

    Complex2D Filtering(const Complex2D& fftArray, const Double2D& filter)
    {
        int M = (int)fftArray.size();
        int N = (int)fftArray[0].size();

        Complex2D filtered(M, Complex1D(N));
        for (int u = 0; u < M; u++)
            for (int v = 0; v < N; v++)
            {
                filtered[u][v] = fftArray[u][v] * filter[u][v];
            }

        return filtered;
    }

    Complex2D Filtering(const Complex2D& fftArray, const Complex2D& filter)
    {
        int M = (int)fftArray.size();
        int N = (int)fftArray[0].size();

        Complex2D filtered(M, Complex1D(N));
        for (int u = 0; u < M; u++)
            for (int v = 0; v < N; v++)
            {
                filtered[u][v] = fftArray[u][v] * filter[u][v];
            }

        return filtered;
    }

    Complex2D InverseFiltering(const Complex2D& fftArray, const Complex2D& filter, double radius)
    {
        int M = (int)fftArray.size();
        int N = (int)fftArray[0].size();

        Complex2D filtered(M, Complex1D(N));
        for (int u = 0; u < M; u++)
            for (int v = 0; v < N; v++)
            {
                if (radius >= Duv(u, v, M, N))
                    filtered[u][v] = fftArray[u][v] / filter[u][v];
                else
                    filtered[u][v] = fftArray[u][v];
            }

        return filtered;
    }

    Complex2D WienerFiltering(const Complex2D& fftArray, const Complex2D& filter, double K)
    {
        int M = (int)fftArray.size();
        int N = (int)fftArray[0].size();

        Complex2D filtered(M, Complex1D(N));
        for (int u = 0; u < M; u++)
            for (int v = 0; v < N; v++)
            {
                Complex H = filter[u][v];
                double absH2 = std::pow(std::abs(H), 2);
                filtered[u][v] = (1.0 / H) * (absH2 / (absH2 + K)) * fftArray[u][v];
            }

        return filtered;
    }

    Complex2D ConstrainedLeastSquareFilter(const Complex2D& fftArray, const Complex2D& filter, const Complex2D& Laplacian, double gamma)
    {
        int M = (int)fftArray.size();
        int N = (int)fftArray[0].size();

        Complex2D filtered(M, Complex1D(N));
        for (int u = 0; u < M; u++)
            for (int v = 0; v < N; v++)
            {
                Complex H_conj = std::conj(filter[u][v]);
                double absH2 = std::pow(std::abs(H_conj), 2);
                double gammaP = gamma * std::pow(std::abs(Laplacian[u][v]), 2);
                filtered[u][v] = (H_conj / (absH2 + gammaP)) * fftArray[u][v];
            }

        return filtered;
    }

    Byte2D drawFFT(const Complex2D& data)
    {
        int rows = (int)data.size();
        int cols = (int)data[0].size();

        Double2D powMagnitude(rows, Double1D(cols, 0));
        double max = std::pow(std::abs(data[0][0]), 0.2);
        double min = max;
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols; c++)
            {
                double pow = std::pow(std::abs(data[r][c]), 0.2);
                powMagnitude[r][c] = pow;
                if (max < pow)
                    max = pow;
                if (min > pow)
                    min = pow;
            }

        Byte2D result(rows, Byte1D(cols, 0));
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols; c++)
            {
                result[r][c] = (uint8_t)std::round((powMagnitude[r][c] - min) / (max - min) * 255);
            }

        return result;
    }

    double Duv(double u, double v, double M, double N)
    {
        double r = u - M * 0.5;
        double c = v - N * 0.5;

        return std::sqrt(r * r + c * c);
    }

    Double2D makeGaussianLPF(double D0, int M, int N)
    {
        Double2D GLPF(M, Double1D(N));

        double M_f = M, N_f = N;
        for (int u = 0; u < M; u++)
            for (int v = 0; v < N; v++)
            {
                double D = Duv(u, v, M_f, N_f);
                GLPF[u][v] = std::exp(-1 * D * D / (2 * D0 * D0));
            }

        return GLPF;
    }

    Double2D makeGaussianHPF(double D0, int M, int N)
    {
        Double2D GHPF(M, Double1D(N));

        double M_f = M, N_f = N;
        for (int u = 0; u < M; u++)
            for (int v = 0; v < N; v++)
            {
                double D = Duv(u, v, M_f, N_f);
                GHPF[u][v] = 1.0 - std::exp(-1 * D * D / (2 * D0 * D0));
            }

        return GHPF;
    }

    Double2D makeButterworthLPF(double D0, double n, int M, int N)
    {
        Double2D BLPF(M, Double1D(N));

        double M_f = M, N_f = N;
        for (int u = 0; u < M; u++)
            for (int v = 0; v < N; v++)
            {
                double D = Duv(u, v, M_f, N_f);
                BLPF[u][v] = 1.0 / (1.0 + std::pow(D / D0, 2 * n));
            }

        return BLPF;
    }

    Double2D makeButterworthHPF(double D0, double n, int M, int N)
    {
        Double2D BHPF(M, Double1D(N));

        double M_f = M, N_f = N;
        for (int u = 0; u < M; u++)
            for (int v = 0; v < N; v++)
            {
                double D = Duv(u, v, M_f, N_f);
                BHPF[u][v] = 1.0 - 1.0 / (1.0 + std::pow(D / D0, 2 * n));
            }

        return BHPF;
    }

    Double2D makeIdealLPF(double D0, int M, int N)
    {
        Double2D ILPF(M, Double1D(N));

        double M_f = M, N_f = N;
        for (int u = 0; u < M; u++)
            for (int v = 0; v < N; v++)
            {
                double D = Duv(u, v, M_f, N_f);
                ILPF[u][v] = (D <= D0) ? 1.0 : 0.0;
            }

        return ILPF;
    }

    Double2D makeIdealHPF(double D0, int M, int N)
    {
        Double2D IHPF(M, Double1D(N));

        double M_f = M, N_f = N;
        for (int u = 0; u < M; u++)
            for (int v = 0; v < N; v++)
            {
                double D = Duv(u, v, M_f, N_f);
                IHPF[u][v] = (D <= D0) ? 0.0 : 1.0;
            }

        return IHPF;
    }

    Complex2D makeSobelFilter(int angle, int M, int N)
    {
        Complex2D sobel(M, Complex1D(N, Complex(0.0, 0.0)));
        if (angle == 0)
        {
            sobel[0][0] = Complex(-1.0, 0.0); sobel[0][1] = Complex(0.0, 0.0); sobel[0][2] = Complex(1.0, 0.0);
            sobel[1][0] = Complex(-2.0, 0.0); sobel[1][1] = Complex(0.0, 0.0); sobel[1][2] = Complex(2.0, 0.0);
            sobel[2][0] = Complex(-1.0, 0.0); sobel[2][1] = Complex(0.0, 0.0); sobel[2][2] = Complex(1.0, 0.0);
        }
        else if (angle == 90)
        {
            sobel[0][0] = Complex(1.0, 0.0);  sobel[0][1] = Complex(2.0, 0.0);  sobel[0][2] = Complex(1.0, 0.0);
            sobel[1][0] = Complex(0.0, 0.0);  sobel[1][1] = Complex(0.0, 0.0);  sobel[1][2] = Complex(0.0, 0.0);
            sobel[2][0] = Complex(-1.0, 0.0); sobel[2][1] = Complex(-2.0, 0.0); sobel[2][2] = Complex(-1.0, 0.0);
        }
        else
            return sobel;

        sobel = fft2D(sobel, true, false);

        return sobel;
    }

    void makeSobelFilter(int angle, int M, int N, Complex2D& output)
    {
        Complex2D sobel(M, Complex1D(N, Complex(0.0, 0.0)));
        if (angle == 0)
        {
            sobel[0][0] = Complex(-1.0, 0.0); sobel[0][1] = Complex(0.0, 0.0); sobel[0][2] = Complex(1.0, 0.0);
            sobel[1][0] = Complex(-2.0, 0.0); sobel[1][1] = Complex(0.0, 0.0); sobel[1][2] = Complex(2.0, 0.0);
            sobel[2][0] = Complex(-1.0, 0.0); sobel[2][1] = Complex(0.0, 0.0); sobel[2][2] = Complex(1.0, 0.0);
        }
        else if (angle == 90)
        {
            sobel[0][0] = Complex(1.0, 0.0);  sobel[0][1] = Complex(2.0, 0.0);  sobel[0][2] = Complex(1.0, 0.0);
            sobel[1][0] = Complex(0.0, 0.0);  sobel[1][1] = Complex(0.0, 0.0);  sobel[1][2] = Complex(0.0, 0.0);
            sobel[2][0] = Complex(-1.0, 0.0); sobel[2][1] = Complex(-2.0, 0.0); sobel[2][2] = Complex(-1.0, 0.0);
        }
        else
            output = sobel;

        sobel = fft2D(sobel, true, false);

        output = sobel;
    }

    Complex2D makeUniformLinearMotionFilter(double T, double a, double b, int M, int N)
    {
        Complex2D ULMF(M, Complex1D(N));
        for (int u = 0; u < M; u++)
            for (int v = 0; v < N; v++)
            {
                double uu = u - M * 0.5;
                double vv = v - N * 0.5;
                double pi_ua_vb = M_PI * (uu * a + vv * b);
                if (pi_ua_vb == 0.0)
                    ULMF[u][v] = 1.0;
                else
                    ULMF[u][v] = T / pi_ua_vb * std::sin(pi_ua_vb) * std::polar(1.0, -1 * pi_ua_vb);
            }

        return ULMF;
    }

    Complex2D makeLaplacianFilter(int neighbor, int M, int N)
    {
        Complex2D laplacian(M, Complex1D(N, Complex(0.0, 0.0)));

        if (neighbor == 4)
        {
            laplacian[0][0] = 0;  laplacian[0][1] = -1; laplacian[0][2] = 0;
            laplacian[1][0] = -1; laplacian[1][1] = 4;  laplacian[1][2] = -1;
            laplacian[2][0] = 0;  laplacian[2][1] = -1; laplacian[2][2] = 0;
        }

        if (neighbor == 8)
        {
            laplacian[0][0] = -1; laplacian[0][1] = -1; laplacian[0][2] = -1;
            laplacian[1][0] = -1; laplacian[1][1] = 8;  laplacian[1][2] = -1;
            laplacian[2][0] = -1; laplacian[2][1] = -1; laplacian[2][2] = -1;
        }

        laplacian = fft2D(laplacian, true, false);

        return laplacian;
    }
}

namespace dct
{
    Double1D fct1D(const Double1D& data, bool inverse)
    {
        return radix2FCT1D(data, inverse);
    }

    Double2D fct2D(const Double2D& data, bool inverse)
    {
        return radix2FCT2D(data, inverse);
    }

    Double1D radix2FCT1D(const Double1D& data, bool inverse)
    {
        int N = (int)data.size();
        Complex1D Y(N);
        if (inverse == false)
            for (int k = 0; k < N / 2; k++)
            {
                Y[k] = data[2 * k];
                Y[N - 1 - k] = data[2 * k + 1];
            }
        else
        {
            Y[0] = data[0] * std::sqrt(1.0 / 2.0);
            for (int r = 1; r < N; r++)
                Y[r] = data[r] * std::polar(1.0, M_PI * r / 2.0 / N);
        }

        // FFT
        Y = dft::radix2FFT1D(Y, false);

        Double1D result(N);
        if (inverse == false)
        {
            result[0] = Y[0].real();
            for (int r = 1; r <= N / 2; r++)
            {
                Complex H = std::polar(1.0, M_PI * r / 2.0 / N) * Y[r];
                result[r] = H.real();
                result[N - r] = H.imag();
            }
        }
        else
        {
            for (int k = 0; k < N / 2; k++)
            {
                result[2 * k] = Y[k].real();
                result[2 * k + 1] = Y[N - 1 - k].real();
            }
        }

        if (inverse == false)
        {
            result[0] *= (2 * std::sqrt(1.0 / 2.0) / N);
            for (int r = 1; r < N; r++)
                result[r] *= (2.0 / N);
        }

        return result;
    }

    Double2D radix2FCT2D(const Double2D& data, bool inverse)
    {
        int rows = (int)data.size();
        int cols = (int)data[0].size();
        Double2D result = data;

        // FCT rows and colums
        for (auto& row : result)
            row = radix2FCT1D(row, inverse);
        for (int c = 0; c < cols; c++)
        {
            Double1D col(rows);
            for (int r = 0; r < rows; r++)
                col[r] = result[r][c];
            col = radix2FCT1D(col, inverse);
            for (int r = 0; r < rows; r++)
                result[r][c] = col[r];
        }

        return result;
    }

    Byte2D drawFCT(const Double2D& data)
    {
        int rows = (int)data.size();
        int cols = (int)data[0].size();

        Double2D powMagnitude(rows, Double1D(cols, 0));
        double max = std::pow(std::abs(data[0][0]), 0.2);
        double min = max;
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols; c++)
            {
                double pow = std::pow(std::abs(data[r][c]), 0.2);
                powMagnitude[r][c] = pow;
                if (max < pow)
                    max = pow;
                if (min > pow)
                    min = pow;
            }

        Byte2D result(rows, Byte1D(cols, 0));
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols; c++)
            {
                result[r][c] = (uint8_t)std::round((powMagnitude[r][c] - min) / (max - min) * 255);
            }

        return result;
    }
}
