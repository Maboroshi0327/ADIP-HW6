#include "HW_2.h"

#include <sstream>

#include "../ImageProcessing/ImageProcessing.hpp"
#include "../FourierTransform/FourierTransform.hpp"

void HW_2_1()
{
    dip::Image lena_256, lena_256_blur, lena_256_blur_noise, img;
    lena_256.openRAW("lena_256.raw", dip::Size(256, 256));
    lena_256_blur.openRAW("lena_256_blur.raw", dip::Size(256, 256));
    lena_256_blur_noise.openRAW("lena_256_blur_noise.raw", dip::Size(256, 256));

    Complex2D ULMF = dft::makeUniformLinearMotionFilter(1.0, -0.07, 0.07, 256, 256);

    Complex2D lena_FFT = vector2Complex2D(lena_256_blur.img);
    lena_FFT = dft::fft2D(lena_FFT, true, false);
    for (int i = 10; i <= 110; i += 20)
    {
        if (i == 110)
            i = 1000;

        Complex2D lena_256_invF = dft::InverseFiltering(lena_FFT, ULMF, i);
        lena_256_invF = dft::fft2D(lena_256_invF, true, true);
        Byte2D result = minMaxNormalization<Byte>(lena_256_invF, 0, 255);

        std::stringstream fileName;
        fileName << "HW_2_1_" << i << ".raw";
        img.setImg(result);
        img.saveRAW(fileName.str().c_str(), true);

        double mse = dip::mse(lena_256, img);
        double psnr = dip::psnr(lena_256, img);
        std::cout << "MSE: " << mse << std::endl;
        std::cout << "PSNR: " << psnr << std::endl;
    }

    lena_FFT = vector2Complex2D(lena_256_blur_noise.img);
    lena_FFT = dft::fft2D(lena_FFT, true, false);
    for (int i = 10; i <= 110; i += 20)
    {
        if (i == 110)
            i = 1000;

        Complex2D lena_256_invF = dft::InverseFiltering(lena_FFT, ULMF, i);
        lena_256_invF = dft::fft2D(lena_256_invF, true, true);
        Byte2D result = minMaxNormalization<Byte>(lena_256_invF, 0, 255);

        std::stringstream fileName;
        fileName << "HW_2_1_noise_" << i << ".raw";
        img.setImg(result);
        img.saveRAW(fileName.str().c_str(), true);

        double mse = dip::mse(lena_256, img);
        double psnr = dip::psnr(lena_256, img);
        std::cout << "MSE: " << mse << std::endl;
        std::cout << "PSNR: " << psnr << std::endl;
    }
}

void HW_2_2()
{
    dip::Image lena_256, lena_256_blur, lena_256_blur_noise, img;
    lena_256.openRAW("lena_256.raw", dip::Size(256, 256));
    lena_256_blur.openRAW("lena_256_blur.raw", dip::Size(256, 256));
    lena_256_blur_noise.openRAW("lena_256_blur_noise.raw", dip::Size(256, 256));

    Complex2D ULMF = dft::makeUniformLinearMotionFilter(1.0, -0.07, 0.07, 256, 256);

    Complex2D lena_FFT = vector2Complex2D(lena_256_blur.img);
    lena_FFT = dft::fft2D(lena_FFT, true, false);
    for (double k = 0.0001; k <= 0.1; k *= 10.0)
    {
        Complex2D lena_256_invF = dft::WienerFiltering(lena_FFT, ULMF, k);
        lena_256_invF = dft::fft2D(lena_256_invF, true, true);
        Byte2D result = minMaxNormalization<Byte>(lena_256_invF, 0, 255);

        std::stringstream fileName;
        fileName << "HW_2_2_" << k << ".raw";
        img.setImg(result);
        img.saveRAW(fileName.str().c_str(), true);

        double mse = dip::mse(lena_256, img);
        double psnr = dip::psnr(lena_256, img);
        std::cout << "MSE: " << mse << std::endl;
        std::cout << "PSNR: " << psnr << std::endl;
    }

    lena_FFT = vector2Complex2D(lena_256_blur_noise.img);
    lena_FFT = dft::fft2D(lena_FFT, true, false);
    for (double k = 0.0001; k <= 0.1; k *= 10.0)
    {
        Complex2D lena_256_invF = dft::WienerFiltering(lena_FFT, ULMF, k);
        lena_256_invF = dft::fft2D(lena_256_invF, true, true);
        Byte2D result = minMaxNormalization<Byte>(lena_256_invF, 0, 255);

        std::stringstream fileName;
        fileName << "HW_2_2_noise_" << k << ".raw";
        img.setImg(result);
        img.saveRAW(fileName.str().c_str(), true);

        double mse = dip::mse(lena_256, img);
        double psnr = dip::psnr(lena_256, img);
        std::cout << "MSE: " << mse << std::endl;
        std::cout << "PSNR: " << psnr << std::endl;
    }
}

void HW_2_3()
{
    dip::Image lena_256, lena_256_blur, lena_256_blur_noise, img;
    lena_256.openRAW("lena_256.raw", dip::Size(256, 256));
    lena_256_blur.openRAW("lena_256_blur.raw", dip::Size(256, 256));
    lena_256_blur_noise.openRAW("lena_256_blur_noise.raw", dip::Size(256, 256));

    Complex2D ULMF = dft::makeUniformLinearMotionFilter(1.0, -0.07, 0.07, 256, 256);
    Complex2D Laplacian = dft::makeLaplacianFilter(4, 256, 256);

    Complex2D lena_FFT = vector2Complex2D(lena_256_blur.img);
    lena_FFT = dft::fft2D(lena_FFT, true, false);
    for (double gamma = 0.0001; gamma <= 0.1; gamma *= 10.0)
    {
        Complex2D lena_256_invF = dft::ConstrainedLeastSquareFilter(lena_FFT, ULMF, Laplacian, gamma);
        lena_256_invF = dft::fft2D(lena_256_invF, true, true);
        Byte2D result = minMaxNormalization<Byte>(lena_256_invF, 0, 255);

        std::stringstream fileName;
        fileName << "HW_2_3_" << gamma << ".raw";
        img.setImg(result);
        img.saveRAW(fileName.str().c_str(), true);

        double mse = dip::mse(lena_256, img);
        double psnr = dip::psnr(lena_256, img);
        std::cout << "MSE: " << mse << std::endl;
        std::cout << "PSNR: " << psnr << std::endl;
    }

    lena_FFT = vector2Complex2D(lena_256_blur_noise.img);
    lena_FFT = dft::fft2D(lena_FFT, true, false);
    for (double gamma = 10; gamma <= 10000; gamma *= 10)
    {
        Complex2D lena_256_invF = dft::ConstrainedLeastSquareFilter(lena_FFT, ULMF, Laplacian, gamma);
        lena_256_invF = dft::fft2D(lena_256_invF, true, true);
        Byte2D result = minMaxNormalization<Byte>(lena_256_invF, 0, 255);

        std::stringstream fileName;
        fileName << "HW_2_3_noise_" << gamma << ".raw";
        img.setImg(result);
        img.saveRAW(fileName.str().c_str(), true);

        double mse = dip::mse(lena_256, img);
        double psnr = dip::psnr(lena_256, img);
        std::cout << "MSE: " << mse << std::endl;
        std::cout << "PSNR: " << psnr << std::endl;
    }
}
