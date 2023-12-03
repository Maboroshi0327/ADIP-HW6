#include "ImageProcessing.hpp"

namespace dip
{
    Size::Size()
    {
        this->rows = 0;
        this->columns = 0;
    }

    Size::Size(int rows, int columns)
    {
        this->rows = rows;
        this->columns = columns;
    }

    Size::Size(const Size& source)
    {
        this->rows = source.rows;
        this->columns = source.columns;
    }

    int Size::area() const
    {
        return this->rows * this->columns;
    }

    Size Size::operator=(const Size& other)
    {
        this->rows = other.rows;
        this->columns = other.columns;
        return *this;
    }

    Size Size::operator+(const int& other) const
    {
        return Size(this->rows + other, this->columns + other);
    }

    Size Size::operator-(const int& other) const
    {
        return Size(this->rows - other, this->columns - other);
    }

    Size Size::operator*(const double& other) const
    {
        return Size(int(this->rows * other), int(this->columns * other));
    }

    Size Size::operator/(const double& other) const
    {
        return Size(int(this->rows / other), int(this->columns / other));
    }

    bool Size::operator==(const Size& other) const
    {
        bool result = this->rows == other.rows;
        result &= this->columns == other.columns;
        return result;
    }

    bool Size::operator!=(const Size& other) const
    {
        bool result = this->rows != other.rows;
        result |= this->columns != other.columns;
        return result;
    }

    void Image::openRAW(const char* FileName, Size size)
    {
        uint8_t* buffer = new uint8_t[size.area()];

        FILE* input_file;
        if (fopen_s(&input_file, FileName, "rb") != 0)
        {
            printf("Couldn't find file %s.\n", FileName);
            exit(1);
        }
        fread_s(buffer, size.area(), sizeof(uint8_t), size.area(), input_file);

        this->size = size;
        this->img = std::vector<std::vector<uint8_t>>(size.rows, std::vector<uint8_t>(size.columns, 0));
        for (int r = 0; r < size.rows; r++)
            for (int c = 0; c < size.columns; c++)
                this->img[r][c] = buffer[c + r * size.columns];

        delete[] buffer;
        fclose(input_file);
    }

    void Image::saveRAW(const char* FileName, bool printName)
    {
        uint8_t* buffer = new uint8_t[this->size.area()];

        for (int r = 0; r < this->size.rows; r++)
            for (int c = 0; c < this->size.columns; c++)
                buffer[c + r * this->size.columns] = this->img[r][c];

        FILE* output_file;
        fopen_s(&output_file, FileName, "wb");
        fwrite(buffer, sizeof(uint8_t), this->size.area(), output_file);

        if (printName == true)
            std::cout << "Output file: " << FileName << std::endl;

        delete[] buffer;
        fclose(output_file);
    }

    void Image::setImg(const Image& source)
    {
        this->size = source.size;
        this->img = source.img;
    }

    void Image::setImg(Size size, uint8_t value)
    {
        this->size = size;
        this->img = std::vector<std::vector<uint8_t>>(size.rows, std::vector<uint8_t>(size.columns, value));
    }

    void Image::setImg(std::vector<std::vector<uint8_t>>& img)
    {
        this->size.rows = (int)img.size();
        this->size.columns = (int)img[0].size();
        this->img = img;
    }

    Image::Image()
    {
        this->size = Size(0, 0);
    }

    Image::Image(const Image& source)
    {
        this->size = source.size;
        this->img = source.img;
    }

    Image::Image(Size size, uint8_t value)
    {
        this->size = size;
        this->img = std::vector<std::vector<uint8_t>>(size.rows, std::vector<uint8_t>(size.columns, value));
    }

    Image::Image(std::vector<std::vector<uint8_t>>& img)
    {
        this->size.rows = (int)img.size();
        this->size.columns = (int)img[0].size();
        this->img = img;
    }

    Image Image::operator=(const Image& other)
    {
        this->size = other.size;
        this->img = other.img;
        return *this;
    }

    std::vector<uint8_t>& Image::operator[](int index)
    {
        return this->img[index];
    }

    std::vector<uint8_t> Image::operator[](int index) const
    {
        return this->img[index];
    }

    double mse(const Image& img1, const Image& img2)
    {
        if (img1.size != img2.size)
        {
            std::cout << "mse: img1.size != img2.size" << std::endl;
        }

        int sum = 0;
        for (int r = 0; r < img1.size.rows; r++)
            for (int c = 0; c < img1.size.columns; c++)
                sum = sum + (img1.img[r][c] - img2.img[r][c]) * (img1.img[r][c] - img2.img[r][c]);

        double value = sum / (double)img1.size.area();

        return value;
    }

    double psnr(const Image& img1, const Image& img2)
    {
        double mse_v = mse(img1, img2);
        double param = 255 * 255 / mse_v;
        double value = 10 * std::log10(param);

        return value;
    }

    void nearest(const Image& imgIn, Image& imgOut, Size size)
    {
        Size sizeOld(imgIn.size);
        Image imgTemp(size, 0);

        for (int r = 0; r < size.rows; r++)
        {
            for (int c = 0; c < size.columns; c++)
            {
                double r_mapping = (double)r / (size.rows - 1) * (sizeOld.rows - 1);
                double c_mapping = (double)c / (size.columns - 1) * (sizeOld.columns - 1);
                imgTemp.img[r][c] = imgIn.img[(int)round(r_mapping)][(int)round(c_mapping)];
            }
        }

        imgOut.setImg(imgTemp);
    }

    void bicubic(const Image& imgIn, Image& imgOut, Size size)
    {
        int pad = 2;
        Image imgPad;
        padding(imgIn, imgPad, pad);

        Size sizeOld(size);
        Image imgTemp(size, 0);

        for (int r = 0; r < size.rows; r++)
        {
            for (int c = 0; c < size.columns; c++)
            {
                double r_mapping = (double)r / (size.rows - 1) * (sizeOld.rows - 1) + pad;
                double c_mapping = (double)c / (size.columns - 1) * (sizeOld.columns - 1) + pad;
                double r1 = floor(r_mapping);
                double c1 = floor(c_mapping);

                double sum = 0;
                for (int i = 0; i < 4; i++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                        sum += imgPad.img[(int)r1 - 1 + i][(int)c1 - 1 + j]
                            * bicubic_param(r_mapping - (r1 - 1 + i), -0.5)
                            * bicubic_param(c_mapping - (c1 - 1 + j), -0.5);
                    }
                }

                imgTemp.img[r][c] = (uint8_t)round(sum);
            }
        }

        imgOut.setImg(imgTemp);
    }

    double bicubic_param(double x, double a)
    {
        x = abs(x);
        if (x <= 1)
            return (a + 2) * x * x * x - (a + 3) * x * x + 1;
        else if (1 < x && x < 2)
            return a * x * x * x - 5 * a * x * x + 8 * a * x - 4 * a;
        else
            return 0;
    }

    void padding(const Image& imgIn, Image& imgOut, int pad)
    {
        int rows = imgIn.size.rows;
        int columns = imgIn.size.columns;
        Image imgPad(imgIn.size + 2 * pad, 0);

        for (int r = 0; r < rows; r++)
            for (int c = 0; c < columns; c++)
                imgPad.img[r + pad][c + pad] = imgIn.img[r][c];

        for (int r = 0; r < rows + pad * 2; r++)
        {
            for (int c = 0; c < columns + pad * 2; c++)
            {
                if (r < pad && c < pad)
                    imgPad.img[r][c] = imgPad.img[2 * pad - 1 - r][2 * pad - 1 - c];
                else if (r < pad && pad <= c && c < columns + pad)
                    imgPad.img[r][c] = imgPad.img[2 * pad - 1 - r][c];
                else if (r < pad && c >= columns + pad)
                    imgPad.img[r][c] = imgPad.img[2 * pad - 1 - r][2 * (columns + pad) - 1 - c];
                else if (pad <= r && r < rows + pad && c < pad)
                    imgPad.img[r][c] = imgPad.img[r][2 * pad - 1 - c];
                else if (pad <= r && r < rows + pad && c >= columns + pad)
                    imgPad.img[r][c] = imgPad.img[r][2 * (columns + pad) - 1 - c];
                else if (r >= rows + pad && c < pad)
                    imgPad.img[r][c] = imgPad.img[2 * (rows + pad) - 1 - r][2 * pad - 1 - c];
                else if (r >= rows + pad && pad <= c && c < columns + pad)
                    imgPad.img[r][c] = imgPad.img[2 * (rows + pad) - 1 - r][c];
                else if (r >= rows + pad && c >= columns + pad)
                    imgPad.img[r][c] = imgPad.img[2 * (rows + pad) - 1 - r][2 * (columns + pad) - 1 - c];
            }
        }

        imgOut.setImg(imgPad);
    }

    void padding(const Image& imgIn, Image& imgOut, int rowPad, int colPad)
    {
        int rows = imgIn.size.rows;
        int columns = imgIn.size.columns;
        Image imgPad(Size(rows + 2 * rowPad, columns + 2 * colPad), 0);

        for (int r = 0; r < rows; r++)
            for (int c = 0; c < columns; c++)
                imgPad.img[r + rowPad][c + colPad] = imgIn.img[r][c];

        for (int r = 0; r < rows + rowPad * 2; r++)
        {
            for (int c = 0; c < columns + colPad * 2; c++)
            {
                if (r < rowPad && c < colPad)
                    imgPad.img[r][c] = imgPad.img[2 * rowPad - 1 - r][2 * colPad - 1 - c];
                else if (r < rowPad && colPad <= c && c < columns + colPad)
                    imgPad.img[r][c] = imgPad.img[2 * rowPad - 1 - r][c];
                else if (r < rowPad && c >= columns + colPad)
                    imgPad.img[r][c] = imgPad.img[2 * rowPad - 1 - r][2 * (columns + colPad) - 1 - c];
                else if (rowPad <= r && r < rows + rowPad && c < colPad)
                    imgPad.img[r][c] = imgPad.img[r][2 * colPad - 1 - c];
                else if (rowPad <= r && r < rows + rowPad && c >= columns + colPad)
                    imgPad.img[r][c] = imgPad.img[r][2 * (columns + colPad) - 1 - c];
                else if (r >= rows + rowPad && c < colPad)
                    imgPad.img[r][c] = imgPad.img[2 * (rows + rowPad) - 1 - r][2 * colPad - 1 - c];
                else if (r >= rows + rowPad && colPad <= c && c < columns + colPad)
                    imgPad.img[r][c] = imgPad.img[2 * (rows + rowPad) - 1 - r][c];
                else if (r >= rows + rowPad && c >= columns + colPad)
                    imgPad.img[r][c] = imgPad.img[2 * (rows + rowPad) - 1 - r][2 * (columns + colPad) - 1 - c];
            }
        }

        imgOut.setImg(imgPad);
    }

    void threshold(const Image& imgIn, Image& imgOut, uint8_t th)
    {
        Image imgTemp(imgIn);

        for (auto& i : imgTemp.img)
            for (auto& j : i)
                j = (j < th) ? 0 : 255;

        imgOut.setImg(imgTemp);
    }

    void replaceBitPlane(const Image& imgIn, Image& imgOut, const Image& plane, int bit)
    {
        Image imgTemp(imgIn.size, 0);

        if (imgIn.size == plane.size)
        {
            for (int r = 0; r < imgIn.size.rows; r++)
                for (int c = 0; c < imgIn.size.columns; c++)
                {
                    bool binary = (plane.img[r][c] == 255) ? 1 : 0;
                    BitPlane bitplane;
                    bitplane.pixel = imgIn.img[r][c];
                    bitplane.b0 = (bit == 0) ? binary : bitplane.b0;
                    bitplane.b1 = (bit == 1) ? binary : bitplane.b1;
                    bitplane.b2 = (bit == 2) ? binary : bitplane.b2;
                    bitplane.b3 = (bit == 3) ? binary : bitplane.b3;
                    bitplane.b4 = (bit == 4) ? binary : bitplane.b4;
                    bitplane.b5 = (bit == 5) ? binary : bitplane.b5;
                    bitplane.b6 = (bit == 6) ? binary : bitplane.b6;
                    bitplane.b7 = (bit == 7) ? binary : bitplane.b7;
                    imgTemp.img[r][c] = bitplane.pixel;
                }
        }
        else
        {
            std::cout << "Image::size != img_in.size" << std::endl;
            return;
        }

        imgOut.setImg(imgTemp);
    }

    void extractBitPlane(const Image& imgIn, Image& imgOut, int bit)
    {
        Image imgTemp(imgIn.size, 0);

        for (int r = 0; r < imgIn.size.rows; r++)
            for (int c = 0; c < imgIn.size.columns; c++)
            {
                BitPlane bitplane;
                bitplane.pixel = imgIn.img[r][c];
                imgTemp.img[r][c] = (bit == 0) ? (bitplane.b0 * 255) : imgTemp.img[r][c];
                imgTemp.img[r][c] = (bit == 1) ? (bitplane.b1 * 255) : imgTemp.img[r][c];
                imgTemp.img[r][c] = (bit == 2) ? (bitplane.b2 * 255) : imgTemp.img[r][c];
                imgTemp.img[r][c] = (bit == 3) ? (bitplane.b3 * 255) : imgTemp.img[r][c];
                imgTemp.img[r][c] = (bit == 4) ? (bitplane.b4 * 255) : imgTemp.img[r][c];
                imgTemp.img[r][c] = (bit == 5) ? (bitplane.b5 * 255) : imgTemp.img[r][c];
                imgTemp.img[r][c] = (bit == 6) ? (bitplane.b6 * 255) : imgTemp.img[r][c];
                imgTemp.img[r][c] = (bit == 7) ? (bitplane.b7 * 255) : imgTemp.img[r][c];
            }

        imgOut.setImg(imgTemp);
    }

    void logTrans(const Image& imgIn, Image& imgOut)
    {
        Image imgTemp(imgIn);

        const double c = 255 / log(1 + 255);
        for (auto& i : imgTemp.img)
            for (auto& j : i)
            {
                j = (uint8_t)round(c * std::log(1 + j));
            }

        imgOut.setImg(imgTemp);
    }

    void invLogTrans(const Image& imgIn, Image& imgOut)
    {
        Image imgTemp(imgIn);

        const double c = 255 / log(1 + 255);
        for (auto& i : imgTemp.img)
            for (auto& j : i)
            {
                j = (uint8_t)std::round(std::exp(j / c) - 1);
            }

        imgOut.setImg(imgTemp);
    }

    void gammaTrans(const Image& imgIn, Image& imgOut, double gamma)
    {
        Image imgTemp(imgIn);

        const double c = 255 / pow(255, gamma);
        for (auto& i : imgTemp.img)
            for (auto& j : i)
            {
                j = (uint8_t)std::round(c * std::pow(j, gamma));
            }

        imgOut.setImg(imgTemp);
    }

    void negative(const Image& imgIn, Image& imgOut)
    {
        Image imgTemp(imgIn);

        for (auto& i : imgTemp.img)
            for (auto& j : i)
            {
                j = 255 - j;
            }

        imgOut.setImg(imgTemp);
    }

    void histEqualize(const Image& imgIn, Image& imgOut)
    {
        Image imgTemp(imgIn);

        std::vector<int> hist(256, 0);
        for (auto& i : imgTemp.img)
            for (auto& j : i)
                hist[j] += 1;
        double area = (double)imgTemp.size.area();

        std::vector<uint8_t> cdf(256, 0);
        double sum = 0;
        for (int i = 0; i < 256; i++)
        {
            sum += hist[i] / area;
            cdf[i] = (uint8_t)round(sum * 255);
        }

        for (auto& i : imgTemp.img)
            for (auto& j : i)
                j = cdf[j];

        imgOut.setImg(imgTemp);
    }

    int moments(const Image& imgIn, int i, int j)
    {
        int M = 0;
        for (int x = 0; x < imgIn.size.columns; x++)
            for (int y = 0; y < imgIn.size.rows; y++)
                M += intPow(x, i) * intPow(y, j) * imgIn.img[y][x];
        return M;
    }

    int intPow(int x, int p)
    {
        int result = 1;
        for (int i = 0; i < p; i++)
            result *= x;

        return result;
    }

    Coordinate<double> centroid(const Image& imgIn)
    {
        int M00, M01, M10;
        M00 = moments(imgIn, 0, 0);
        M01 = moments(imgIn, 0, 1);
        M10 = moments(imgIn, 1, 0);
        return Coordinate<double>((double)M01 / M00, (double)M10 / M00);
    }

    double centralMoments(const Image& imgIn, int p, int q)
    {
        Coordinate<double> cen = centroid(imgIn);
        double result = 0.0;
        for (int x = 0; x < imgIn.size.columns; x++)
            for (int y = 0; y < imgIn.size.rows; y++)
                result += pow(x - cen.column, p) * pow(y - cen.row, q) * imgIn.img[y][x];
        return result;
    }

    void boxedBlur(const Image& imgIn, Image& imgOut, Size kernelSize)
    {
        Kernel<int> kernel = makeBoxKernel(kernelSize);
        imgOut.img = convolution<uint8_t, int>(imgIn, kernel, true);
        imgOut.size = imgIn.size;
    }

    void gaussianBlur(const Image& imgIn, Image& imgOut, Size kernelSize, double sigma)
    {
        Kernel<int> kernel = makeGaussianKernel(kernelSize, sigma);
        imgOut.img = convolution<uint8_t, int>(imgIn, kernel, true);
        imgOut.size = imgIn.size;
    }

    void robertsFiltering(const Image& imgIn, Image& imgOut)
    {
        Kernel<int> kernel;
        std::vector<std::vector<int>> vectorN45, vector45;
        kernel = makeRobertsKernel(-45);
        vectorN45 = convolution<int, int>(imgIn, kernel, false);
        kernel = makeRobertsKernel(45);
        vector45 = convolution<int, int>(imgIn, kernel, false);

        Image imgTemp(imgIn.size, 0);
        for (int r = 0; r < imgIn.size.rows; r++)
            for (int c = 0; c < imgIn.size.columns; c++)
            {
                //int abs = std::abs(vectorN45[r][c]) + std::abs(vector45[r][c]);
                int abs = (int)std::sqrt(vectorN45[r][c] * vectorN45[r][c] + vector45[r][c] + vector45[r][c]);
                imgTemp.img[r][c] = (abs > 255) ? 255 : abs;
            }

        imgOut.setImg(imgTemp);
    }

    void robertsFiltering(const Image& imgIn, Image& imgOut, int angle)
    {
        if (angle == -45 || angle == 45)
        {
            Kernel<int> kernel = makeRobertsKernel(angle);
            std::vector<std::vector<int>> conv = convolution<int, int>(imgIn, kernel, false);

            Image imgTemp(imgIn.size, 0);
            for (int r = 0; r < imgIn.size.rows; r++)
                for (int c = 0; c < imgIn.size.columns; c++)
                {
                    int abs = std::abs(conv[r][c]);
                    imgTemp.img[r][c] = (abs > 255) ? 255 : abs;
                }

            imgOut.setImg(imgTemp);
        }
        else
        {
            std::cout << "robertsFiltering: angle != -45 or angle != 45" << std::endl;
        }
    }

    void prewittFiltering(const Image& imgIn, Image& imgOut)
    {
        Kernel<int> kernel;
        std::vector<std::vector<int>> vector0, vector90;

        Image imgTemp(imgIn.size, 0);
        kernel = makePrewittKernel(0);
        vector0 = convolution<int, int>(imgIn, kernel, false);
        kernel = makePrewittKernel(90);
        vector90 = convolution<int, int>(imgIn, kernel, false);

        for (int r = 0; r < imgIn.size.rows; r++)
            for (int c = 0; c < imgIn.size.columns; c++)
            {
                //int abs = std::abs(vector0[r][c]) + std::abs(vector90[r][c]);
                int abs = (int)std::sqrt(vector0[r][c] * vector0[r][c] + vector90[r][c] * vector90[r][c]);
                imgTemp.img[r][c] = (abs > 255) ? 255 : abs;
            }

        imgOut.setImg(imgTemp);
    }

    void prewittFiltering(const Image& imgIn, Image& imgOut, int angle)
    {
        if (angle == 0 || angle == 90)
        {
            Kernel<int> kernel = makePrewittKernel(angle);
            std::vector<std::vector<int>> conv = convolution<int, int>(imgIn, kernel, false);

            Image imgTemp(imgIn.size, 0);
            for (int r = 0; r < imgIn.size.rows; r++)
                for (int c = 0; c < imgIn.size.columns; c++)
                {
                    int abs = std::abs(conv[r][c]);
                    imgTemp.img[r][c] = (abs > 255) ? 255 : abs;
                }

            imgOut.setImg(imgTemp);
        }
        else
        {
            std::cout << "prewittFiltering: angle != 0 or angle != 90" << std::endl;
        }
    }

    void sobelFiltering(const Image& imgIn, Image& imgOut)
    {
        Kernel<int> kernel;
        std::vector<std::vector<int>> vector0, vector90;

        Image imgTemp(imgIn.size, 0);
        kernel = makeSobelKernel(0);
        vector0 = convolution<int, int>(imgIn, kernel, false);
        kernel = makeSobelKernel(90);
        vector90 = convolution<int, int>(imgIn, kernel, false);

        for (int r = 0; r < imgIn.size.rows; r++)
            for (int c = 0; c < imgIn.size.columns; c++)
            {
                //int abs = std::abs(vector0[r][c]) + std::abs(vector90[r][c]);
                int abs = (int)std::sqrt(vector0[r][c] * vector0[r][c] + vector90[r][c] * vector90[r][c]);
                imgTemp.img[r][c] = (abs > 255) ? 255 : abs;
            }

        imgOut.setImg(imgTemp);
    }

    void sobelFiltering(const Image& imgIn, Image& imgOut, int angle)
    {
        if (angle == -45 || angle == 0 || angle == 45 || angle == 90)
        {
            Kernel<int> kernel = makeSobelKernel(angle);
            std::vector<std::vector<int>> conv = convolution<int, int>(imgIn, kernel, false);

            Image imgTemp(imgIn.size, 0);
            for (int r = 0; r < imgIn.size.rows; r++)
                for (int c = 0; c < imgIn.size.columns; c++)
                {
                    int abs = std::abs(conv[r][c]);
                    imgTemp.img[r][c] = (abs > 255) ? 255 : abs;
                }

            imgOut.setImg(imgTemp);
        }
        else
        {
            std::cout << "sobelFiltering: angle != -45 or angle != 0 or angle != 45 or angle != 90" << std::endl;
        }
    }

    void laplacianFiltering(const Image& imgIn, Image& imgOut, int neighbor)
    {
        if (neighbor == 4 || neighbor == 8)
        {
            Kernel<int> kernel = makeLaplacianKernel(neighbor);
            std::vector<std::vector<int>> conv = convolution<int, int>(imgIn, kernel, false);

            Image imgTemp(imgIn.size, 0);
            for (int r = 0; r < imgIn.size.rows; r++)
                for (int c = 0; c < imgIn.size.columns; c++)
                {
                    int abs = std::abs(conv[r][c]);
                    imgTemp.img[r][c] = (abs > 255) ? 255 : abs;
                }

            imgOut.setImg(imgTemp);
        }
        else
        {
            std::cout << "sobelFiltering: neighbor != 4 or neighbor != 8" << std::endl;
        }
    }

    void laplacianSharpening(const Image& imgIn, Image& imgOut, int neighbor)
    {
        if (neighbor == 4 || neighbor == 8)
        {
            Kernel<int> kernel = makeLaplacianKernel(neighbor);
            kernel[1][1] += 1;
            std::vector<std::vector<int>> conv = convolution<int, int>(imgIn, kernel, false);

            Image imgTemp(imgIn.size, 0);
            for (int r = 0; r < imgIn.size.rows; r++)
                for (int c = 0; c < imgIn.size.columns; c++)
                {
                    int abs = std::abs(conv[r][c]);
                    imgTemp.img[r][c] = (abs > 255) ? 255 : abs;
                }

            imgOut.setImg(imgTemp);
        }
        else
        {
            std::cout << "sobelFiltering: neighbor != 4 or neighbor != 8" << std::endl;
        }
    }

    Kernel<int> makeBoxKernel(Size size)
    {
        Kernel<int> kernel(size, 1);
        return kernel;
    }

    Kernel<int> makeGaussianKernel(Size size, double sigma)
    {
        Kernel<int> kernel(size);

        int rowsKernel = kernel.size.rows;
        int colsKernel = kernel.size.columns;
        int rowOrigin = kernel.origin.row;
        int colOrigin = kernel.origin.column;

        double min = 1;
        //double denominator = 2 * M_PI * sigma * sigma;
        std::vector<std::vector<double>> kernel_temp(rowsKernel, std::vector<double>(colsKernel, 0.0));
        for (int i = -rowOrigin; i < rowsKernel - rowOrigin; i++)
            for (int j = -colOrigin; j < colsKernel - colOrigin; j++)
            {
                double power = -(i * i + j * j) / (2 * sigma * sigma);
                double exp_v = std::exp(power);
                kernel_temp[i + rowOrigin][j + colOrigin] = exp_v;
                if (min > exp_v)
                    min = exp_v;
            }

        for (auto& i : kernel_temp)
            for (auto& j : i)
                j = j / min;

        for (int i = 0; i < rowsKernel; i++)
            for (int j = 0; j < colsKernel; j++)
                kernel[i][j] = (int)std::round(kernel_temp[i][j]);

        return kernel;
    }

    Kernel<double> makeGaussianKernel(Size size, double sigma, bool normalize)
    {
        Kernel<double> kernel(size);
        int rowsKernel = kernel.size.rows;
        int colsKernel = kernel.size.columns;
        int rowOrigin = kernel.origin.row;
        int colOrigin = kernel.origin.column;

        double sum = 0, min = 1;
        //double denominator = 2 * M_PI * sigma * sigma;
        for (int i = -rowOrigin; i < rowsKernel - rowOrigin; i++)
            for (int j = -colOrigin; j < colsKernel - colOrigin; j++)
            {
                double power = -(i * i + j * j) / (2 * sigma * sigma);
                double exp_v = std::exp(power);
                kernel[i + rowOrigin][j + colOrigin] = exp_v;
                sum += exp_v;
                if (min > exp_v)
                    min = exp_v;
            }

        if (normalize)
        {
            for (auto& i : kernel.kernel)
                for (auto& j : i)
                    j = j / sum;
        }
        else
        {
            for (auto& i : kernel.kernel)
                for (auto& j : i)
                    j = j / min;
        }

        return kernel;
    }

    Kernel<int> makeRobertsKernel(int angle)
    {
        Kernel<int> kernel(Size(2, 2), Coordinate<int>(0, 0));

        if (angle == -45)
        {
            kernel[0][0] = -1; kernel[0][1] = 0;
            kernel[1][0] = 0;  kernel[1][1] = 1;
        }

        if (angle == 45)
        {
            kernel[0][0] = 0; kernel[0][1] = -1;
            kernel[1][0] = 1; kernel[1][1] = 0;
        }

        return kernel;
    }

    Kernel<int> makePrewittKernel(int angle)
    {
        Kernel<int> kernel(Size(3, 3));

        if (angle == 0)
        {
            kernel[0][0] = -1; kernel[0][1] = 0; kernel[0][2] = 1;
            kernel[1][0] = -1; kernel[1][1] = 0; kernel[1][2] = 1;
            kernel[2][0] = -1; kernel[2][1] = 0; kernel[2][2] = 1;
        }

        if (angle == 90)
        {
            kernel[0][0] = 1;  kernel[0][1] = 1;  kernel[0][2] = 1;
            kernel[1][0] = 0;  kernel[1][1] = 0;  kernel[1][2] = 0;
            kernel[2][0] = -1; kernel[2][1] = -1; kernel[2][2] = -1;
        }

        return kernel;
    }

    Kernel<int> makeSobelKernel(int angle)
    {
        Kernel<int> kernel(Size(3, 3));

        if (angle == -45)
        {
            kernel[0][0] = -2; kernel[0][1] = -1; kernel[0][2] = 0;
            kernel[1][0] = -1; kernel[1][1] = 0;  kernel[1][2] = 1;
            kernel[2][0] = 0;  kernel[2][1] = 1;  kernel[2][2] = 2;
        }

        if (angle == 0)
        {
            kernel[0][0] = -1; kernel[0][1] = 0; kernel[0][2] = 1;
            kernel[1][0] = -2; kernel[1][1] = 0; kernel[1][2] = 2;
            kernel[2][0] = -1; kernel[2][1] = 0; kernel[2][2] = 1;
        }

        if (angle == 45)
        {
            kernel[0][0] = 0;  kernel[0][1] = 1;  kernel[0][2] = 2;
            kernel[1][0] = -1; kernel[1][1] = 0;  kernel[1][2] = 1;
            kernel[2][0] = -2; kernel[2][1] = -1; kernel[2][2] = 0;
        }

        if (angle == 90)
        {
            kernel[0][0] = 1;  kernel[0][1] = 2;  kernel[0][2] = 1;
            kernel[1][0] = 0;  kernel[1][1] = 0;  kernel[1][2] = 0;
            kernel[2][0] = -1; kernel[2][1] = -2; kernel[2][2] = -1;
        }

        return kernel;
    }

    Kernel<int> makeLaplacianKernel(int neighbor)
    {
        Kernel<int> kernel(Size(3, 3));

        if (neighbor == 4)
        {
            kernel[0][0] = 0;  kernel[0][1] = -1; kernel[0][2] = 0;
            kernel[1][0] = -1; kernel[1][1] = 4;  kernel[1][2] = -1;
            kernel[2][0] = 0;  kernel[2][1] = -1; kernel[2][2] = 0;
        }

        if (neighbor == 8)
        {
            kernel[0][0] = -1; kernel[0][1] = -1; kernel[0][2] = -1;
            kernel[1][0] = -1; kernel[1][1] = 8;  kernel[1][2] = -1;
            kernel[2][0] = -1; kernel[2][1] = -1; kernel[2][2] = -1;
        }

        return kernel;
    }

    void AdaptiveLocalNoiseReductionFiltering(const Image& imgIn, Image& imgOut, Size kernelSize)
    {
        int rows = imgIn.size.rows;
        int columns = imgIn.size.columns;
        int rowsKernel = kernelSize.rows;
        int colsKernel = kernelSize.columns;
        int rowOrigin = rowsKernel / 2;
        int colOrigin = colsKernel / 2;
        double invArea = 1.0 / kernelSize.area();

        Image imgPad;
        int rowPad = rowsKernel - rowOrigin - 1;
        int colPad = colsKernel - colOrigin - 1;
        padding(imgIn, imgPad, rowPad, colPad);

        double noiseVariance = 0.0;
        std::vector<std::vector<double>> localMean(rows, std::vector<double>(columns));
        std::vector<std::vector<double>> localVariance(rows, std::vector<double>(columns));

        Image imgTemp(imgIn.size, 0);
        for (int r = rowPad; r < rows + rowPad; r++)
            for (int c = colPad; c < columns + colPad; c++)
            {
                double localMean_v = 0.0;
                for (int i = -rowOrigin; i < rowsKernel - rowOrigin; i++)
                    for (int j = -colOrigin; j < colsKernel - colOrigin; j++)
                        localMean_v += imgPad[r + i][c + j];
                localMean_v *= invArea;
                localMean[r - rowPad][c - colPad] = localMean_v;

                double localVariance_v = 0.0;
                for (int i = -rowOrigin; i < rowsKernel - rowOrigin; i++)
                    for (int j = -colOrigin; j < colsKernel - colOrigin; j++)
                        localVariance_v += std::pow(imgPad[r + i][c + j] - localMean_v, 2);
                localVariance_v *= invArea;
                localVariance[r - rowPad][c - colPad] = localVariance_v;

                noiseVariance += localVariance_v;
            }
        noiseVariance /= rows * columns;

        for (int r = 0; r < rows; r++)
            for (int c = 0; c < columns; c++)
            {
                double ALNRF;
                double localVariance_v = localVariance[r][c];
                double localMean_v = localMean[r][c];

                if (noiseVariance > localVariance_v)
                    ALNRF = localMean_v;
                else
                    ALNRF = imgIn[r][c] - (noiseVariance / localVariance_v) * (imgIn[r][c] - localMean_v);

                imgTemp[r][c] = (uint8_t)std::round(ALNRF);
            }

        imgOut.setImg(imgTemp);
    }

    void AlphaTrimmedMeanFiltering(const Image& imgIn, Image& imgOut, Size kernelSize, int d)
    {
        int rows = imgIn.size.rows;
        int columns = imgIn.size.columns;
        int rowsKernel = kernelSize.rows;
        int colsKernel = kernelSize.columns;
        int rowOrigin = rowsKernel / 2;
        int colOrigin = colsKernel / 2;
        int areaKernel = kernelSize.area();

        if (areaKernel <= d)
        {
            std::cout << "AlphaTrimmedMeanFiltering: Kernel <= d" << std::endl;
            d = areaKernel - 1;
        }

        Image imgPad;
        int rowPad = rowsKernel - rowOrigin - 1;
        int colPad = colsKernel - colOrigin - 1;
        padding(imgIn, imgPad, rowPad, colPad);

        Image imgTemp(imgIn.size, 0);
        std::vector<uint8_t> local(kernelSize.area());
        for (int r = rowPad; r < rows + rowPad; r++)
            for (int c = colPad; c < columns + colPad; c++)
            {
                int index = 0;
                for (int i = -rowOrigin; i < rowsKernel - rowOrigin; i++)
                    for (int j = -colOrigin; j < colsKernel - colOrigin; j++)
                    {
                        local[index] = imgPad[r + i][c + j];
                        index++;
                    }
                std::sort(local.begin(), local.end());

                double mean = 0.0;
                for (index = 0; index < local.size() - d; index++)
                    mean += local[index + d / 2];
                mean /= index;

                imgTemp[r - rowPad][c - colPad] = (uint8_t)std::round(mean);
            }

        imgOut.setImg(imgTemp);
    }

    void PerspectiveTransformation(const Image& imgIn, Image& imgOut, Coordinate<int> xy[4], Coordinate<int> uv[4])
    {
        // Case 1
        double a, b, c, d, e, f, g, h;
        parameterPerspective(a, b, c, d, e, f, g, h, uv);

        // Case 2
        double aa, bb, cc, dd, ee, ff, gg, hh;
        parameterPerspective(aa, bb, cc, dd, ee, ff, gg, hh, xy);

        int newRows = std::max(uv[1].row, uv[2].row) + 1;
        int newCols = std::max(uv[2].column, uv[3].column) + 1;

        Image imageOut(Size(newRows, newCols), 0);
        for (int u = 0; u < newRows; u++)
            for (int v = 0; v < newCols; v++)
            {
                // Case 1
                std::vector<std::vector<double>> invMatrix(3, std::vector<double>(3));
                invMatrix[0][0] = a - h * u;		invMatrix[0][1] = b - h * u;		invMatrix[0][2] = c;
                invMatrix[1][0] = d - g * v;		invMatrix[1][1] = e - h * v;		invMatrix[1][2] = f;
                invMatrix[2][0] = 0.0;				invMatrix[2][1] = 0.0;				invMatrix[2][2] = 1.0;
                invMatrix = inverseMatrix(invMatrix);

                std::vector<std::vector<double>> uvMatrix(3, std::vector<double>(1));
                uvMatrix[0][0] = u;
                uvMatrix[1][0] = v;
                uvMatrix[2][0] = 1.0;

                std::vector<std::vector<double>> uuvvMatrix = multipleMatrix(invMatrix, uvMatrix);
                double uu = uuvvMatrix[0][0];
                double vv = uuvvMatrix[1][0];
                if (uu > 1 || uu < 0 || vv > 1 || vv < 0)
                    continue;

                // Case 2
                int x = (int)std::round((aa * uu + bb * vv + cc) / (gg * uu + hh * vv + 1));
                int y = (int)std::round((dd * uu + ee * vv + ff) / (gg * uu + hh * vv + 1));

                imageOut[u][v] = imgIn[x][y];
            }

        imgOut.setImg(imageOut);
    }

    void parameterPerspective(double& a, double& b, double& c, double& d, double& e, double& f, double& g, double& h, Coordinate<int> xy[4])
    {
        int x0 = xy[0].row, y0 = xy[0].column;
        int x1 = xy[1].row, y1 = xy[1].column;
        int x2 = xy[2].row, y2 = xy[2].column;
        int x3 = xy[3].row, y3 = xy[3].column;
        int xSigma = x0 - x1 + x2 - x3;
        int ySigma = y0 - y1 + y2 - y3;
        int xDelta1 = x1 - x2;
        int xDelta2 = x3 - x2;
        int yDelta1 = y1 - y2;
        int yDelta2 = y3 - y2;
        g = (double)(xSigma * yDelta2 - ySigma * xDelta2) / (xDelta1 * yDelta2 - yDelta1 * xDelta2);
        h = (double)(xDelta1 * ySigma - yDelta1 * xSigma) / (xDelta1 * yDelta2 - yDelta1 * xDelta2);
        a = x1 - x0 + g * x1;
        b = x3 - x0 + h * x3;
        c = x0;
        d = y1 - y0 + g * y1;
        e = y3 - y0 + h * y3;
        f = y0;
    }
}
