template<typename T>
std::vector<std::vector<T>> minMaxNormalization(const Double2D& data, double min_v, double max_v)
{
    int rows = (int)data.size();
    int cols = (int)data[0].size();

    double max = data[0][0];
    double min = max;
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
        {
            if (max < data[r][c])
                max = data[r][c];
            if (min > data[r][c])
                min = data[r][c];
        }

    std::vector<std::vector<T>> result(rows, std::vector<T>(cols, 0));
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
        {
            result[r][c] = (T)((data[r][c] - min) / (max - min) * (max_v - min_v) + min_v);
        }

    return result;
}

template<typename T>
std::vector<std::vector<T>> minMaxNormalization(const Complex2D& data, double min_v, double max_v)
{
    int rows = (int)data.size();
    int cols = (int)data[0].size();

    double max = std::abs(data[0][0]);
    double min = max;
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
        {
            double abs = std::abs(data[r][c]);
            if (max < abs)
                max = abs;
            if (min > abs)
                min = abs;
        }

    std::vector<std::vector<T>> result(rows, std::vector<T>(cols, 0));
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
        {
            double abs = std::abs(data[r][c]);
            result[r][c] = (T)((abs - min) / (max - min) * (max_v - min_v) + min_v);
        }

    return result;
}

template<typename T>
Complex2D vector2Complex2D(const std::vector<std::vector<T>>& data)
{
    int rows = (int)data.size();
    int cols = (int)data[0].size();

    Complex2D result(rows, Complex1D(cols));
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
        {
            result[r][c] = Complex(data[r][c], 0);
        }

    return result;
}

template<typename T>
Double2D vector2Double2D(const std::vector<std::vector<T>>& data)
{
    int rows = (int)data.size();
    int cols = (int)data[0].size();

    Double2D result(rows, Double1D(cols));
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
        {
            result[r][c] = (double)data[r][c];
        }

    return result;
}


template<typename T>
std::vector<std::vector<T>> complex2dReal2vector(const Complex2D& data)
{
    int rows = (int)data.size();
    int cols = (int)data[0].size();
    std::vector<std::vector<T>> result(rows, std::vector<T>(cols));
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
        {
            result[r][c] = (T)std::abs(data[r][c].real());
        }

    return result;
}

template<typename T>
std::vector<std::vector<T>> complex2dAbs2vector(const Complex2D& data)
{
    int rows = (int)data.size();
    int cols = (int)data[0].size();
    std::vector<std::vector<T>> result(rows, std::vector<T>(cols));
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
        {
            result[r][c] = (T)std::abs(data[r][c]);
        }

    return result;
}

template<typename T>
std::vector<std::vector<T>> double2d2vector(const Double2D& data)
{
    int rows = (int)data.size();
    int cols = (int)data[0].size();
    std::vector<std::vector<T>> result(rows, std::vector<T>(cols));
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
        {
            result[r][c] = (T)data[r][c];
        }

    return result;
}
