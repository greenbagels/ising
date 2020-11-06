#include <iostream>
#include <cmath>
#include <random>
#include <stdexcept>
#include <omp.h>
// #include <png.h>
#include <png++/png.hpp>
// #include <zlib.h>

const auto n = 512lu;
const auto T = 0.1;

void initialize(int spins[][n])
{
    std::random_device rd;
    std::mt19937 engine(rd());
//    #pragma omp parallel for
    for (auto i = 0; i < n; i++)
    {
        for (auto j = 0; j < n; j++)
        {
            spins[i][j] = 1 - 2 * (engine() % 2);
        }
    }
}

char get_color(int spins[][n], int i, int j)
{
    if (spins[i][j] == 1)
    {
        return 'x';
    }
    else if (spins[i][j] == -1)
    {
        return '.';
    }
    else
    {
        throw std::runtime_error("Invalid spin state detected!\n");
    }
}

double calc_deltaU(int spins[][n], int i, int j)
{
    int top, bottom, left, right;

    switch (i)
    {
        case 0:
            top = spins[n - 1][j];
            bottom = spins[i + 1][j];
            break;
        case n-1:
            top = spins[i - 1][j];
            bottom = spins[0][j];
            break;
        default:
            top = spins[i - 1][j];
            bottom = spins[i + 1][j];
            break;
    }

    switch (j)
    {
        case 0:
            left = spins[i][n-1];
            right = spins[i][j+1];
            break;
        case n-1:
            left = spins[i][j-1];
            right = spins[i][0];
            break;
        default:
            left = spins[i][j-1];
            right = spins[i][j+1];
            break;
    }
    return 2 * spins[i][j] * (top + bottom + left + right);
}

void snapshot(const char* fname, int spins[][n])
{
    /*
    FILE *fp = fopen(fname, "wb");

    if (!fp)
        throw std::runtime_error("File creation failed!");

    png_structp png_ptr = png_create_write_struct
        (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr)
        throw std::runtime_error("PNG Pointer creation failed!");

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
    {
        png_destroy_write_struct(&png_ptr, (png_infopp)nullptr);
        throw std::runtime_error("Info pointer creation failed!");
    }

    if (setjmp(png_jmpbuf(png_ptr)))
    {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fp);
        throw std::runtime_error("setjmp failed!");
    }

    png_init_io(png_ptr, fp);

    png_set_IHDR(png_ptr, info_ptr, n, n, 1, PNG_COLOR_TYPE_GRAY,
            PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
            PNG_FILTER_TYPE_DEFAULT);

    png_write_info(png_ptr, info_ptr);
    */

    png::image<png::gray_pixel_1> img(n, n);
//    #pragma omp parallel for
    for (auto y = 0; y < n; ++y)
    {
        for (auto x = 0; x < n; x++)
        {
            img[y][x] = png::gray_pixel_1(spins[y][x] > 0);
        }
    }
    img.write(fname);
}

int main()
{
    int (*grid)[n] = new int[n][n];
    if (grid == nullptr)
    {
        throw std::runtime_error("Grid allocation failed!");
    }
    initialize(grid);
    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<int> dist(0, n - 1);
    std::uniform_real_distribution<double> floatdist;
    auto limit = 100000u*n*n;
    // We aim for 30 seconds of 60Hz footage, so 1800 frames. So let's bump
    // the limit up as little as necessary:
    auto step = limit / 20;
    std::cout << "Performing " << limit << " iterations...\n";
    for (auto t = 0; t < limit; t++)
    {
        auto i = dist(engine);
        auto j = dist(engine);
        auto dU = calc_deltaU(grid, i, j);
        if (dU <= 0.)
        {
            grid[i][j] = -grid[i][j];
        }
        else
        {
            if (floatdist(engine) < std::exp(-dU / T))
            {
                grid[i][j] = -grid[i][j];
            }
        }
        if (t % step == 0)
        {
            char filename[256];
            snprintf(filename, 256, "snapshot_%lux%lu_%.8lu.png", n, n, t / step);
            snapshot(filename, grid);
        }
    }

    /*
    for (auto i = 0; i < n; i++)
    {
        for (auto j = 0; j < n; j++)
        {
            std::cout << get_color(grid, i, j);
        }
        std::cout << std::endl;
    }
    */
}
