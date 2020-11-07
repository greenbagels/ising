// C++ Standard Libraries
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

// System Libraries
#include <png++/png.hpp>

// Project Libraries
#include "ising.hpp"

ising::ising(std::size_t iters, std::size_t width, unsigned neighbors,
             std::size_t nimg, unsigned scale, double temp, double fstr,
             std::string backend)
{
    this->iters = iters;
    this->width = width;
    this->neighbors = neighbors;
    this->scale = scale;

    switch(neighbors)
    {
        case 2:
            dim = 1;
            spins = std::make_unique<char[]>(width);
            break;
        case 4:
            dim = 2;
            spins = std::make_unique<char[]>(width * width);
            break;
        case 6:
        case 8:
        case 12:
            throw std::runtime_error("Not implemented!");
            break;
        default:
            throw std::runtime_error("Option parsing missed invalid neighbor quantity!");
    }

    initialize_spins();

    this->nimg = nimg;
    this->temp = temp;
    this->field_strength = fstr;
    if (backend == "cpu-serial")
    {
        this->backend = 0;
    }
    else if (backend == "cpu-parallel")
    {
        this->backend = 1;
    }
    else if (backend == "gpu")
    {
        this->backend = 2;
    }
}

double ising::calc_deltaU(unsigned i, unsigned j)
{
    int top, bottom, left, right;

    if (i == 0)
    {
        top = spins[(width - 1) * width + j];
        bottom = spins[(i + 1) * width + j];
    }
    else if (i == width - 1)
    {
        top = spins[(i - 1) * width + j];
        bottom = spins[0 * width + j];
    }
    else
    {
        top = spins[(i - 1) * width + j];
        bottom = spins[(i + 1) * width + j];
    }

    if (j == 0)
    {
        left = spins[i * width + width - 1];
        right = spins[i * width + j + 1];
    }
    else if (j == width - 1)
    {
        left = spins[i * width + j - 1];
        right = spins[i * width + 0];
    }
    else
    {
        left = spins[i * width + j - 1];
        right = spins[i * width + j + 1];
    }

    return 2 * spins[i * width + j] * (top + bottom + left + right);
}


void ising::set_display_mode(unsigned char mode)
{
    if (mode < 2)
    {
        display_mode = mode;
    }
    else
    {
        throw std::runtime_error("Invalid display mode supplied!");
    }
}

void ising::set_benchmark_mode(unsigned char mode)
{
    if (mode < 2)
    {
        benchmark_mode = mode;
    }
    else
    {
        throw std::runtime_error("Invalid benchmark mode supplied!");
    }
}

void ising::initialize_spins()
{
    std::random_device rd;
    std::mt19937 engine(rd());

    std::size_t size;

    if (dim == 1)
    {
        size = width;
    }
    else if (dim == 2)
    {
        size = width * width;
    }
    else if (dim == 3)
    {
        size = width * width * width;
        throw std::runtime_error("Dim-3 not implemented!");
    }
    else
    {
        throw std::runtime_error("Dimensionality constraint violated!");
    }

    for (auto i = 0; i < size; i++)
    {
        spins[i] = 1 - 2 * (engine() % 2);
    }
}

void ising::run()
{
    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<int> dist(0, width - 1);
    std::uniform_real_distribution<double> floatdist;
    // We aim for 30 seconds of 60Hz footage, so 1800 frames. So let's bump
    // the limit up as little as necessary:
    iters += nimg - (iters % nimg);
    auto step = iters / nimg;
    std::cerr << "Performing " << iters << " iterations...\n";
    for (auto t = 0; t < iters; t++)
    {
        auto i = dist(engine);
        auto j = dist(engine);
        auto dU = calc_deltaU(i, j);
        if (dU <= 0.)
        {
            flip_spin(i,j);
        }
        else
        {
            if (floatdist(engine) < std::exp(-dU / temp))
            {
                flip_spin(i,j);
            }
        }
        if (t % (iters / nimg) == 0)
        {
            if (display_mode == 0)
            {
                print_snapshot();
            }
            else
            {
                char filename[256];
                snprintf(filename, 256, "snapshot_%lux%lu_%.8lu.png", width, width, t / step);
                save_png_snapshot(filename);
            }
        }
    }
}

inline void ising::flip_spin(std::size_t i, std::size_t j)
{
    spins[i * width + j] = -spins[i * width + j];
}

void ising::save_png_snapshot(const char* fname)
{
    // TODO: generalize to more dims
    png::image<png::gray_pixel_1> img(scale * width, scale * width);
    for (auto y = 0; y < width; ++y)
    {
        for (auto i = 0; i < scale; i++)
        {
            for (auto x = 0; x < width; x++)
            {
                for (auto j = 0; j < scale; j++)
                {
                    img[scale*y + i][scale*x + j] = png::gray_pixel_1(spins[y * width + x] > 0);
                }
            }
        }
    }
    img.write(fname);
}

void ising::print_snapshot()
{
    for (auto y = 0; y < width; ++y)
    {
        for (auto i = 0; i < scale; i++)
        {
            for (auto x = 0; x < width; x++)
            {
                for (auto j = 0; j < scale; j++)
                {
                    std::cout << (spins[y * width + x] > 0 ? '.' : '+');
                }
            }
            std::cout << '\n';
        }
    }
    std::cout << "\n\n";
}
