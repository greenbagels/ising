// C++ Standard Libraries
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

// C Libraries
#include <stdio.h>

// System Libraries
#include <png++/png.hpp>
#include <omp.h>
#include <cuda.h>

// Project Libraries
#include "ising.hpp"
#include "data.hpp"

ising::ising(std::size_t sweeps, std::size_t width, unsigned num_neighbors,
             std::size_t nimg, unsigned scale, double temp, double fstr,
             std::string backend)
{
    this->sweeps = sweeps;
    data.E.reserve(sweeps);
    data.E_fluc.reserve(sweeps);
    data.M.reserve(sweeps);
    data.M_fluc.reserve(sweeps);

    this->width = width;
    this->num_neighbors = num_neighbors;
    this->scale = scale;

    switch(num_neighbors)
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
    neighbors nb = get_neighbors(i, j);

    return 2 * spins[i * width + j] * (nb.top + nb.bottom + nb.left + nb.right + field_strength);
}

neighbors ising::get_neighbors(unsigned i, unsigned j)
{
    neighbors nb;

    if (i == 0)
    {
        nb.top = get_spin(width - 1, j);
        nb.bottom = get_spin(i + 1, j);

    }
    else if (i == width - 1)
    {
        nb.top = get_spin(i - 1, j);
        nb.bottom = get_spin(0, j);
    }
    else
    {
        nb.top = get_spin(i-1, j);
        nb.bottom = get_spin(i+1, j);
    }

    if (j == 0)
    {
        nb.left = get_spin(i, width - 1);
        nb.right = get_spin(i, j + 1);
    }
    else if (j == width - 1)
    {
        nb.left = get_spin(i, j - 1);
        nb.right = get_spin(i, 0);
    }
    else
    {
        nb.left = get_spin(i, j - 1);
        nb.right = get_spin(i, j + 1);
    }

    int check = nb.left * nb.right * nb.top * nb.bottom;
    check = check * check;
    if (check != 1)
    {
        std::cerr << "Error: some neighbors have garbage values! " + std::to_string(check) << std::endl;
    }
    return nb;
}

__global__ void cudaCalcTotalU(int width, int *grid, double *interactions, double h)
{
    int id = threadIdx.x;
    if (id < width * width)
    {
        int y = id / width;
        int x = id % width;

        int top, bottom, left, right;

        int spin = grid[id];

        if (y == 0)
        {
            top = grid[(width-1)*width + x];
            bottom = grid[width + x];
        }
        else if (y == width - 1)
        {
            top = grid[(y-1) * width + x];
            bottom = grid[x];
        }
        else
        {
            top = grid[(y-1)*width + x];
            bottom = grid[(y+1)*width + x];
        }

        if (x == 0)
        {
            left = grid[y*width + width-1];
            right = grid[y*width + 1];
        }
        else if (x == width - 1)
        {
            left = grid[y*width+width-2];
            right = grid[y*width];
        }
        else
        {
            left = grid[y*width + x - 1];
            right = grid[y*width + x + 1];
        }
        printf("%d %d %d %d\n", top, bottom, left, right);
        assert(top * top + bottom * bottom + left * left * right * right == 4);
        printf("%d %d %d %d\n", top, bottom, left, right);

        interactions[id] = grid[id] * (h + top + left + (y == width-1) * bottom + (x == width-1) * right);
    }

}

__global__ void reduce(double *in)
{
    int id = threadIdx.x;

    auto stride = 1;
    int threads = blockDim.x;

    while (threads > 0)
    {
        if (id < threads)
        {
            auto idx1 = id * stride * 2;
            auto idx2 = idx1 + stride;
            in[idx1] += in[idx2];
        }

        stride <<= 1;
        threads >>= 1;
    }
}

double ising::calc_totalU()
{
    auto total = 0.;

    // Our Hamiltonian is H = -epsilon*Sum[(s_i)(s_j)] - h sum[s_i]
    // So at every cell site, you sum neighbor interactions

    // Normal Stencil:
    //    |
    //  --+
    //
    //  Right Edge Stencil:
    //    |
    //  --+--
    //
    //  Bottom Edge Stencil:
    //    |
    //  --+
    //    |
    //
//    #pragma omp parallel for
    for (auto i = 0; i < width; i++)
    {
        for (auto j = 0; j < width; j++)
        {
            double interactions = 0.;
            neighbors nb = get_neighbors(i,j);

            interactions += nb.top + nb.left + field_strength;

            if (i == width-1)
            {
                interactions += nb.bottom;
            }
            if (j == width-1)
            {
                interactions += nb.right;
            }

            total -= get_spin(i,j) * interactions;
        }
    }
    return total;
}

double ising::calc_totalM()
{
    auto total = 0.;

    for (auto i = 0; i < width; i++)
    {
        for (auto j = 0; j < width; j++)
        {
            total += get_spin(i,j);
        }
    }
    return total;
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

//    #pragma omp parallel for
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
    // We want to print `nimg` images, so we print every `iter/nimg` step. But
    // this isn't always an integer, so let's increase iter until it is.
    for (auto T = 0.; T < 6; T += 0.01)
    {
        initialize_spins();
        auto iters = sweeps * width;

        double avg_E = 0.;
        double avg_E_square = 0.;
        double avg_M = 0.;
        double avg_M_square = 0.;

        data.T.push_back(T);
        std::cerr << "Performing " << iters << " iterations...\n";
        for (auto t = 0; t < iters; t++)
        {
            auto i = dist(engine);
            auto j = dist(engine);
            auto dU = calc_deltaU(i, j);

            if (dU <= 0.)
            {
                //auto u1 = calc_totalU();
                flip_spin(i,j);
                //auto u2 = calc_totalU();

                //auto totald = (u2-u1) - dU;
                // TODO: diagnose why we are getting errors as bad as 2...
                //if (totald > 4)
                {
                //    throw std::runtime_error("Energy " + std::to_string(totald) + " doesn't match!");
                }
            }
            else
            {
                // If T approaches zero, then the boltzmann factor becomes infinitely small
                if (T != 0.)
                {
                    if (floatdist(engine) < std::exp(-dU / T))
                    {
                        flip_spin(i,j);
                    }
                }
            }

            // Equilibrium sweeps contribute to equilibirum averages
            if (t > 10*width*width)
            {
                cudaError_t err;
                int *d_spins;
                double *d_interactions;
                err = cudaMalloc(&d_spins, width*width*sizeof(int));
                if (err != cudaSuccess) std::cerr << cudaGetErrorString(err) << std::endl;
                err =cudaMalloc(&d_interactions, width*width*sizeof(double));
                if (err != cudaSuccess) std::cerr << cudaGetErrorString(err) << std::endl;
                err =cudaMemcpy(d_spins, spins.get(), width*sizeof(int), cudaMemcpyHostToDevice);
                if (err != cudaSuccess) std::cerr << cudaGetErrorString(err) << std::endl;
                err =cudaDeviceSynchronize();
                if (err != cudaSuccess) std::cerr << cudaGetErrorString(err) << std::endl;
                cudaCalcTotalU<<<1, width*width>>>(width, d_spins, d_interactions, field_strength);
                err =cudaDeviceSynchronize();
                if (err != cudaSuccess) std::cerr << cudaGetErrorString(err) << std::endl;
                // reduce<<<1, width*width>>>(d_interactions);
                err = cudaDeviceSynchronize();
                if (err != cudaSuccess) std::cerr << cudaGetErrorString(err) << std::endl;
                double E;
                err = cudaMemcpy(&E, d_interactions, sizeof(double), cudaMemcpyDeviceToHost);
                if (err != cudaSuccess) std::cerr << cudaGetErrorString(err) << std::endl;
                err =cudaFree(d_spins);
                if (err != cudaSuccess) std::cerr << cudaGetErrorString(err) << std::endl;
                err =cudaFree(d_interactions);
                if (err != cudaSuccess) std::cerr << cudaGetErrorString(err) << std::endl;
                // std::cerr << "E has value " << E << std::endl;
                avg_E += E;
                avg_E_square += E * E;
                auto M = calc_totalM();
                avg_M += M;
                avg_M_square += M * M;
            }

            /* Now, handle visualization
            if ((t+1) % width == 0)
            {
                data.E.push_back(calc_totalU());
                if (display_mode == 0)
                {
                    print_snapshot();
                }
                else
                {
                    char filename[256];
                    snprintf(filename, 256, "snapshot_%lux%lu_%.8lu.png", width, width, t / width);
                    save_png_snapshot(filename);
                }
                std::cerr << "Total energy at sweep " << t / width << " is " << data.E.at(t / width) << std::endl;
            }
            */
        }

        auto count = iters - 10 * width * width;
        avg_E /= count;
        avg_E_square /= count;
        avg_M /= count;
        avg_M_square /= count;

        data.E.push_back(avg_E);
        data.E_fluc.push_back(avg_E_square - avg_E * avg_E);
        data.M.push_back(avg_M);
        data.M_fluc.push_back(avg_M_square - avg_M * avg_M);
        std::cerr << "Equilibrium E: " << data.E.back() << std::endl;
    }

    std::ofstream f1("output_" + std::to_string(width) + ".dat");

    for (auto i = 0; i < data.T.size(); i++)
    {
        f1 << static_cast<double>(i)*0.05 << " " << data.E[i] << " " << data.E_fluc[i] << " " << data.M[i] << " " << data.M_fluc[i] << "\n";
    }
}

inline void ising::flip_spin(std::size_t i, std::size_t j)
{
    auto spin = get_spin(i,j);

    if (spin * spin != 1)
    {
        throw std::runtime_error("Invalid spin detected");
    }
    spins[i * width + j] = -spin;
}

int ising::get_spin(std::size_t i, std::size_t j) const
{
    return spins[i * width + j];
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
                    img[scale*y + i][scale*x + j] = png::gray_pixel_1(get_spin(y,x) > 0);
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
                    std::cout << (get_spin(y,x) > 0 ? '.' : '+');
                }
            }
            std::cout << '\n';
        }
    }
    std::cout << "\n\n";
}
