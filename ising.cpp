// C++ Standard Libraries
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

// System Libraries
#include <png++/png.hpp>
#include <omp.h>

// Project Libraries
#include "ising.hpp"
#include "data.hpp"

ising::ising(std::size_t sweeps, std::size_t width, unsigned num_neighbors,
             std::size_t nimg, unsigned scale, double temp, double fstr,
             std::string backend)
{
    this->sweeps = sweeps;
    /*
    data.E.reserve(sweeps);
    data.E_fluc.reserve(sweeps);
    data.M.reserve(sweeps);
    data.M_fluc.reserve(sweeps);
    */

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
    total_U = calc_totalU();
    total_M = calc_totalM();

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

double ising::calc_totalU()
{
    auto total = 0.;

    // Our Hamiltonian is H = -epsilon*Sum[(s_i)(s_j)] - h sum[s_i]
    // So at every cell site, you sum neighbor interactions
    // Let's ignore boundary conditions and only calculate grid energies

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
        // Reuse previous state:
        // initialize_spins();
        auto iters = sweeps * width;

        double avg_E = 0.;
        double avg_E_square = 0.;
        double avg_M = 0.;
        double avg_M_square = 0.;

        std::size_t count = 0;

        data.cf.push_back(std::vector<double>(width / 2, 0.));
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
                total_U += dU;
                total_M += 2*get_spin(i,j);

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
                        total_U += dU;
                        total_M += 2*get_spin(i,j);
                    }
                }
            }

            if (!(t % width * width) && t > 100*width * width)
            {

                // Equilibrium sweeps contribute to equilibirum averages
                avg_E += total_U;
                avg_E_square += total_U * total_U;
                avg_M += total_M;
                avg_M_square += total_M * total_M;

                auto avg_spin = 0.;

                for (auto i = 0u; i < width; i++)
                {
                    for (auto j = 0u; j < width; j++)
                    {
                        avg_spin += get_spin(i,j);
                    }
                }

                avg_spin /= width * width;

                for (auto i = 0; i < width/2; i++)
                {
                    double correlation = 0.;
                    int count = 0;
                    for (auto y = 0; y < width/2 - i; y++)
                    {
                        for (auto x = 0; x < width/2 - i; x++)
                        {
                            correlation += get_spin(y,x) * (get_spin(y,x+i) + get_spin(y+i, x)); 
                            count++;
                        }
                    }
                    correlation = correlation / count - avg_spin * avg_spin;
                    data.cf.back()[i] += correlation;
                }

                /*
                if (t >  iters - 5000 * width * width)
                {
                    data.M_sweeps.push_back(M);
                    data.E_sweeps.push_back(E);
                }
                */
                count++;
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

        avg_E /= count;
        avg_E_square /= count;
        avg_M /= count;
        avg_M_square /= count;
        for (auto i = 0; i != data.cf.back().size(); i++)
        {
            data.cf.back()[i] /= count;
        }

        /*
        auto start_i = 0;
        unsigned corr_len;
        for (auto i = start_i; i != data.cf.back().size(); i++)
        {
            if (data.cf.back()[start_i] == 0.)
            {
                start_i++;
                continue;
            }
            double diff = data.cf.back()[i]/data.cf.back()[start_i];

            if (diff >= 2.718)
            {
               corr_len = (i - start_i)+1;
            }
        }*/

        // The entropy at T=0 is zero, but at any other state we get it by forward-derivatives:
        data.S.push_back(T == 0. ? 0. : data.S.back() + (avg_E - data.avg_E.back())/ T); 
        data.avg_E.push_back(avg_E);
        data.E_fluc.push_back(avg_E_square - avg_E * avg_E);
        data.avg_M.push_back(avg_M);
        data.M_fluc.push_back(avg_M_square - avg_M * avg_M);
        // The heat capacity is the energy fluctiation divided by Temp^2:
        data.C_v.push_back( T == 0. ? NAN : data.M_fluc.back() / (T * T));
        // data.cf_len.push_back(corr_len);
        std::cerr << "Equilibrium E: " << data.avg_E.back() << std::endl;
    }

    std::ofstream f1("output_" + std::to_string(width) + ".dat");

    for (auto i = 0; i < data.T.size(); i++)
    {
        f1 << static_cast<double>(i)*0.01 << " " << data.avg_E[i] << " " << data.E_fluc[i] << " " << data.avg_M[i] << " " << data.M_fluc[i] << " " << data.S[i] << " " << data.C_v[i] << "\n";
    }

    std::ofstream f2("sweeps_" + std::to_string(width) + ".dat");

    for (auto i = 0; i < data.E_sweeps.size(); i++)
    {
        f2 << i << " " << data.E_sweeps[i] << " " << data.M_sweeps[i] << "\n";
    }

    std::ofstream f3("cfs_" + std::to_string(width) + ".dat");

    for (auto j = 0; j < data.cf.back().size(); j++)
    {
       f3 << j << " " << data.cf[15][j] << "\n";
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
