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

double ising::calc_deltaU(long x, long y)
{
    neighbors nb = get_neighbors(x, y);

    return 2 * get_spin(x,y) * (nb.top + nb.bottom + nb.left + nb.right + field_strength);
}

neighbors ising::get_neighbors(long x, long y)
{
    neighbors nb;

    // std::cerr << "Getting neighbors for (" << x << ',' << y << "):\n";
    nb.left = get_spin(x - 1, y);
    nb.right = get_spin(x + 1, y);
    nb.top = get_spin(x, y - 1);
    nb.bottom = get_spin(x, y + 1);
    // std::cerr << nb.left << " " << nb.right << " " << nb.top << " " << nb.bottom << "\n\n";

    int check = std::abs(nb.left * nb.right * nb.top * nb.bottom);
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
    // BOUNDARY CONDITIONS: remember we're on a 2d torus, so we shouldn't
    // double count the horizontal and vertical edges!

    // Stencil:
    //    |
    //  --+
    for (auto y = 0; y < width; y++)
    {
        for (auto x = 0; x < width; x++)
        {
            int interactions = 0;
            neighbors nb = get_neighbors(x,y);
            interactions += nb.top + nb.left;

            total -= get_spin(x,y) * (static_cast<double>(interactions) + field_strength);
        }
    }
    if (total < - (2 + field_strength) * width * width)
        throw std::runtime_error("Total energy calculation is wrong!");
    return total;
}

/*
double ising::calc_S(std::size_t n, double dT)
{
    double ll, l, r, rr;

    if (n == 0)
    {
        return 0.;
    }
}
*/

double ising::calc_totalM()
{
    auto total = 0.;

    for (auto y = 0; y < width; y++)
    {
        for (auto x = 0; x < width; x++)
        {
            total += get_spin(x,y);
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
    total_U = calc_totalU();
    total_M = calc_totalM();
    std::cerr << "Starting E = " << total_U << " M = " << total_M << "\n"; 
}

void ising::run()
{
    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<int> dist(0, width - 1);
    std::uniform_real_distribution<double> floatdist;
    // We want to print `nimg` images, so we print every `iter/nimg` step. But
    // this isn't always an integer, so let's increase iter until it is.
    for (auto T = 0.; T < 6; T += 0.05)
    {
        // Reuse previous state:
        if (sanitize)
            initialize_spins();

        auto iters = sweeps * width * width;

        double cU = 0.;

        double avg_E = 0.;
        double cE = 0.;
        double avg_E_square = 0.;

        double avg_M = 0.;
        double cM = 0.;
        double avg_M_square = 0.;

        std::size_t count = 0;

        data.cf.push_back(std::vector<double>(width / 2, 0.));
        data.T.push_back(T);
        std::cerr << "Performing " << iters << " iterations...\n";
        std::vector<double> temp1;
        std::vector<double> temp2;
        for (auto t = 0; t < iters; t++)
        {
            auto x = dist(engine);
            auto y = dist(engine);
            auto dU = calc_deltaU(x, y);

            if (std::abs(dU) > 2 * (4 + field_strength))
            {
                throw std::runtime_error("Energy growing too fast!");
            }

            if (dU <= 0.)
            {
                //auto u1 = calc_totalU();
                flip_spin(x,y);
                //auto u2 = calc_totalU();
                double yU = dU - cU;
                double tU = total_U + yU;
                cU = (tU - total_U) - yU;
                total_U = tU;
                total_M += 2*get_spin(x,y);

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
                        flip_spin(x,y);
                        double yU = dU - cU;
                        double tU = total_U + yU;
                        cU = (tU - total_U) - yU;
                        total_U = tU;
                        total_M += 2*get_spin(x,y);
                    }
                }
            }

            if (!(t % width * width) && t > 100*width * width)
            {
                temp1.push_back(total_U);
                temp2.push_back(total_M);
                // Equilibrium sweeps contribute to equilibirum averages
                auto yE = total_U - cE;
                auto tE = avg_E + yE;
                cE = (tE - avg_E) - yE;
                avg_E = tE;

                auto yM = total_M - cM;
                auto tM = avg_M + yM;
                cM = (tM - avg_M) - yM;
                avg_M = tM;

                /*
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
        std::cerr << "avg_E is now " << avg_E << '\n';
        avg_M /= count;


        cE = 0.;
        cM = 0.;
        for (auto i = 0; i < temp1.size(); i++)
        {
            double a = temp1[i] - avg_E;
            double b = temp2[i] - avg_M;
            a *= a;
            b *= b;
            if (a < 0. || b < 0.)
                throw std::runtime_error("squaring overflow?");

            auto yE = a - cE;
            auto tE = avg_E_square + yE;
            cE = (tE - avg_E_square) - yE;
            avg_E_square = tE;

            auto yM = b - cM;
            auto tM = avg_M_square + yM;
            cM = (tM - avg_M_square) - yM;
            avg_M_square = tM;
        }
        avg_E_square /= count;
        avg_M_square /= count;

        /*
        for (auto i = 0; i != data.cf.back().size(); i++)
        {
            data.cf.back()[i] /= count;
        }
        */

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
        data.E_fluc.push_back(avg_E_square);
        data.avg_M.push_back(avg_M);
        data.M_fluc.push_back(avg_M_square);
        // The heat capacity is the energy fluctiation divided by Temp^2:
        data.C_v.push_back( T == 0. ? NAN : data.M_fluc.back() / (T * T));
        // data.cf_len.push_back(corr_len);
        std::cerr << "Equilibrium E: " << avg_E << std::endl;
        avg_E = 0.;
    }

    std::ofstream f1("output_" + std::to_string(width) + ".dat");

    for (auto i = 0; i < data.T.size(); i++)
    {
        f1 << static_cast<double>(i)*0.05 << " " << data.avg_E[i] << " " << data.E_fluc[i] << " " << data.avg_M[i] << " " << data.M_fluc[i] << " " << data.S[i] << " " << data.C_v[i] << "\n";
    }

    std::ofstream f2("sweeps_" + std::to_string(width) + ".dat");

    for (auto i = 0; i < data.E_sweeps.size(); i++)
    {
        f2 << i << " " << data.E_sweeps[i] << " " << data.M_sweeps[i] << "\n";
    }

    std::ofstream f3("cfs_" + std::to_string(width) + ".dat");

    for (auto j = 0; j < data.cf.back().size(); j++)
    {
       f3 << j * 0.05;
       for (auto k = 0; k < data.cf.size(); k++)
       {
           f3 << " " << data.cf[k][j];
       }
       f3 << "\n";
    }
}

inline void ising::flip_spin(long x, long y)
{
    auto spin = get_spin(x,y);

    if (spin * spin != 1)
    {
        throw std::runtime_error("Invalid spin detected");
    }
    spins[y * width + x] = -spin;
}

char ising::get_spin(long x, long y) const
{
    x = x % width;
    y = y % width;
    // We will manually enforce periodic BCs here:
    if (x < 0)
    {
        // std::cerr << "Correcting value " << x;
        x += width;
        // std::cerr << " to " << x << '\n';
    }
    if (y < 0)
    {
        // std::cerr << "Correcting value " << y;
        y += width;
        // std::cerr << " to " << y << '\n';
    }

    return spins[y * width + x];
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

void ising::set_sanitizer(bool flag)
{
    sanitize = flag;
}
