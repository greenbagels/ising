/*! @file sim.cpp
 *  @brief Ising model simulation class implementation
 *  @author Sameed Perviaz
 *  @copyright GPLv3
 *  @date 2021-05-04
 */

#include <iostream>
#include <cmath>
#include <queue>
#include <tuple>
#include <utility>
#include <png++/png.hpp>
#include "misc.hpp"
#include "sim.hpp"

namespace ising
{
    simulation::simulation(int d, int w, double temp, double field, bool randomize, int backend)
        : dim(d), width(w) spin_grid(, T(temp), engine(randomize ? std::random_device()() : 6)
    {
        // The class should be initialized in a usable state. The ctor arguments
        // ask for everything required to make this possible.
        auto gridsize = width;
        for (auto i = 1; i < dim; i++)
        {
            gridsize *= width;
        }
        grid = std::vector<int>(gridsize);
        randomize_grid();
        calculate_M();
        calculate_U();
        // std::cerr << "Starting M: " << M << "\nStarting U: " << U << std::endl;
        if (backend < 0 || backend >= NUM_BACKENDS)
        {
            throw std::runtime_error("Error: backend not supported!");
        }
        algo_backend = backend;
    }

    void simulation::randomize_grid()
    {
        for (auto i = spin_grid.begin(); i != spin_grid.end(); i++)
        {
            *i = rand_int(0,1);
        }
    }

    void simulation::set_T(double new_temp)
    {
        T = new_temp;
    }

    void simulation::calculate_M()
    {
        // This is just the total of the spins.
        auto total_M = 0.;
        for (auto i = 0; i < width; i++)
        {
            for (auto j = 0; j < width; j++)
            {
                // To optimize, we could keep our (0,1) spin representation,
                // and just total the 1s, and add -1(width*width-total). But
                // this would hide any errors in our spin flipping routines,
                // so for now, we explicitly calculate the sum over all spins.
                total_M += 2 * spin_grid[{i, j}] - 1;
            }
        }
        M = total_M;
    }

    void simulation::calculate_U()
    {
        // Because of periodic boundary conditions, we only need to check the energy
        // of the right neighbor and bottom neighbor interactions (since our top and
        // left neighbors will already include us this way!)
        double total_U = 0.;
        for (auto i = 0; i < width; i++)
        {
            for (auto j = 0; j < width; j++)
            {
                // In this case, we're using pixel coordinates [(0,0) = top left corner]
                // So checking bottom and right corresponds to (x+1,y) and (x,y+1).
                // We're also mapping (0,1) to (-1,1) with f(x)=2*x-1
                auto site_spin = 2 * spin_grid[{i, j}] - 1;
                auto right_spin = 2 * spin_grid[{i + 1, j}] - 1;
                auto bot_spin = 2 * spin_grid[{i, j + 1}] - 1;
                total_U -= site_spin * (right_spin + bot_spin);
            }
        }
        U = total_U;
    }

    double simulation::calc_cf(int r)
    {
        // The average spin is just the Magnetization / (width * width), so
        // we calculate <spin1 * spin2> - <spin>, where spin1 and spin2 are
        // a distance r apart:
        auto num_pairs = 2 * width * width;
        double avg_spin_sq = M / (width * width);
        avg_spin_sq *= avg_spin_sq;
        double cf = 0.;
        for (auto y = 0; y < width; y++)
        {
            for (auto x = 0; x < width; x++)
            {
                auto spin = 2*(*get_spin(x,y))-1;
                auto spinR = 2*(*get_spin(x+r,y))-1;
                // auto spinL = 2*(*get_spin(x-r,y))-1;
                // auto spinU = 2*(*get_spin(x,y-r))-1;
                auto spinD = 2*(*get_spin(x,y+r))-1;
                cf += spin * (spinR + /*spinL + spinU +*/ spinD);
            }
        }
        return cf / num_pairs - avg_spin_sq;
    }

    double simulation::get_M() const
    {
        return M;
    }

    double simulation::get_U() const
    {
        return U;
    }

    double simulation::get_T() const
    {
        return T;
    }

    int simulation::rand_int(int a, int b)
    {
        std::uniform_int_distribution<int> dist(a,b);
        return dist(engine);
    }

    double simulation::rand_dbl(double a, double b)
    {
        std::uniform_real_distribution<double> dist(a, b);
        return dist(engine);
    }

    int* simulation::get_spin(int x, int y)
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

        return &grid[y * width + x];
    }

    int simulation::single_spin_flip()
    {
        int randx = rand_int(0, width-1);
        int randy = rand_int(0, width-1);
        int parent_spin = 2*(*get_spin(randx, randy))-1;
        auto deltaU = 2 * parent_spin * (2 * (*get_spin(randx - 1, randy)
                                            + *get_spin(randx + 1, randy)
                                            + *get_spin(randx,     randy - 1)
                                            + *get_spin(randx,     randy + 1)) - 4);
        if (deltaU <= 0)
        {
            *get_spin(randx, randy) ^= 1;
            return 1;
        }
        else
        {
            if (rand_dbl(0., 1.) < std::exp(-deltaU / T))
            {
                *get_spin(randx, randy) ^= 1;
                return 1;
            }
        }

        return 0;
    }

    int simulation::wolff_cluster_flip()
    {
        // We need to store coordinates rather than pointers, so we can check our neighbors
        // without fearing segfaults
        std::queue<std::pair<int, int>> Q;
        int spins_flipped = 0;

        int randx = rand_int(0, width-1);
        int randy = rand_int(0, width-1);
        // std::cerr << "Generated random coordinates (" << randx << "," << randy << ").\n";
        int parent_spin = *get_spin(randx, randy);
        //std::cerr << "Parent spin state is " << parent_spin << ". Starting probabilistic"
        //    " flood fill now.";
        // We always flip the initial spin
        Q.emplace(randx, randy);
        while (!Q.empty())
        {
            int x, y;
            // These pairs correspond to (right, down, left, up)
            int dx[4] = {+1, 0,-1, 0};
            int dy[4] = { 0,+1, 0,-1};
            //const char* dirs[4] = {"east", "south", "west", "north"};
            std::tie(x, y) = Q.front();
            //std::cerr << "Popped coordinate (" << x << "," << y << ")." << std::endl;
            Q.pop();
            auto spin = get_spin(x,y);
            if (*spin == parent_spin)
            {
                // check our four neighbors iteratively
                for (auto i = 0; i < 4; i++)
                {
                    //std::cerr << "Checking " << dirs[i] << " neighbor.\n";
                    auto neighborx = x+dx[i];
                    auto neighbory = y+dy[i];
                    while (neighborx < 0)
                    {
                        neighborx += width;
                    }
                    while (neighbory < 0)
                    {
                        neighbory += width;
                    }
                    while (neighborx >= width)
                    {
                        neighborx -= width;
                    }
                    while (neighbory >= width)
                    {
                        neighbory -= width;
                    }
                    auto neighborspin = get_spin(neighborx, neighbory);
                    // First, update the energy change, since we know the current
                    // child site is being flipped.
                    // std::cerr << "Spin is " << 2*(*spin)-1
                    // << ". Neighbor is " << 2*(*neighborspin)-1
                    // << ". Change in energy is "
                    // << 2 * (2*(*spin)-1) * (2*(*neighborspin)-1) << "\n";
                    U += 2 * (2*(*spin)-1) * (2*(*neighborspin)-1);
                    // Now, add to the cluster (and flip) with the proper probability
                    // beta = 1/(kB * T);
                    // T is defined in terms of J/kB, so the numerical value of 1/T
                    // ends up being the value of beta in these units.
                    auto Jbeta = 1. / T;

                    // The argument of exp needs the product of spin sites, assuming
                    // +1/-1 spins. But our spins are 0 or 1, so we have to compare
                    // the spins instead of just multiplying them.
                    auto p = 1 - std::exp(-2*Jbeta*(*spin==*neighborspin));
                    // God doesn't play dice, but that won't stop us!
                    auto roll = rand_dbl(0., 1.);
                    if (roll < p)
                    {
                        //std::cerr << "Pushing coords (" << neighborx << "," << neighbory << ").\n";
                        Q.emplace(neighborx, neighbory);
                    }
                }
                // Last but not least, don't forget to update M and flip our state.
                *spin = !parent_spin;
                spins_flipped++;
                // If we move from -1 to 1, our M changes by 2. Similarly, going from
                // 1 to -1, our M changes by -2. Both of these changes are equal to
                // -2*before = 2*after. Let's just go with the after one.
                M += 2*(2*(*spin)-1);
            }
        }
        return spins_flipped;
    }

    int simulation::iterate()
    {
        return spin_algos[algo_backend]();
    }
    void simulation::run(int sweeps)
    {
        while (sweeps--)
        {
            iterate();
        }
    }

    void simulation::print_to_file(std::string filename, int scale)
    {
        png::image<png::rgb_pixel> img(scale * width, scale * width);

        for (auto y = 0; y < scale * width; y++)
        {
            for (auto x = 0; x < scale * width; x++)
            {
                // Our spin location to pixel mapping is 1:scale*scale; so the state
                // at (x,y) gets scaled to to the block
                // [x*scale, x*scale + scale) x [y*scale, y*scale+scale).

                // As a (potentially premature) optimization, we're going to write
                // the image data one horizontal line at a time.
                auto spin = *get_spin(x / scale, y / scale);
                img[y][x] = png::rgb_pixel(255 * spin, 255 * spin, 255 * spin);
            }
        }
        img.write(filename);
    }
}
