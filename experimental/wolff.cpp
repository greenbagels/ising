
#include <iostream>
#include <cmath>
#include <queue>
#include <tuple>
#include <utility>
#include <png++/png.hpp>
#include "wolff.hpp"

wolff::wolff(int n, int w, int pix_scale, double temp)
    : width(w), sweeps(n), scale(pix_scale), T(temp), engine(6)
{
    grid = new int[width*width];
    randomize_grid();
}

void wolff::randomize_grid()
{
    for (auto y = 0; y < width; y++)
    {
        for (auto x = 0; x < width; x++)
        {
            *get_spin(x, y) = rand_int(0,1);
        }
    }
}

int wolff::rand_int(int a, int b)
{
    std::uniform_int_distribution<int> dist(a,b);
    return dist(engine);
}

double wolff::rand_dbl(double a, double b)
{
    std::uniform_real_distribution<double> dist(a, b);
    return dist(engine);
}

int* wolff::get_spin(int x, int y)
{
    // todo: optimize this later
    while (x < 0)
    {
        x += width;
    }
    while (y < 0)
    {
        y += width;
    }
    while (x >= width)
    {
        x -= width;
    }
    while (y > width)
    {
        y -= width;
    }

    return &grid[y * width + x];
}

void wolff::iterate()
{
    std::queue<std::pair<int, int>> Q;

    int randx = rand_int(0, width-1);
    int randy = rand_int(0, width-1);
    // std::cerr << "Generated random coordinates (" << randx << "," << randy << ").\n";
    int parent_spin = *get_spin(randx, randy);
    // std::cerr << "Parent spin state is " << parent_spin << ". Starting probabilistic"
    //    " flood fill now.";
    // We always flip the initial spin
    Q.emplace(randx, randy);

    while (!Q.empty())
    {
        int x, y;
        // These pairs correspond to (right, down, left, up)
        int dx[4] = {+1, 0, -1, 0};
        int dy[4] = {0, +1, 0, -1};
        //const char* dirs[4] = {"east", "south", "west", "north"};
        std::tie(x, y) = Q.front();
        //std::cerr << "Popped coordinate (" << x << "," << y << ").\n";
        Q.pop();
        // the same node might get queued multiple times before it's handled
        // we can't just negate the nodes value, since we could negate it multiple
        // times. so instead, we always assign it the proper negated spin using
        // the parent value!
        *get_spin(x,y) = !parent_spin;
        for (auto i = 0; i < 4; i++)
        {
            //std::cerr << "Checking " << dirs[i] << " neighbor.\n";
            auto neighborx = x+dx[i];
            auto neighbory = y+dy[i];
            auto neighborspin = *get_spin(neighborx, neighbory);
            if (neighborspin == parent_spin)
            {
                //std::cerr << "Neighbor matched parent spin!\n";
                // Now, flip with the proper probability
                // beta = 1/(kB * T);
                // T is defined in terms of J/kB, so the numerical value of 1/T
                // ends up being the value of beta in these units.
                auto Jbeta = 1. / T;

                // T is defined in terms of J/kB, so the numerical value of 1/T
                // ends up being the value of beta in these units.

                // The argument of exp needs the product of spin sites,
                // but this is always going to be 1 if the only the active links
                // are to be handled
                auto p = 1 - std::exp(-2*Jbeta);
                if (rand_dbl(0., 1.) < p)
                {
                    Q.emplace(neighborx, neighbory);
                }
            }
        }
    }
}


void wolff::run()
{
    while (sweeps--)
    {
        iterate();
    }
}

void wolff::print_to_file(std::string filename)
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
