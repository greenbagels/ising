/*! @file grid.hpp
 *  @brief Grid class implementation file
 *  @author Sameed Pervaiz (pervaiz.8@osu.edu)
 *  @copyright GPLv3
 *  @date 2021-05-08
 */

#include "grid.hpp"
#include <functional>
#include <numeric>

namespace ising
{
    grid::grid(std::vector<int> width_list)
        : widths(width_list), data(std::accumulate(widths.begin(), widths.end(), 1, std::multiplies))
    {
    }

    int grid::dim() const
    {
        return widths.size();
    }

    int grid::width(int i) const
    {
        return widths.at(i);
    }

    int grid::size() const
    {
        return data.size();
    }
    
    int& grid::operator[](std::vector<int> coords)
    {
        /* Suppose the coords are {x, y, z} and the grid dims are {dimx, dimy, dimz}.
         * Then the 1D index for the element at position (x,y,z) is equal to
         * x + y * dimx + z * dimx * dimy. In general, this is just the sum of
         * coord_i * (dim_0 * ... dim_(i-1)) over i, which we can iteratively calculate. */
        auto index = 0;
        auto dim_weight = 1;
        for (auto i = 0; i < coords.size(); i++)
        {
            index += coords[i] * dim_weight;
            scalar *= widths[i];
        }
        return data[index]; 
    }

    const int& grid::operator[](std::vector<int> coords) const
    {
        auto index = 0;
        auto dim_weight = 1;
        for (auto i = 0; i < coords.size(); i++)
        {
            index += coords[i] * dim_weight;
            scalar *= widths[i];
        }
        return data[index]; 
    }

    int& grid::at(std::vector<int> coords)
    {
        auto index = 0;
        auto dim_weight = 1;
        for (auto i = 0; i < coords.size(); i++)
        {
            index += coords.at(i) * dim_weight;
            scalar *= widths.at(i);
        }
        return data.at(i); 
    }
    
    const int& grid::at(std::vector<int> coords) const
    {
        auto index = 0;
        auto dim_weight = 1;
        for (auto i = 0; i < coords.size(); i++)
        {
            index += coords.at(i) * dim_weight;
            scalar *= widths.at(i);
        }
        return data.at(i); 
    }
}
