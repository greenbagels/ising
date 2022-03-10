/*! @file grid.hpp
 *  @brief Grid class interface file
 *  @author Sameed Pervaiz (pervaiz.8@osu.edu)
 *  @copyright GPLv3
 *  @date 2021-05-08
 */

#ifndef GRID_HPP
#define GRID_HPP

#include <initializer_list>
#include <vector>

namespace ising
{
    class grid
    {
        // Internal typedefs
        // Underlying data types
        typedef int value_type;
        typedef int& reference;
        typedef (const int&) const_reference;
        // The underlying container
        typedef std::vector<value_type> container;
        // Iterators
        typedef container::iterator iterator;
        typedef container::const_iterator const_iterator;
        typedef container::reverse_iterator reverse_iterator;
        typedef container::const_reverse_iterator const_reverse_iterator;

        public:
            grid(std::vector<int> width_list);
            value_type dim() const;
            value_type width(int i) const;
            value_type size() const;
            // TODO: implement resizing
            // void resize(std::initializer_list<int> new_widths);

            /*! Element-access member functions */
            reference operator[](std::vector<int> coords);
            const_reference operator[](std::vector<int> coords) const;
            reference at(std::vector<int> coords);
            const_reference at(std::vector<int> coords) const;

            /*! Iterator member functions */
            /*! Forward iterators */
            iterator begin() { return data.begin(); }
            const_iterator begin() const { return data.begin(); }
            const_iterator cbegin() const { return data.cbegin(); }
            iterator end() { return data.end(); }
            const_iterator end() const { return data.end(); }
            const_iterator cend() const { return data.cend(); }

            /*! Reverse iterators */
            iterator rbegin() { return data.rbegin(); }
            const_iterator rbegin() const { return data.rbegin(); }
            const_iterator crbegin() const { return data.crbegin(); }
            iterator rend() { return data.rend(); }
            const_iterator rend() const { return data.rend(); }
            const_iterator crend() const { return data.crend(); }
        private:
            /* Contains the width of each of the grid's dimensions */
            container widths;
            /* Contains the actual grid elements, of size widths[0] * ... * widths[n-1] */
            container data;
    };
}
#endif
