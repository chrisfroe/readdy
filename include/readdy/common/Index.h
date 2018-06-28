/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * This header file contains definitions for the Index class, which makes the access to contiguously stored multi-
 * dimensional arrays easier.
 *
 * @file Index.h
 * @brief Definitions of the Index class, making contiguously stored multi-dimensional array access easier.
 * @author clonker
 * @date 11.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <initializer_list>
#include <array>
#include <numeric>
#include <tuple>
#include <cmath>
#include "macros.h"
#include "tuple_utils.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)

NAMESPACE_BEGIN(detail)
template <typename T1, typename... T>
struct variadic_first {
    /**
     * type of the first element of a variadic type tuple
     */
    using type = typename std::decay<T1>::type;
};
template<bool B, typename T = void> using disable_if = std::enable_if<!B, T>;

template<std::size_t N, typename T, typename Tuple>
struct all_of_type_impl {
    /**
     * type of the n-th element of a packed variadic type tuple
     */
    using nth_elem = typename std::decay<typename std::tuple_element<N-1, Tuple>::type>::type;
    /**
     * if this is true, all elements in the parameter pack were of type T
     */
    static constexpr bool value = std::is_same<nth_elem, T>::value && all_of_type_impl<N-1, T, Tuple>::value;
};

template<typename T, typename Tuple>
struct all_of_type_impl<0, T, Tuple> {
    /**
     * this is true for empty parameter packs
     */
    static constexpr bool value = true;
};

template<typename T, typename ...Args>
struct all_of_type {
    /**
     * if this is true, all elements of the parameter pack are of type T
     */
    static constexpr bool value = all_of_type_impl<sizeof...(Args), T, typename std::tuple<Args...>>::value;
};

NAMESPACE_END(detail)

template<std::size_t Dims>
class Index {
    static_assert(Dims > 0, "Dims has to be > 0");
public:
    /**
     * Type that holds the dimensions of the index grid
     */
    using GridDims = std::array<std::size_t, Dims>;
    /**
     * The value type, inherited from GridDims::value_type
     */
    using value_type = typename GridDims::value_type;

    /**
     * Constructs an empty index object of specified dimensionality. Not of much use, really.
     */
    Index() : _size(), n_elems(0) {}

    /**
     * we can copy this
     */
    Index<Dims>(const Index<Dims> &) = default;
    /**
     * we can copy assign this
     */
    Index<Dims> &operator=(const Index<Dims> &) = default;
    /**
     * we can move this
     */
    Index<Dims>(Index<Dims> &&) = default;
    /**
     * we can move assign this
     */
    Index<Dims> &operator=(Index<Dims> &&) = default;

    /**
     * Constructs an index object with a number of size_t arguments that must coincide with the number of dimensions,
     * specifying the grid.
     * @tparam Args the argument types, must all be size_t
     * @param args the arguments
     */
    template<typename ...Args, typename = typename std::enable_if<detail::all_of_type<std::size_t, Args...>::value>::type>
    Index(Args &&...args) : _size({std::forward<Args>(args)...}),
                            n_elems(std::accumulate(_size.begin(), _size.end(), static_cast<std::size_t>(1), std::multiplies<value_type>())) {}

    /**
     * the number of elements in this index, exactly the product of the grid dimensions
     * @return the number of elements
     */
    value_type size() const {
        return n_elems;
    }

    /**
     * the number of elements in this index
     * @return the number of elements
     */
    value_type nElements() const {
        return n_elems;
    }

    /**
     * Retrieve size of N-th axis
     * @tparam N the axis
     * @return size of N-th axis
     */
    template<int N>
    constexpr value_type get() const {
        return _size[N];
    }

    /**
     * retrieve size of N-th axis
     * @param N N
     * @return size of N-th axis
     */
    template<typename T>
    constexpr value_type operator[](T N) const {
        return _size[N];
    }

    /**
     * map Dims-dimensional index to 1D index
     * @tparam Ix the d-dimensional index template param type
     * @param ix the d-dimensional index
     * @return the 1D index
     */
    template<typename... Ix>
    constexpr value_type operator()(Ix &&... ix) const {
        static_assert(sizeof...(ix) == Dims, "wrong input dim");
        /*if(n_elems > 0)*/ {
            std::array<typename detail::variadic_first<Ix...>::type, Dims> indices{std::forward<Ix>(ix)...};
            std::size_t result = 0;
            auto prefactor = n_elems / _size[0];
            for (std::size_t d = 0; d < Dims - 1; ++d) {
                result += prefactor * indices[d];
                prefactor /= _size[d + 1];
            }
            result += indices[Dims - 1];
            return result;
        }
        return 0;
    }

    /**
     * Inverse mapping 1D index to Dims-dimensional tuple
     * @param idx
     * @return
     */
    GridDims inverse(std::size_t idx) const {
        GridDims result;
        auto prefactor = n_elems / _size[0];
        for(std::size_t d = 0; d < Dims-1; ++d) {
            auto x = std::floor(idx / prefactor);
            result[d] = x;
            idx -= x * prefactor;
            prefactor /= _size[d+1];
        }
        result[Dims-1] = idx;
        return result;
    }

private:

    GridDims _size;
    value_type n_elems;
};

/**
 * Definition of a 1D index object, basically an identity operator
 */
using Index1D = Index<1>;
/**
 * Definition of a 2D index object
 */
using Index2D = Index<2>;
/**
 * Definition of a 3D index object
 */
using Index3D = Index<3>;

NAMESPACE_END(util)
NAMESPACE_END(readdy)
