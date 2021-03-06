/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          * 
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
 * << detailed description >>
 *
 * @file tuple_utils.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <functional>

#include "macros.h"
#include "index_sequence.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(detail)

template<typename R, template<typename...> class Params, typename... Args, std::size_t... I>
R call_helper(std::function<R(Args...)> const &func, Params<Args...> const &params, std::index_sequence<I...>) {
    return func(std::get<I>(params)...);
}

template<typename R, template<typename...> class Params, typename... Args, std::size_t... I>
R call_helper(R (*func)(Args...), Params<Args...> const &params, std::index_sequence<I...>) {
    return func(std::get<I>(params)...);
}

template<typename T, typename TT = typename std::remove_reference<T>::type, size_t... I>
auto reverse_impl(T &&t, std::index_sequence<I...>) -> std::tuple<typename std::tuple_element<
        sizeof...(I) - 1 - I, TT>::type...> {
    return std::make_tuple(std::get<sizeof...(I) - 1 - I>(std::forward<T>(t))...);
}

template<class F, class...Ts, std::size_t...Is>
void for_each_in_tuple_impl(const std::tuple<Ts...> &tuple, F func, std::index_sequence<Is...>) {
    using expander = int[];
    (void) expander {0, ((void) func(std::get<Is>(tuple)), 0)...};
}

template<class F, class...Ts, std::size_t...Is>
void for_each_in_tuple_reverse_impl(const std::tuple<Ts...> &tuple, F func, std::index_sequence<Is...>) {
    using expander = int[];
    (void) expander {0, ((void) func(std::get<sizeof...(Is) - 1 - Is>(tuple)), 0)...};
}

template<class Ch, class Tr, class Tuple, std::size_t... Is>
void print_tuple_impl(std::basic_ostream<Ch, Tr> &os, const Tuple &t, std::index_sequence<Is...>) {
    using swallow = int[]; // guarantees left to right order
    (void) swallow{0, (void(os << (Is == 0 ? "" : ", ") << std::get<Is>(t)), 0)...};
}


NAMESPACE_END(detail)

/**
 * Invoke the callable object func with a tuple of arguments.
 * @tparam R the return type
 * @tparam Params argument pack
 * @tparam Args arguments
 * @param func the function
 * @param params parameter tuple
 * @return the return value of func(Args...)
 */
template<typename R, template<typename...> class Params, typename... Args>
R call(std::function<R(Args...)> const &func, Params<Args...> const &params) {
    return detail::call_helper(func, params, std::index_sequence_for<Args...>{});
}

/**
 * Invoke the callable object func with a tuple of arguments.
 * @tparam R the return type
 * @tparam Params argument pack
 * @tparam Args arguments
 * @param func the function
 * @param params parameter tuple
 * @return the return value of func(Args...)
 */
template<typename R, template<typename...> class Params, typename... Args>
R call(R (*func)(Args...), Params<Args...> const &params) {
    return detail::call_helper(func, params, std::index_sequence_for<Args...>{});
}

/**
 * Reverse a tuple
 * @tparam T the tuple type
 * @tparam TT tuple type stripped of reference
 * @param t the tuple
 * @return reversed tuple
 */
template<typename T, typename TT = typename std::remove_reference<T>::type>
auto reverse(T &&t) -> decltype(detail::reverse_impl(std::forward<T>(t),
                                                     std::make_index_sequence<std::tuple_size<TT>::value>())) {
    return detail::reverse_impl(std::forward<T>(t),
                                std::make_index_sequence<std::tuple_size<TT>::value>());
}

/**
 * map a function to each element in a tuple starting with the last element
 * @tparam F the function type
 * @tparam Ts the tuple type
 * @param tuple the tuple
 * @param func the function
 */
template<class F, class... Ts>
void for_each_in_tuple_reverse(const std::tuple<Ts...> &tuple, F func) {
    return detail::for_each_in_tuple_reverse_impl(tuple, func, std::make_index_sequence<sizeof...(Ts)>());
}


/**
 * map function to each element in a tuple
 * @tparam F the function type
 * @tparam Ts the tuple type
 * @param tuple the tuple
 * @param func the function
 */
template<class F, class...Ts>
void for_each_in_tuple(const std::tuple<Ts...> &tuple, F func) {
    detail::for_each_in_tuple_impl(tuple, func, std::make_index_sequence<sizeof...(Ts)>());
}

NAMESPACE_END(util)
NAMESPACE_END(readdy)

/**
 * pretty print a tuple
 */
template<class Ch, class Tr, class... Args>
auto operator<<(std::basic_ostream<Ch, Tr> &os, const std::tuple<Args...> &t) -> std::basic_ostream<Ch, Tr>&{
    os << "(";
    readdy::util::detail::print_tuple_impl(os, t, std::make_index_sequence<sizeof...(Args)>());
    return os << ")";
}
