/*
  Copyright 2026 Dr Alan Tominey.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef OPM_FIXED_CAPACITY_LIST_HPP
#define OPM_FIXED_CAPACITY_LIST_HPP

#include <array>
#include <cstddef>
#include <stdexcept>

namespace Opm
{

/*!
 * \brief A list of at most N elements stored in place, so filling it never allocates.
 *
 * \tparam T Element type, default constructible
 * \tparam N Capacity
 */
template <class T, std::size_t N>
class FixedCapacityList
{
public:
    void push_back(const T& value)
    {
        if (size_ == N) {
            throw std::length_error("FixedCapacityList is full");
        }
        data_[size_++] = value;
    }

    [[nodiscard]] std::size_t size() const noexcept
    {
        return size_;
    }

    [[nodiscard]] bool empty() const noexcept
    {
        return size_ == 0;
    }

    [[nodiscard]] static constexpr std::size_t capacity() noexcept
    {
        return N;
    }

    const T& operator[](std::size_t i) const
    {
        return data_[i];
    }

    T& operator[](std::size_t i)
    {
        return data_[i];
    }

    auto begin() noexcept
    {
        return data_.begin();
    }

    auto end() noexcept
    {
        return data_.begin() + static_cast<std::ptrdiff_t>(size_);
    }

    [[nodiscard]] auto begin() const noexcept
    {
        return data_.begin();
    }

    [[nodiscard]] auto end() const noexcept
    {
        return data_.begin() + static_cast<std::ptrdiff_t>(size_);
    }

private:
    std::array<T, N> data_ {};
    std::size_t size_ = 0;
};

} // namespace Opm

#endif // OPM_FIXED_CAPACITY_LIST_HPP
