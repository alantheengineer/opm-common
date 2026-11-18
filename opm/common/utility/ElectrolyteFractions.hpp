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
#ifndef OPM_ELECTROLYTE_FRACTIONS_HPP
#define OPM_ELECTROLYTE_FRACTIONS_HPP

#include <array>
#include <cstddef>

namespace Opm
{

/// Number of electrolytes in SaltElectrolytes
inline constexpr std::size_t numBrineElectrolytes = 17;

/// Mass fraction of each electrolyte per unit salt mass fraction, in SaltElectrolytes
/// order. Float is enough for ion ratios and halves the per-cell storage.
using ElectrolyteFractions = std::array<float, numBrineElectrolytes>;

} // namespace Opm

#endif
