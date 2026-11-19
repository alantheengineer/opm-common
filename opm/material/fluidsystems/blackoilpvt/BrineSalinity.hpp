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
#ifndef OPM_BRINE_SALINITY_HPP
#define OPM_BRINE_SALINITY_HPP

#include <opm/common/utility/SaltArray.hpp>
#include <opm/common/utility/gpuDecorators.hpp>

#include <cstddef>

namespace Opm::BrineSalinity
{

// The brine module stores salt per surface volume of brine, and b_w is taken against the
// region's surface brine density, so salt mass over brine mass is the concentration over
// that density at any pressure and temperature.

/// Total salt mass fraction from the BRINE salt concentration [kg/sm3]
template <class LhsEval, class Scalar>
OPM_HOST_DEVICE LhsEval
massFraction(const LhsEval& saltConcentration, Scalar surfaceDensity)
{
    return saltConcentration / surfaceDensity;
}

/// Salt mass fractions per ion, in the ratios of the region's composition, or as NaCl
/// when the composition is empty
template <class LhsEval, class Scalar>
OPM_HOST_DEVICE SaltArray<LhsEval, SaltMassFraction>
massFractions(const LhsEval& saltConcentration,
              Scalar surfaceDensity,
              const SaltArray<Scalar, SaltMassFraction>& composition)
{
    const LhsEval total = massFraction(saltConcentration, surfaceDensity);
    SaltArray<LhsEval, SaltMassFraction> salinity;
    const Scalar compositionTotal = composition.sum();
    if (compositionTotal > 0.0) {
        for (std::size_t i = 0; i < SaltArray<Scalar, SaltMassFraction>::ncomp; ++i) {
            const auto ion = static_cast<SaltIndex>(i);
            salinity[ion] = total * (composition[ion] / compositionTotal);
        }
    } else {
        const Scalar na = saltMolarMass<Scalar>(SaltIndex::NA);
        const Scalar cl = saltMolarMass<Scalar>(SaltIndex::CL);
        salinity[SaltIndex::NA] = total * (na / (na + cl));
        salinity[SaltIndex::CL] = total * (cl / (na + cl));
    }
    return salinity;
}

} // namespace Opm::BrineSalinity

#endif // OPM_BRINE_SALINITY_HPP
