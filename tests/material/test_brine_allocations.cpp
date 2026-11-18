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
#include "config.h"

#define BOOST_TEST_MODULE BrineAllocationTests
#include <boost/test/unit_test.hpp>

#include <opm/common/utility/SaltArray.hpp>
#include <opm/common/utility/SaltElectrolytes.hpp>
#include <opm/material/components/BrineDynamic.hpp>
#include <opm/material/components/SimpleHuDuanH2O.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <cstddef>
#include <cstdlib>
#include <new>

// Brine properties are evaluated for every cell in every Newton iteration, so they must
// not allocate. This binary replaces the global allocator to count allocations.
namespace
{
bool counting = false;
std::size_t allocations = 0;
} // namespace

void*
operator new(std::size_t size)
{
    if (counting) {
        ++allocations;
    }
    if (void* p = std::malloc(size == 0 ? 1 : size)) {
        return p;
    }
    throw std::bad_alloc();
}

void
operator delete(void* p) noexcept
{
    std::free(p);
}

void
operator delete(void* p, std::size_t) noexcept
{
    std::free(p);
}

namespace
{
using Scalar = double;
using Evaluation = Opm::DenseAd::Evaluation<Scalar, 3>;
using Brine = Opm::BrineDynamic<Scalar, Opm::SimpleHuDuanH2O<Scalar>>;

Opm::SaltArray<Evaluation, Opm::SaltMassFraction>
mixedBrine()
{
    Opm::SaltArray<Evaluation, Opm::SaltMolality> molality;
    molality[Opm::SaltIndex::NA] = Evaluation(2.0, 0);
    molality[Opm::SaltIndex::K] = 0.1;
    molality[Opm::SaltIndex::CA] = 0.4;
    molality[Opm::SaltIndex::MG] = 0.1;
    molality[Opm::SaltIndex::CL] = 3.0;
    molality[Opm::SaltIndex::SO4] = 0.05;
    return molality.convert_to<Opm::SaltMassFraction>();
}
} // namespace

BOOST_AUTO_TEST_CASE(BrinePropertiesDoNotAllocate)
{
    const auto salinity = mixedBrine();
    const Evaluation T(353.15, 1);
    const Evaluation p(2.0e7, 2);

    // Once first, so one-time initialisation is not counted.
    auto rho = Brine::liquidDensityMulticompSalt(T, p, salinity);
    auto mu = Brine::liquidViscosityMulticompSalt(T, salinity);
    auto cp = Brine::liquidHeatCapacityLaliberte(T, salinity);
    auto h = Brine::liquidEnthalpyLaliberte(T, p, salinity);

    counting = true;
    for (int i = 0; i < 10; ++i) {
        rho += Brine::liquidDensityMulticompSalt(T, p, salinity);
        mu += Brine::liquidViscosityMulticompSalt(T, salinity);
        cp += Brine::liquidHeatCapacityLaliberte(T, salinity);
        h += Brine::liquidEnthalpyLaliberte(T, p, salinity);
    }
    const auto electrolytes = Opm::SaltElectrolytes<Scalar>::fromSaltComponents(salinity);
    counting = false;

    BOOST_CHECK_EQUAL(allocations, 0U);
    BOOST_CHECK(rho.value() > 0.0);
    BOOST_CHECK(mu.value() > 0.0);
    BOOST_CHECK(cp.value() > 0.0);
    BOOST_CHECK(h.value() > 0.0);
    BOOST_CHECK(!electrolytes.empty());
}
