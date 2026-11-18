/*
  Copyright 2026 NORCE.

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
#include <config.h>

#include <opm/common/utility/SaltElectrolytes.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/components/CaIon.hpp>
#include <opm/material/components/ClIon.hpp>
#include <opm/material/components/KIon.hpp>
#include <opm/material/components/MgIon.hpp>
#include <opm/material/components/NaIon.hpp>
#include <opm/material/components/SO4Ion.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <type_traits>
#include <vector>

#define BOOST_TEST_MODULE SaltElectrolytesTests
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

using namespace Opm;

using EvalTypes = boost::mpl::list<float,
                                   double,
                                   DenseAd::Evaluation<float, 3>,
                                   DenseAd::Evaluation<double, 3> >;

BOOST_AUTO_TEST_CASE_TEMPLATE(MolarMassTests, Eval, EvalTypes)
{
    using Scalar = typename MathToolbox<Eval>::Scalar;
    using Electrolyte = typename SaltElectrolytes<Scalar>::Electrolyte;

    const auto tol = static_cast<Scalar>(1.0e-4);
    BOOST_CHECK_CLOSE(SaltElectrolytes<Scalar>::molarMass(Electrolyte::NaCl),
                       NaIon<Scalar>::molarMass() + ClIon<Scalar>::molarMass(),
                       tol);
    BOOST_CHECK_CLOSE(SaltElectrolytes<Scalar>::molarMass(Electrolyte::KCl),
                       KIon<Scalar>::molarMass() + ClIon<Scalar>::molarMass(),
                       tol);
    BOOST_CHECK_CLOSE(SaltElectrolytes<Scalar>::molarMass(Electrolyte::CaCl2),
                       CaIon<Scalar>::molarMass() + 2.0 * ClIon<Scalar>::molarMass(),
                       tol);
    BOOST_CHECK_CLOSE(SaltElectrolytes<Scalar>::molarMass(Electrolyte::MgCl2),
                       MgIon<Scalar>::molarMass() + 2.0 * ClIon<Scalar>::molarMass(),
                       tol);
    BOOST_CHECK_CLOSE(SaltElectrolytes<Scalar>::molarMass(Electrolyte::MgSO4),
                       MgIon<Scalar>::molarMass() + SO4Ion<Scalar>::molarMass(),
                       tol);
    BOOST_CHECK_CLOSE(SaltElectrolytes<Scalar>::molarMass(Electrolyte::K2SO4),
                       2.0 * KIon<Scalar>::molarMass() + SO4Ion<Scalar>::molarMass(),
                       tol);
    BOOST_CHECK_CLOSE(SaltElectrolytes<Scalar>::molarMass(Electrolyte::Na2SO4),
                       2.0 * NaIon<Scalar>::molarMass() + SO4Ion<Scalar>::molarMass(),
                       tol);
    BOOST_CHECK_CLOSE(SaltElectrolytes<Scalar>::molarMass(Electrolyte::LiCl),
                      static_cast<Scalar>(42.39e-3),
                      static_cast<Scalar>(0.01));
    BOOST_CHECK_CLOSE(SaltElectrolytes<Scalar>::molarMass(Electrolyte::SrCl2),
                      static_cast<Scalar>(158.53e-3),
                      static_cast<Scalar>(0.01));
    BOOST_CHECK_CLOSE(SaltElectrolytes<Scalar>::molarMass(Electrolyte::K2CO3),
                      static_cast<Scalar>(138.21e-3),
                      static_cast<Scalar>(0.01));
    BOOST_CHECK_CLOSE(SaltElectrolytes<Scalar>::molarMass(Electrolyte::BOH3),
                      static_cast<Scalar>(61.83e-3),
                      static_cast<Scalar>(0.01));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ZeroSalinityTests, Eval, EvalTypes)
{
    using Scalar = typename MathToolbox<Eval>::Scalar;

    // Default-constructed SaltArray is all zero: no electrolyte can be formed
    SaltArray<Eval, SaltMassFraction> salinity;
    auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity);
    BOOST_CHECK(electrolytes.empty());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(SingleElectrolyteTests, Eval, EvalTypes)
{
    using Scalar = typename MathToolbox<Eval>::Scalar;
    using Electrolyte = typename SaltElectrolytes<Scalar>::Electrolyte;

    // Stoichiometrically balanced cation/anion molalities (no leftover after pairing) for
    // each supported electrolyte
    struct Case
    {
        SaltIndex cationIdx;
        Scalar cationMolal;
        SaltIndex anionIdx;
        Scalar anionMolal;
        Electrolyte electrolyte;
    };
    // clang-format off
    const std::vector<Case> cases = {
        {SaltIndex::NA, 0.7, SaltIndex::CL, 0.7, Electrolyte::NaCl},
        {SaltIndex::K,  0.7, SaltIndex::CL, 0.7, Electrolyte::KCl},
        {SaltIndex::CA, 0.7, SaltIndex::CL, 1.4, Electrolyte::CaCl2},
        {SaltIndex::MG, 0.7, SaltIndex::CL, 1.4, Electrolyte::MgCl2},
        {SaltIndex::MG, 0.7, SaltIndex::SO4, 0.7, Electrolyte::MgSO4},
        {SaltIndex::NA, 1.4, SaltIndex::SO4, 0.7, Electrolyte::Na2SO4},
        {SaltIndex::K,  1.4, SaltIndex::SO4, 0.7, Electrolyte::K2SO4},
        {SaltIndex::LI, 0.7, SaltIndex::CL, 0.7, Electrolyte::LiCl},
        {SaltIndex::LI, 1.4, SaltIndex::SO4, 0.7, Electrolyte::Li2SO4},
        {SaltIndex::SR, 0.7, SaltIndex::CL, 1.4, Electrolyte::SrCl2},
        {SaltIndex::NA, 0.7, SaltIndex::BR, 0.7, Electrolyte::NaBr},
        {SaltIndex::K,  0.7, SaltIndex::BR, 0.7, Electrolyte::KBr},
        {SaltIndex::NA, 0.7, SaltIndex::HCO3, 0.7, Electrolyte::NaHCO3},
        {SaltIndex::NA, 1.4, SaltIndex::CO3, 0.7, Electrolyte::Na2CO3},
        {SaltIndex::K,  0.7, SaltIndex::HCO3, 0.7, Electrolyte::KHCO3},
        {SaltIndex::K,  1.4, SaltIndex::CO3, 0.7, Electrolyte::K2CO3},
        {SaltIndex::B,  0.7, SaltIndex::B, 0.0, Electrolyte::BOH3},
    };
    // clang-format on
    BOOST_CHECK_EQUAL(cases.size(), SaltElectrolytes<Scalar>::numElectrolytes);

    for (const auto& c : cases) {
        SaltArray<Eval, SaltMolality> molal;
        molal[c.cationIdx] += c.cationMolal;
        molal[c.anionIdx] += c.anionMolal;
        SaltArray<Eval, SaltMassFraction> salinity = molal.template convert_to<SaltMassFraction>();

        auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity);
        BOOST_REQUIRE_EQUAL(electrolytes.size(), 1U);
        BOOST_CHECK(electrolytes[0].first == c.electrolyte);

        // No residual molality: the single electrolyte's mass fraction reconstructs the
        // full dissolved-salt mass fraction
        BOOST_CHECK_CLOSE(electrolytes[0].second, salinity.sum(), static_cast<Scalar>(1.0e-4));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(GreedyMultipleElectrolytesOrderingTests, Eval, EvalTypes)
{
    using Scalar = typename MathToolbox<Eval>::Scalar;
    constexpr auto greedy = SaltElectrolytes<Scalar>::PairingMethod::Greedy;
    using Electrolyte = typename SaltElectrolytes<Scalar>::Electrolyte;

    // Mg2+ pairs fully with SO4^2- into MgSO4, and K+ pairs fully with Cl- into KCl: two
    // electrolytes, no leftover molality
    SaltArray<Eval, SaltMolality> molal;
    molal[SaltIndex::MG] = 1.0;
    molal[SaltIndex::SO4] = 1.0;
    molal[SaltIndex::K] = 2.0;
    molal[SaltIndex::CL] = 2.0;
    SaltArray<Eval, SaltMassFraction> salinity = molal.template convert_to<SaltMassFraction>();

    auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity, greedy);
    BOOST_REQUIRE_EQUAL(electrolytes.size(), 2U);

    // Electrolytes are generated anion-first in decreasing ion strength: SO4 (and its
    // MgSO4 pairing) before Cl (and its KCl pairing)
    BOOST_CHECK(electrolytes[0].first == Electrolyte::MgSO4);
    BOOST_CHECK(electrolytes[1].first == Electrolyte::KCl);

    // No residual molality: the two electrolyte mass fractions sum to the total
    // dissolved-salt mass fraction
    BOOST_CHECK_CLOSE(electrolytes[0].second + electrolytes[1].second,
                       salinity.sum(),
                       static_cast<Scalar>(1.0e-4));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(UnpairedIonsAreIgnoredTests, Eval, EvalTypes)
{
    using Scalar = typename MathToolbox<Eval>::Scalar;

    // Cations only, no anion present to pair with: no electrolyte can be formed and the
    // dissolved salt is ignored (with a log warning), not thrown as an error
    SaltArray<Eval, SaltMolality> molal;
    molal[SaltIndex::NA] = 1.0;
    molal[SaltIndex::CA] = 2.0;
    SaltArray<Eval, SaltMassFraction> salinity = molal.template convert_to<SaltMassFraction>();
    BOOST_REQUIRE(salinity.any_nonzero());

    auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity);
    BOOST_CHECK(electrolytes.empty());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(PartialRemainderIsIgnoredTests, Eval, EvalTypes)
{
    using Scalar = typename MathToolbox<Eval>::Scalar;
    using Electrolyte = typename SaltElectrolytes<Scalar>::Electrolyte;

    // Excess Na+ beyond what the available Cl- can pair with; only the paired part forms
    // NaCl, the leftover Na+ molality is ignored
    SaltArray<Eval, SaltMolality> molal;
    molal[SaltIndex::NA] = 3.0;
    molal[SaltIndex::CL] = 1.0;
    SaltArray<Eval, SaltMassFraction> salinity = molal.template convert_to<SaltMassFraction>();

    auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity);
    BOOST_REQUIRE_EQUAL(electrolytes.size(), 1U);
    BOOST_CHECK(electrolytes[0].first == Electrolyte::NaCl);

    // Only 1 of the 3 molal Na+ could pair with Cl-, so the returned NaCl mass fraction is
    // strictly less than the total dissolved-salt mass fraction
    BOOST_CHECK(electrolytes[0].second < salinity.sum());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(GreedyIncompatiblePairIsSkippedAndRetriedTests, Eval, EvalTypes)
{
    using Scalar = typename MathToolbox<Eval>::Scalar;
    constexpr auto greedy = SaltElectrolytes<Scalar>::PairingMethod::Greedy;
    using Electrolyte = typename SaltElectrolytes<Scalar>::Electrolyte;

    // CaCl2 + K2SO4 mixture: since CaSO4 has not been implemented (yet) but Ca++ and SO4-- are
    // used first to build an electrolyte, they should be skipped. The algorithm should move on to
    // K+ and SO4--, and then to Ca++ and Cl-
    SaltArray<Eval, SaltMolality> molal;
    molal[SaltIndex::CA] = 1.0;
    molal[SaltIndex::CL] = 2.0;
    molal[SaltIndex::K] = 2.0;
    molal[SaltIndex::SO4] = 1.0;
    SaltArray<Eval, SaltMassFraction> salinity = molal.template convert_to<SaltMassFraction>();

    auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity, greedy);
    BOOST_REQUIRE_EQUAL(electrolytes.size(), 2U);
    BOOST_CHECK(electrolytes[0].first == Electrolyte::K2SO4);
    BOOST_CHECK(electrolytes[1].first == Electrolyte::CaCl2);

    // No residual molality: the two electrolyte mass fractions sum to the total
    // dissolved-salt mass fraction
    BOOST_CHECK_CLOSE(electrolytes[0].second + electrolytes[1].second,
                       salinity.sum(),
                       static_cast<Scalar>(1.0e-4));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(GreedyBasinBrineTests, Eval, EvalTypes)
{
    using Scalar = typename MathToolbox<Eval>::Scalar;
    constexpr auto greedy = SaltElectrolytes<Scalar>::PairingMethod::Greedy;
    using Electrolyte = typename SaltElectrolytes<Scalar>::Electrolyte;

    // A charge-balanced Na-Ca-Cl brine with every ion [mol/kg water]. Anions pair in order
    // SO4, CO3, Cl, Br, HCO3, each with cations in order Ca, Mg, Sr, Na, K, Li.
    SaltArray<Eval, SaltMolality> molal;
    molal[SaltIndex::NA] = 3.0;
    molal[SaltIndex::K] = 0.1;
    molal[SaltIndex::LI] = 0.03;
    molal[SaltIndex::CA] = 0.8;
    molal[SaltIndex::MG] = 0.1;
    molal[SaltIndex::SR] = 0.02;
    molal[SaltIndex::SO4] = 0.01;
    molal[SaltIndex::CO3] = 0.001;
    molal[SaltIndex::CL] = 4.913;
    molal[SaltIndex::BR] = 0.03;
    molal[SaltIndex::HCO3] = 0.005;
    molal[SaltIndex::B] = 0.01;
    const auto salinity = molal.template convert_to<SaltMassFraction>();

    const auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity, greedy);
    // clang-format off
    const std::vector<Electrolyte> expected = {
        Electrolyte::MgSO4, Electrolyte::Na2CO3, Electrolyte::CaCl2, Electrolyte::MgCl2,
        Electrolyte::SrCl2, Electrolyte::NaCl, Electrolyte::KCl, Electrolyte::KBr,
        Electrolyte::BOH3};
    // clang-format on
    BOOST_REQUIRE_EQUAL(electrolytes.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        BOOST_CHECK(electrolytes[i].first == expected[i]);
    }

    // Cl runs out before Li, which has no Br or HCO3 pairing. Boron needs no partner. Left
    // unpaired: 0.03 Li, 0.025 Br and 0.005 HCO3.
    const auto unpaired = SaltElectrolytes<Scalar>::unpairedMolality(salinity, greedy);
    BOOST_CHECK_CLOSE(scalarValue(unpaired), static_cast<Scalar>(0.06), static_cast<Scalar>(0.1));
}

// Molality of one electrolyte in a list of electrolyte mass fractions
template <class Scalar, class List, class Electrolyte>
Scalar
molalityOf(const List& electrolytes, const Electrolyte which)
{
    Scalar sumW = 0.0;
    for (const auto& e : electrolytes) {
        sumW += scalarValue(e.second);
    }
    for (const auto& e : electrolytes) {
        if (e.first == which) {
            return scalarValue(e.second) / SaltElectrolytes<Scalar>::molarMass(which) / (1 - sumW);
        }
    }
    return 0.0;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ProportionalPairingTests, Eval, EvalTypes)
{
    using Scalar = typename MathToolbox<Eval>::Scalar;
    using Electrolyte = typename SaltElectrolytes<Scalar>::Electrolyte;

    // Equal charge on Mg and K, and on SO4 and Cl: each pair takes a quarter of the 4
    // equivalents [mol/kg water]
    SaltArray<Eval, SaltMolality> molal;
    molal[SaltIndex::MG] = 1.0;
    molal[SaltIndex::SO4] = 1.0;
    molal[SaltIndex::K] = 2.0;
    molal[SaltIndex::CL] = 2.0;
    const auto salinity = molal.template convert_to<SaltMassFraction>();
    const auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity);
    BOOST_REQUIRE_EQUAL(electrolytes.size(), 4U);
    const auto molalOf = [&](Electrolyte e) { return molalityOf<Scalar>(electrolytes, e); };
    const Scalar tol = 1e-3;
    BOOST_CHECK_CLOSE(molalOf(Electrolyte::MgSO4), Scalar(0.5), tol);
    BOOST_CHECK_CLOSE(molalOf(Electrolyte::MgCl2), Scalar(0.5), tol);
    BOOST_CHECK_CLOSE(molalOf(Electrolyte::K2SO4), Scalar(0.5), tol);
    BOOST_CHECK_CLOSE(molalOf(Electrolyte::KCl), Scalar(1.0), tol);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ProportionalPairingRefitsMissingPairsTests, Eval, EvalTypes)
{
    using Scalar = typename MathToolbox<Eval>::Scalar;
    using Electrolyte = typename SaltElectrolytes<Scalar>::Electrolyte;

    // The NaCl, MgSO4 and NaBr mixture of Lilek and Zuend (2022), supplement S1. There is
    // no MgBr2, so Br must all go to Na and the other pairs make room for it.
    SaltArray<Eval, SaltMolality> molal;
    molal[SaltIndex::NA] = 5.0;
    molal[SaltIndex::CL] = 2.0;
    molal[SaltIndex::MG] = 9.0;
    molal[SaltIndex::SO4] = 9.0;
    molal[SaltIndex::BR] = 3.0;
    const auto salinity = molal.template convert_to<SaltMassFraction>();

    const auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity);
    BOOST_CHECK_EQUAL(electrolytes.size(), 5U);

    // Na keeps 2 equivalents for SO4 and Cl, split 9 to 1 as for Mg. The fit stops near
    // machine precision, so float gets a wider tolerance [%].
    const auto molalOf = [&](Electrolyte e) { return molalityOf<Scalar>(electrolytes, e); };
    const Scalar tol = std::is_same_v<Scalar, float> ? 1e-2 : 1e-4;
    BOOST_CHECK_CLOSE(molalOf(Electrolyte::NaBr), Scalar(3.0), tol);
    BOOST_CHECK_CLOSE(molalOf(Electrolyte::NaCl), Scalar(0.2), tol);
    BOOST_CHECK_CLOSE(molalOf(Electrolyte::Na2SO4), Scalar(0.9), tol);
    BOOST_CHECK_CLOSE(molalOf(Electrolyte::MgSO4), Scalar(8.1), tol);
    BOOST_CHECK_CLOSE(molalOf(Electrolyte::MgCl2), Scalar(0.9), tol);

    const auto unpaired = SaltElectrolytes<Scalar>::unpairedMolality(salinity);
    BOOST_CHECK_SMALL(scalarValue(unpaired), Scalar(1e-4));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ProportionalBasinBrineTests, Eval, EvalTypes)
{
    using Scalar = typename MathToolbox<Eval>::Scalar;

    // The charge-balanced brine of the greedy test. Every ion has a supported partner, so
    // proportional pairing uses all of them, Li and Br included.
    SaltArray<Eval, SaltMolality> molal;
    molal[SaltIndex::NA] = 3.0;
    molal[SaltIndex::K] = 0.1;
    molal[SaltIndex::LI] = 0.03;
    molal[SaltIndex::CA] = 0.8;
    molal[SaltIndex::MG] = 0.1;
    molal[SaltIndex::SR] = 0.02;
    molal[SaltIndex::SO4] = 0.01;
    molal[SaltIndex::CO3] = 0.001;
    molal[SaltIndex::CL] = 4.913;
    molal[SaltIndex::BR] = 0.03;
    molal[SaltIndex::HCO3] = 0.005;
    molal[SaltIndex::B] = 0.01;
    const auto salinity = molal.template convert_to<SaltMassFraction>();

    const auto unpaired = SaltElectrolytes<Scalar>::unpairedMolality(salinity);
    BOOST_CHECK_SMALL(scalarValue(unpaired), Scalar(1e-4));

    // Paired mass makes up all the dissolved mass
    const auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity);
    Scalar sumW = 0.0;
    for (const auto& e : electrolytes) {
        sumW += scalarValue(e.second);
    }
    BOOST_CHECK_CLOSE(sumW, scalarValue(salinity.sum()), Scalar(1e-3));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ProportionalUnbalancedTests, Eval, EvalTypes)
{
    using Scalar = typename MathToolbox<Eval>::Scalar;
    using Electrolyte = typename SaltElectrolytes<Scalar>::Electrolyte;

    // 3 Na against 1 Cl: 1 NaCl, 2 Na unpaired. Ca and SO4 alone have no supported pair.
    SaltArray<Eval, SaltMolality> molal;
    molal[SaltIndex::NA] = 3.0;
    molal[SaltIndex::CL] = 1.0;
    auto salinity = molal.template convert_to<SaltMassFraction>();
    auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity);
    BOOST_REQUIRE_EQUAL(electrolytes.size(), 1U);
    BOOST_CHECK(electrolytes[0].first == Electrolyte::NaCl);
    BOOST_CHECK_CLOSE(scalarValue(SaltElectrolytes<Scalar>::unpairedMolality(salinity)),
                      Scalar(2.0),
                      Scalar(1e-3));

    SaltArray<Eval, SaltMolality> gypsum;
    gypsum[SaltIndex::CA] = 0.01;
    gypsum[SaltIndex::SO4] = 0.01;
    salinity = gypsum.template convert_to<SaltMassFraction>();
    BOOST_CHECK(SaltElectrolytes<Scalar>::fromSaltComponents(salinity).empty());
}
