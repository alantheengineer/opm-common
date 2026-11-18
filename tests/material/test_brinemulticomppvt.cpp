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

#define BOOST_TEST_MODULE BrineMulticompPvt
#include <boost/test/unit_test.hpp>

#include <opm/common/utility/SaltArray.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/material/components/BrineDynamic.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidsystems/blackoilpvt/BrineMulticompPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WaterPvtMultiplexer.hpp>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using Pvt = Opm::BrineMulticompPvt<double>;
using H2O = Opm::H2O<double>;
using Brine = Opm::BrineDynamic<double, H2O>;
using Salinity = Opm::SaltArray<double, Opm::SaltMassFraction>;
using Opm::SaltIndex;

constexpr double tRef = 288.71;
constexpr double pRef = 1.01325e5;
constexpr double tRes = 360.0;
constexpr double pRes = 3.0e7;

std::string
deckWith(const std::string& saltmc, const std::string& props)
{
    return R"(
RUNSPEC
DIMENS
1 1 1 /
WATER
BRINE
TABDIMS
1 2 /
SALTMC
)" + saltmc
        + R"( /
START
1 JAN 2026 /
GRID
DX
1 /
DY
1 /
DZ
1 /
TOPS
0 /
PORO
0.3 /
PERMX
100 /
PERMY
100 /
PERMZ
100 /
PROPS
)" + props
        + R"(
SCHEDULE
)";
}

// Charge-balanced brines: region 1 is a basin chloride brine, region 2 adds lithium and
// sulphate
const std::string salinitc = R"(
SALINITC
-- NA  K     CA   MG    CL     SO4   LI    SR    BR  HCO3 CO3 B
   2.0 0.05  0.4  0.05  2.97   0     0     0.01  0   0    0   0 /
   1.0 0.2   0    0.3   1.9    0.1   0.3   0     0   0    0   0 /
)";

struct Setup
{
    explicit Setup(const std::string& input)
        : deck(Opm::Parser {}.parseString(input))
        , eclState(deck)
        , schedule(deck, eclState, std::make_shared<Opm::Python>())
    {
    }

    Opm::Deck deck;
    Opm::EclipseState eclState;
    Opm::Schedule schedule;
};

template <class Multiplexer>
const Pvt&
brinePvt(Multiplexer& multiplexer, const Setup& setup)
{
    multiplexer.initFromState(setup.eclState, setup.schedule);
    BOOST_REQUIRE(multiplexer.approach() == Opm::WaterPvtApproach::BrineMulticomp);
    return multiplexer.template getRealPvt<Opm::WaterPvtApproach::BrineMulticomp>();
}

// Ion mass fractions of a region scaled to a new total salt mass fraction
Salinity
scaled(const Salinity& salinity, double total)
{
    Salinity result = salinity;
    const double factor = total / salinity.sum();
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[static_cast<SaltIndex>(i)] *= factor;
    }
    return result;
}

double
ratioToWater(double T, double p)
{
    const double p0 = std::max(1e5, H2O::vaporPressure(T));
    return H2O::liquidViscosity(T, p, true) / H2O::liquidViscosity(T, p0, true);
}

} // namespace

BOOST_AUTO_TEST_CASE(BasinPreset)
{
    const Setup setup(deckWith("", salinitc));
    Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
    const auto& pvt = brinePvt(multiplexer, setup);

    BOOST_CHECK(pvt.densityModel() == Pvt::DensityModel::MaoDuanYoung);
    BOOST_CHECK(pvt.viscosityModel() == Pvt::ViscosityModel::MaoDuan);
    BOOST_CHECK(pvt.enthalpyModel() == Pvt::EnthalpyModel::Driesner);
    BOOST_CHECK(pvt.pairing() == Pvt::PairingMethod::Proportional);
    BOOST_REQUIRE_EQUAL(pvt.numRegions(), 2U);

    const auto& config = setup.eclState.getCo2StoreConfig();
    for (unsigned r = 0; r < 2; ++r) {
        const auto& salinity = config.saltComponents(r);
        const double rhoRef = Brine::liquidDensityMaoDuan(tRef, pRef, salinity);
        BOOST_CHECK_CLOSE(pvt.waterReferenceDensity(r), rhoRef, 1e-10);

        // Total salt from BRINE: the SALINITC total and twice that
        for (const double factor : {1.0, 2.0}) {
            const double total = factor * salinity.sum();
            const double c = total * rhoRef;
            const auto brine = scaled(salinity, total);
            BOOST_CHECK_CLOSE(pvt.saltMassFraction(r, c), total, 1e-10);

            const double b = pvt.inverseFormationVolumeFactor(r, tRes, pRes, 0.0, c);
            BOOST_CHECK_CLOSE(b * rhoRef, Brine::liquidDensityMaoDuan(tRes, pRes, brine), 1e-8);

            const double mu = pvt.viscosity(r, tRes, pRes, 0.0, c);
            BOOST_CHECK_CLOSE(mu, Brine::liquidViscosityMaoDuan(tRes, pRes, brine), 1e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(SalarPreset)
{
    const Setup setup(deckWith("1* 1* SALAR", salinitc));
    Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
    const auto& pvt = brinePvt(multiplexer, setup);

    BOOST_CHECK(pvt.densityModel() == Pvt::DensityModel::LaliberteCooper);
    BOOST_CHECK(pvt.viscosityModel() == Pvt::ViscosityModel::Laliberte);
    BOOST_CHECK(pvt.enthalpyModel() == Pvt::EnthalpyModel::Laliberte);

    const auto& salinity = setup.eclState.getCo2StoreConfig().saltComponents(1);
    const double rhoRef = Brine::liquidDensityMulticompSalt(tRef, pRef, salinity, true, true);
    BOOST_CHECK_CLOSE(pvt.waterReferenceDensity(1), rhoRef, 1e-10);

    const double c = salinity.sum() * rhoRef;
    const double b = pvt.inverseFormationVolumeFactor(1, tRes, pRes, 0.0, c);
    BOOST_CHECK_CLOSE(
        b * rhoRef, Brine::liquidDensityMulticompSalt(tRes, pRes, salinity, true, true), 1e-8);

    // Laliberte viscosity takes pressure through pure water
    const double mu = pvt.viscosity(1, tRes, pRes, 0.0, c);
    BOOST_CHECK_CLOSE(
        mu, Brine::liquidViscosityMulticompSalt(tRes, salinity) * ratioToWater(tRes, pRes), 1e-8);
    BOOST_CHECK(mu > Brine::liquidViscosityMulticompSalt(tRes, salinity));
    BOOST_CHECK_CLOSE(pvt.viscosity(1, tRes, 1e5, 0.0, c),
                      Brine::liquidViscosityMulticompSalt(tRes, salinity),
                      1e-8);
}

BOOST_AUTO_TEST_CASE(LaliberteCooperWater)
{
    const Setup setup(deckWith("YES 1* SALAR", salinitc));
    Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
    const auto& pvt = brinePvt(multiplexer, setup);

    const auto& salinity = setup.eclState.getCo2StoreConfig().saltComponents(0);
    const double rhoRef = Brine::liquidDensityMulticompSalt(tRef, pRef, salinity, false);
    BOOST_CHECK_CLOSE(pvt.waterReferenceDensity(0), rhoRef, 1e-10);

    // That water density has no pressure term
    const double c = salinity.sum() * rhoRef;
    BOOST_CHECK_CLOSE(pvt.inverseFormationVolumeFactor(0, tRes, pRes, 0.0, c),
                      pvt.inverseFormationVolumeFactor(0, tRes, 1e5, 0.0, c),
                      1e-12);
}

BOOST_AUTO_TEST_CASE(ModelOverrides)
{
    const Setup setup(deckWith("1* 1* SALAR MAO_DUAN_YOUNG MAO_DUAN DRIESNER", salinitc));
    Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
    const auto& pvt = brinePvt(multiplexer, setup);

    BOOST_CHECK(pvt.densityModel() == Pvt::DensityModel::MaoDuanYoung);
    BOOST_CHECK(pvt.viscosityModel() == Pvt::ViscosityModel::MaoDuan);
    BOOST_CHECK(pvt.enthalpyModel() == Pvt::EnthalpyModel::Driesner);
}

BOOST_AUTO_TEST_CASE(FixedTotalSaltWithoutBrine)
{
    const Setup setup(deckWith("", salinitc));
    Opm::WaterPvtMultiplexer<double, false, false> multiplexer;
    const auto& pvt = brinePvt(multiplexer, setup);

    // Without brine in the model the SALINITC total is used whatever is passed in
    const auto& salinity = setup.eclState.getCo2StoreConfig().saltComponents(0);
    BOOST_CHECK_CLOSE(pvt.saltMassFraction(0, 0.0), salinity.sum(), 1e-12);
    const double b = pvt.inverseFormationVolumeFactor(0, tRes, pRes, 0.0, 0.0);
    BOOST_CHECK_CLOSE(
        b * pvt.waterReferenceDensity(0), Brine::liquidDensityMaoDuan(tRes, pRes, salinity), 1e-8);
}

BOOST_AUTO_TEST_CASE(SodiumChlorideWithoutSalinitc)
{
    const Setup setup(deckWith("", ""));
    Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
    const auto& pvt = brinePvt(multiplexer, setup);

    const double rhoW = H2O::liquidDensity(tRef, pRef, true);
    BOOST_CHECK_CLOSE(pvt.waterReferenceDensity(0), rhoW, 1e-8);

    // Salt from BRINE counts as NaCl
    const double c = 60.0;
    const double w = c / rhoW;
    const double mmNa = Opm::saltMolarMass<double>(SaltIndex::NA);
    const double mmCl = Opm::saltMolarMass<double>(SaltIndex::CL);
    Salinity nacl;
    nacl[SaltIndex::NA] = w * mmNa / (mmNa + mmCl);
    nacl[SaltIndex::CL] = w * mmCl / (mmNa + mmCl);
    BOOST_CHECK_CLOSE(pvt.inverseFormationVolumeFactor(0, tRes, pRes, 0.0, c) * rhoW,
                      Brine::liquidDensityMaoDuan(tRes, pRes, nacl),
                      1e-8);
}

BOOST_AUTO_TEST_CASE(SaltDerivatives)
{
    const Setup setup(deckWith("", salinitc));
    Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
    const auto& pvt = brinePvt(multiplexer, setup);

    using Evaluation = Opm::DenseAd::Evaluation<double, 1>;
    const Evaluation T(tRes);
    const Evaluation p(pRes);
    const Evaluation rsw(0.0);
    const double c0 = 120.0;
    const Evaluation c(c0, 0);

    const auto b = pvt.inverseFormationVolumeFactor(0, T, p, rsw, c);
    const auto mu = pvt.viscosity(0, T, p, rsw, c);
    BOOST_CHECK(b.derivative(0) > 0.0);
    BOOST_CHECK(mu.derivative(0) > 0.0);

    // Against a central difference
    const double dc = 1e-3;
    const double bUp = pvt.inverseFormationVolumeFactor(0, tRes, pRes, 0.0, c0 + dc);
    const double bDown = pvt.inverseFormationVolumeFactor(0, tRes, pRes, 0.0, c0 - dc);
    BOOST_CHECK_CLOSE(b.derivative(0), (bUp - bDown) / (2 * dc), 1e-4);
    const double muUp = pvt.viscosity(0, tRes, pRes, 0.0, c0 + dc);
    const double muDown = pvt.viscosity(0, tRes, pRes, 0.0, c0 - dc);
    BOOST_CHECK_CLOSE(mu.derivative(0), (muUp - muDown) / (2 * dc), 1e-4);
}

namespace
{

Pvt
directPvt(Pvt::DensityModel density, Pvt::ViscosityModel viscosity, Pvt::EnthalpyModel enthalpy)
{
    Opm::SaltArray<double, Opm::SaltMolality> molality;
    molality[SaltIndex::NA] = 1.5;
    molality[SaltIndex::CA] = 0.3;
    molality[SaltIndex::CL] = 2.1;
    Pvt pvt;
    pvt.setNumRegions(1);
    pvt.setSalinity(0, molality.convert_to<Opm::SaltMassFraction>());
    pvt.setModels(density, viscosity, enthalpy);
    pvt.setEnableSaltConcentration(true);
    pvt.initEnd();
    return pvt;
}

} // namespace

BOOST_AUTO_TEST_CASE(InternalEnergy)
{
    using D = Pvt::DensityModel;
    using V = Pvt::ViscosityModel;
    using E = Pvt::EnthalpyModel;

    for (const auto enthalpy : {E::Driesner, E::Laliberte}) {
        const auto pvt = directPvt(D::MaoDuanYoung, V::MaoDuan, enthalpy);
        const double c = 100.0;
        const double w = pvt.saltMassFraction(0, c);
        Opm::SaltArray<double, Opm::SaltMolality> molality;
        molality[SaltIndex::NA] = 1.5;
        molality[SaltIndex::CA] = 0.3;
        molality[SaltIndex::CL] = 2.1;
        const auto brine = scaled(molality.convert_to<Opm::SaltMassFraction>(), w);

        for (const double T : {320.0, 400.0}) {
            const double h = enthalpy == E::Driesner
                ? Brine::liquidEnthalpyLaliberteDriesner(T, pRes, brine)
                : Brine::liquidEnthalpyLaliberte(T, pRes, brine);
            const double rho = Brine::liquidDensityMaoDuan(T, pRes, brine);
            BOOST_CHECK_CLOSE(pvt.internalEnergy(0, T, pRes, 0.0, c), h - pRes / rho, 1e-8);
        }
    }

    // TABLE takes the curve as it is
    auto pvt = directPvt(D::MaoDuanYoung, V::MaoDuan, E::Table);
    BOOST_CHECK_THROW(pvt.internalEnergy(0, 350.0, pRes, 0.0, 100.0), std::runtime_error);
    Pvt::TabulatedOneDFunction curve;
    curve.setXYContainers(std::vector<double> {300.0, 400.0}, std::vector<double> {1.2e6, 1.6e6});
    pvt.setInternalEnergyCurve(0, curve);
    BOOST_CHECK_CLOSE(pvt.internalEnergy(0, 350.0, pRes, 0.0, 100.0), 1.4e6, 1e-10);
}

BOOST_AUTO_TEST_CASE(CalibrationFactors)
{
    auto pvt = directPvt(Pvt::DensityModel::MaoDuanYoung,
                         Pvt::ViscosityModel::MaoDuan,
                         Pvt::EnthalpyModel::Driesner);
    const double rhoRef = pvt.waterReferenceDensity(0);
    const double b = pvt.inverseFormationVolumeFactor(0, tRes, pRes, 0.0, 0.0);
    const double mu = pvt.viscosity(0, tRes, pRes, 0.0, 0.0);

    pvt.setCalibrationFactors(0, 1.02, 1.1);

    // The density factor moves the surface density and leaves b alone
    BOOST_CHECK_CLOSE(pvt.waterReferenceDensity(0), 1.02 * rhoRef, 1e-12);
    BOOST_CHECK_CLOSE(pvt.inverseFormationVolumeFactor(0, tRes, pRes, 0.0, 0.0), b, 1e-12);
    BOOST_CHECK_CLOSE(pvt.viscosity(0, tRes, pRes, 0.0, 0.0), 1.1 * mu, 1e-12);

    BOOST_CHECK_THROW(pvt.setCalibrationFactors(0, 0.0, 1.0), std::invalid_argument);
    BOOST_CHECK_THROW(pvt.setCalibrationFactors(0, 1.0, -1.0), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(MultiplexerCopy)
{
    const Setup setup(deckWith("", salinitc));
    Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
    const auto& pvt = brinePvt(multiplexer, setup);

    const Opm::WaterPvtMultiplexer<double, true, true> copy(multiplexer);
    BOOST_REQUIRE(copy.approach() == Opm::WaterPvtApproach::BrineMulticomp);
    BOOST_CHECK(copy.getRealPvt<Opm::WaterPvtApproach::BrineMulticomp>() == pvt);
    BOOST_CHECK_CLOSE(copy.waterReferenceDensity(1), multiplexer.waterReferenceDensity(1), 1e-12);
}

namespace
{

// The parts of a fluid state the brine PVT reads
template <class Evaluation>
struct CellFluidState
{
    using ValueType = Evaluation;
    static constexpr unsigned waterPhaseIdx = 0;

    const Evaluation& temperature(unsigned) const
    {
        return T;
    }
    const Evaluation& pressure(unsigned) const
    {
        return p;
    }
    const Evaluation& saltConcentration() const
    {
        return salt;
    }
    const Opm::ElectrolyteFractions* brineComposition() const
    {
        return composition;
    }

    Evaluation T;
    Evaluation p;
    Evaluation salt;
    const Opm::ElectrolyteFractions* composition = nullptr;
};

} // namespace

BOOST_AUTO_TEST_CASE(CellCompositionMatchesRegion)
{
    const Setup setup(deckWith("", salinitc));
    Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
    const auto& pvt = brinePvt(multiplexer, setup);

    // A cell with the region's own ratios gives the region's properties
    const auto cell = pvt.cellFractions(pvt.salinity(1));
    const double c = 150.0;
    CellFluidState<double> fs {tRes, pRes, c, &cell};
    const auto [b, mu] = pvt.inverseFormationVolumeFactorAndViscosity(fs, 1);
    BOOST_CHECK_CLOSE(b, pvt.inverseFormationVolumeFactor(1, tRes, pRes, 0.0, c), 1e-5);
    BOOST_CHECK_CLOSE(mu, pvt.viscosity(1, tRes, pRes, 0.0, c), 1e-5);
    BOOST_CHECK_CLOSE(pvt.internalEnergy(1, tRes, pRes, 0.0, c, &cell),
                      pvt.internalEnergy(1, tRes, pRes, 0.0, c),
                      1e-5);

    // No fractions means the region's
    fs.composition = nullptr;
    const auto [bRegion, muRegion] = pvt.inverseFormationVolumeFactorAndViscosity(fs, 1);
    BOOST_CHECK_EQUAL(bRegion, pvt.inverseFormationVolumeFactor(1, tRes, pRes, 0.0, c));
    BOOST_CHECK_EQUAL(muRegion, pvt.viscosity(1, tRes, pRes, 0.0, c));
}

BOOST_AUTO_TEST_CASE(CellCompositionOwnRatios)
{
    const Setup setup(deckWith("", salinitc));
    Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
    const auto& pvt = brinePvt(multiplexer, setup);

    // Region 1 deck, but this cell holds a lithium-rich chloride brine
    Opm::SaltArray<double, Opm::SaltMolality> molality;
    molality[SaltIndex::NA] = 1.0;
    molality[SaltIndex::LI] = 1.5;
    molality[SaltIndex::CA] = 0.2;
    molality[SaltIndex::CL] = 2.9;
    const auto ions = molality.convert_to<Opm::SaltMassFraction>();
    const auto cell = pvt.cellFractions(ions);

    const double c = 120.0;
    const double w = pvt.saltMassFraction(0, c);
    const auto brine = scaled(ions, w);
    CellFluidState<double> fs {tRes, pRes, c, &cell};
    const auto [b, mu] = pvt.inverseFormationVolumeFactorAndViscosity(fs, 0);

    // Float fractions limit the match to about 1e-7
    const double rhoRef = pvt.waterReferenceDensity(0);
    BOOST_CHECK_CLOSE(b * rhoRef, Brine::liquidDensityMaoDuan(tRes, pRes, brine), 1e-4);
    BOOST_CHECK_CLOSE(mu, Brine::liquidViscosityMaoDuan(tRes, pRes, brine), 1e-4);

    // Lithium chloride is lighter than the region's brine at the same salt
    BOOST_CHECK(b < pvt.inverseFormationVolumeFactor(0, tRes, pRes, 0.0, c));

    // The multiplexer passes the fractions on for internal energy
    const double h = Brine::liquidEnthalpyLaliberteDriesner(tRes, pRes, brine);
    const double rho = Brine::liquidDensityMaoDuan(tRes, pRes, brine);
    BOOST_CHECK_CLOSE(
        multiplexer.internalEnergy(0u, tRes, pRes, 0.0, c, &cell), h - pRes / rho, 1e-4);
}

BOOST_AUTO_TEST_CASE(CellFractionsSumToPaired)
{
    const Setup setup(deckWith("", salinitc));
    Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
    const auto& pvt = brinePvt(multiplexer, setup);

    // Charge-balanced ions pair completely, so the fractions sum to one
    const auto cell = pvt.cellFractions(pvt.salinity(0));
    double sum = 0.0;
    for (const float f : cell) {
        BOOST_CHECK(f >= 0.0f);
        sum += f;
    }
    BOOST_CHECK_CLOSE(sum, 1.0, 1e-4);

    // No salt gives no fractions
    const auto empty = pvt.cellFractions(Salinity {});
    for (const float f : empty) {
        BOOST_CHECK_EQUAL(f, 0.0f);
    }
}
