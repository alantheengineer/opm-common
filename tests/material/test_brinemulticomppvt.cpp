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
#include <opm/common/utility/SaltElectrolytes.hpp>

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

BOOST_AUTO_TEST_CASE(DensityAndViscosityShareWaterDensity)
{
    // The fluid state path computes water density once for both properties; it must give
    // the same values as the separate calls, which compute it each
    using DM = Pvt::DensityModel;
    using VM = Pvt::ViscosityModel;
    const double c = 150.0;
    for (const auto density : {DM::MaoDuanYoung, DM::LaliberteCooper}) {
        for (const auto viscosity : {VM::MaoDuan, VM::Laliberte}) {
            auto pvt = directPvt(density, viscosity, Pvt::EnthalpyModel::Driesner);
            for (const bool clamp : {false, true}) {
                pvt.setClampToFitRange(clamp);
                for (const double T : {290.0, 340.0, 420.0}) {
                    for (const double p : {2e6, 25e6}) {
                        const CellFluidState<double> fs {T, p, c};
                        const auto [b, mu] = pvt.inverseFormationVolumeFactorAndViscosity(fs, 0);
                        BOOST_CHECK_EQUAL(b, pvt.inverseFormationVolumeFactor(0, T, p, 0.0, c));
                        BOOST_CHECK_EQUAL(mu, pvt.viscosity(0, T, p, 0.0, c));
                    }
                }
            }
        }
    }
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

BOOST_AUTO_TEST_CASE(IonicStrengthAndMeanMolarMass)
{
    using Electrolyte = Pvt::Electrolyte;
    const double mNaCl = Opm::SaltElectrolytes<double>::molarMass(Electrolyte::NaCl);
    const double mCaCl2 = Opm::SaltElectrolytes<double>::molarMass(Electrolyte::CaCl2);
    const double mLiCl = Opm::SaltElectrolytes<double>::molarMass(Electrolyte::LiCl);

    // 1 mol/kg NaCl: I = 1; 1 mol/kg CaCl2: I = 3
    Pvt::ElectrolyteList<double> nacl;
    nacl.push_back({Electrolyte::NaCl, mNaCl / (1.0 + mNaCl)});
    BOOST_CHECK_CLOSE(Brine::ionicStrength(nacl), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(Brine::meanSaltMolarMass(nacl), mNaCl, 1e-10);

    Pvt::ElectrolyteList<double> cacl2;
    cacl2.push_back({Electrolyte::CaCl2, mCaCl2 / (1.0 + mCaCl2)});
    BOOST_CHECK_CLOSE(Brine::ionicStrength(cacl2), 3.0, 1e-10);

    // Equal masses of NaCl and LiCl
    Pvt::ElectrolyteList<double> mixed;
    mixed.push_back({Electrolyte::NaCl, 0.05});
    mixed.push_back({Electrolyte::LiCl, 0.05});
    BOOST_CHECK_CLOSE(Brine::meanSaltMolarMass(mixed), 0.1 / (0.05 / mNaCl + 0.05 / mLiCl), 1e-10);

    BOOST_CHECK_EQUAL(Brine::meanSaltMolarMass(Pvt::ElectrolyteList<double> {}), 0.0);
    BOOST_CHECK_EQUAL(Brine::ionicStrength(Pvt::ElectrolyteList<double> {}), 0.0);
}

BOOST_AUTO_TEST_CASE(RegionIonicStrength)
{
    const Setup setup(deckWith("", salinitc));
    Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
    const auto& pvt = brinePvt(multiplexer, setup);

    // Region 1 at its SALINITC salt: I = 1/2 sum m z^2 from the SALINITC molalities
    //   NA 2.0, K 0.05, CA 0.4, MG 0.05, CL 2.97, SR 0.01
    const double expected = 0.5 * (2.0 + 0.05 + 4 * 0.4 + 4 * 0.05 + 2.97 + 4 * 0.01);
    const double c = pvt.salinity(0).sum() * pvt.waterReferenceDensity(0);
    BOOST_CHECK_CLOSE(pvt.ionicStrength(0, c), expected, 1e-6);

    // The same through the cell path
    const auto cell = pvt.cellFractions(pvt.salinity(0));
    BOOST_CHECK_CLOSE(pvt.ionicStrength(0, c, &cell), expected, 1e-4);

    // Mostly NaCl with some CaCl2, so between the two molar masses
    const double mw = pvt.meanSaltMolarMass(0, c);
    BOOST_CHECK(mw > Opm::SaltElectrolytes<double>::molarMass(Pvt::Electrolyte::NaCl));
    BOOST_CHECK(mw < Opm::SaltElectrolytes<double>::molarMass(Pvt::Electrolyte::CaCl2));

    // No salt
    BOOST_CHECK_EQUAL(pvt.ionicStrength(0, 0.0), 0.0);
}

BOOST_AUTO_TEST_CASE(CalibrationFactorGetters)
{
    auto pvt = directPvt(Pvt::DensityModel::MaoDuanYoung,
                         Pvt::ViscosityModel::MaoDuan,
                         Pvt::EnthalpyModel::Driesner);
    BOOST_CHECK_EQUAL(pvt.densityFactor(0, tRes), 1.0);
    BOOST_CHECK_EQUAL(pvt.viscosityFactor(0, tRes), 1.0);
    pvt.setCalibrationFactors(0, 1.03, 0.95);
    BOOST_CHECK_EQUAL(pvt.densityFactor(0, tRes), 1.03);
    BOOST_CHECK_EQUAL(pvt.viscosityFactor(0, tRes), 0.95);
}

BOOST_AUTO_TEST_CASE(LabCalibration)
{
    using Eval = Opm::DenseAd::Evaluation<double, 1>;
    auto pvt = directPvt(Pvt::DensityModel::MaoDuanYoung,
                         Pvt::ViscosityModel::MaoDuan,
                         Pvt::EnthalpyModel::Driesner);
    // The SALINITC brine, whatever the surface density in use
    const auto salt = [&] { return pvt.salinity(0).sum() * pvt.waterReferenceDensity(0); };
    const auto b = [&](double T, double p) {
        return pvt.inverseFormationVolumeFactor(0, T, p, 0.0, salt());
    };
    const auto rho = [&](double T, double p) { return b(T, p) * pvt.waterReferenceDensity(0); };
    const auto mu = [&](double T, double p) { return pvt.viscosity(0, T, p, 0.0, salt()); };

    // Two points at 330 K average to a density factor of 1.02; 380 K gives 0.98
    constexpr double t1 = 330.0;
    constexpr double t2 = 380.0;
    const double rhoRef = pvt.waterReferenceDensity(0);
    const double b2 = b(t2, pRes);
    const double mu2 = mu(t2, pRes);
    const double muMid = mu(355.0, pRes);
    pvt.setLabData(0,
                   {{t2, pRes, 0.98 * rho(t2, pRes), 0.9 * mu2},
                    {t1, 1e7, 1.01 * rho(t1, 1e7), 1.1 * mu(t1, 1e7)},
                    {t1, pRes, 1.03 * rho(t1, pRes), 1.1 * mu(t1, pRes)}});
    pvt.initEnd();

    const auto& c = pvt.calibration(0);
    BOOST_REQUIRE_EQUAL(c.temperature.size(), 2U);
    BOOST_CHECK_CLOSE(c.density[0], 1.02, 1e-10);
    BOOST_CHECK_CLOSE(c.density[1], 0.98, 1e-10);
    BOOST_CHECK_CLOSE(c.viscosity[0], 1.1, 1e-10);
    BOOST_CHECK_CLOSE(c.viscosity[1], 0.9, 1e-10);

    // Linear between the points, constant beyond them
    BOOST_CHECK_CLOSE(pvt.densityFactor(0, 355.0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(pvt.densityFactor(0, 300.0), 1.02, 1e-10);
    BOOST_CHECK_CLOSE(pvt.densityFactor(0, 400.0), 0.98, 1e-10);
    const Eval T = Eval::createVariable(355.0, 0);
    BOOST_CHECK_CLOSE(pvt.densityFactor(0, T).derivative(0), -0.04 / (t2 - t1), 1e-8);

    // The surface density takes the factor at standard temperature, below the lab range
    BOOST_CHECK_CLOSE(pvt.waterReferenceDensity(0), 1.02 * rhoRef, 1e-10);
    BOOST_CHECK_CLOSE(b(t2, pRes), 0.98 / 1.02 * b2, 1e-10);
    BOOST_CHECK_CLOSE(mu(t2, pRes), 0.9 * mu2, 1e-10);
    BOOST_CHECK_CLOSE(mu(355.0, pRes), muMid, 1e-10);

    // Running initEnd again gives the same factors
    const auto first = c;
    pvt.initEnd();
    BOOST_CHECK(pvt.calibration(0) == first);

    BOOST_CHECK_THROW(pvt.setLabData(0, {{t1, pRes, -1.0, 1e-3}}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(LabCalibrationFromDeck)
{
    const Setup plain(deckWith("", salinitc));
    Opm::WaterPvtMultiplexer<double, true, true> plainMultiplexer;
    const auto& model = brinePvt(plainMultiplexer, plain);

    const Setup lab(deckWith("", salinitc + R"(
BRINELAB
-- T (C)  p (bar)  density (kg/m3)  viscosity (cP)
   90     250      1172             0.69
   140    450      1150             0.48 /
   25     1        1100             1.10 /
)"));
    Opm::WaterPvtMultiplexer<double, true, true> labMultiplexer;
    const auto& pvt = brinePvt(labMultiplexer, lab);

    const auto check
        = [&](unsigned region, std::size_t i, double tC, double pBar, double rhoLab, double muLab) {
              const double T = tC + 273.15;
              const double p = pBar * 1e5;
              const double rhoRef = model.waterReferenceDensity(region);
              const double c = model.salinity(region).sum() * rhoRef;
              const double rho = model.inverseFormationVolumeFactor(region, T, p, 0.0, c) * rhoRef;
              const double mu = model.viscosity(region, T, p, 0.0, c);
              const auto& cal = pvt.calibration(region);
              BOOST_CHECK_CLOSE(cal.temperature[i], T, 1e-10);
              BOOST_CHECK_CLOSE(cal.density[i], rhoLab / rho, 1e-8);
              BOOST_CHECK_CLOSE(cal.viscosity[i], muLab * 1e-3 / mu, 1e-8);
          };
    check(0, 0, 90.0, 250.0, 1172.0, 0.69);
    check(0, 1, 140.0, 450.0, 1150.0, 0.48);
    check(1, 0, 25.0, 1.0, 1100.0, 1.10);
}

namespace
{

Pvt
pvtFor(const Opm::SaltArray<double, Opm::SaltMolality>& molality,
       Pvt::DensityModel density,
       Pvt::ViscosityModel viscosity)
{
    Pvt pvt;
    pvt.setNumRegions(1);
    pvt.setSalinity(0, molality.convert_to<Opm::SaltMassFraction>());
    pvt.setModels(density, viscosity, Pvt::EnthalpyModel::Driesner);
    pvt.setEnableSaltConcentration(true);
    pvt.initEnd();
    return pvt;
}

// Salt concentration giving the region's ion ratios at the ionic strength
double
concentrationAt(const Pvt& pvt, double ionicStrength)
{
    const double c = pvt.salinity(0).sum() * pvt.waterReferenceDensity(0);
    const double i0 = pvt.ionicStrength(0, c);
    // Ionic strength is not linear in mass fraction, so step to it
    double lo = 0.0;
    double hi = 4.0 * c * ionicStrength / i0;
    for (int k = 0; k < 100; ++k) {
        const double mid = 0.5 * (lo + hi);
        (pvt.ionicStrength(0, mid) < ionicStrength ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

constexpr double kelvin = 273.15;

} // namespace

BOOST_AUTO_TEST_CASE(BasinFitRange)
{
    // NaCl and CaCl2 both over 5% of the salt mass
    const auto pvt = directPvt(Pvt::DensityModel::MaoDuanYoung,
                               Pvt::ViscosityModel::MaoDuan,
                               Pvt::EnthalpyModel::Driesner);

    // Density: NaCl to 573 K, 100 MPa, 6 mol/kg; CaCl2 to 523 K, 60 MPa, 6 mol/kg (I = 18)
    const auto& density = pvt.fitRange(0, Pvt::FitProperty::Density);
    BOOST_CHECK_CLOSE(density.minTemperature, kelvin, 1e-10);
    BOOST_CHECK_CLOSE(density.maxTemperature, 523.0, 1e-10);
    BOOST_CHECK_CLOSE(density.maxPressure, 60e6, 1e-10);
    BOOST_CHECK_CLOSE(density.maxIonicStrength, 6.0, 1e-10);
    BOOST_CHECK_EQUAL(density.maxSaltMassFraction, Pvt::FitRange::open);

    // Viscosity: NaCl to 623 K, 100 MPa, 6 mol/kg; the CaCl2 fit to 575 K, 60 MPa
    const auto& viscosity = pvt.fitRange(0, Pvt::FitProperty::Viscosity);
    BOOST_CHECK_CLOSE(viscosity.maxTemperature, 575.0, 1e-10);
    BOOST_CHECK_CLOSE(viscosity.maxPressure, 60e6, 1e-10);
    BOOST_CHECK_CLOSE(viscosity.maxIonicStrength, 6.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(SalarFitRange)
{
    const auto pvt = directPvt(Pvt::DensityModel::LaliberteCooper,
                               Pvt::ViscosityModel::Laliberte,
                               Pvt::EnthalpyModel::Laliberte);

    // Laliberte NaCl and CaCl2: density 0 to 140 and 0 to 100 C, w 0.26 and 0.513;
    // viscosity 5 to 154 and 0 to 100 C, w 0.264 and 0.51. No pressure limit.
    const auto& density = pvt.fitRange(0, Pvt::FitProperty::Density);
    BOOST_CHECK_CLOSE(density.minTemperature, kelvin, 1e-10);
    BOOST_CHECK_CLOSE(density.maxTemperature, 100.0 + kelvin, 1e-10);
    BOOST_CHECK_CLOSE(density.maxSaltMassFraction, 0.260, 1e-10);
    BOOST_CHECK_EQUAL(density.maxPressure, Pvt::FitRange::open);
    BOOST_CHECK_EQUAL(density.maxIonicStrength, Pvt::FitRange::open);

    const auto& viscosity = pvt.fitRange(0, Pvt::FitProperty::Viscosity);
    BOOST_CHECK_CLOSE(viscosity.minTemperature, 5.0 + kelvin, 1e-10);
    BOOST_CHECK_CLOSE(viscosity.maxTemperature, 100.0 + kelvin, 1e-10);
    BOOST_CHECK_CLOSE(viscosity.maxSaltMassFraction, 0.264, 1e-10);
}

BOOST_AUTO_TEST_CASE(FitRangeSkipsTraceSalts)
{
    // K2SO4 is about 1.5% of the salt mass, so its narrow Laliberte range is not applied
    Opm::SaltArray<double, Opm::SaltMolality> molality;
    molality[SaltIndex::NA] = 2.0;
    molality[SaltIndex::CL] = 2.0;
    molality[SaltIndex::K] = 0.02;
    molality[SaltIndex::SO4] = 0.01;
    const auto pvt
        = pvtFor(molality, Pvt::DensityModel::LaliberteCooper, Pvt::ViscosityModel::Laliberte);

    const auto& density = pvt.fitRange(0, Pvt::FitProperty::Density);
    BOOST_CHECK_CLOSE(density.maxTemperature, 140.0 + kelvin, 1e-10);
    BOOST_CHECK_CLOSE(density.maxSaltMassFraction, 0.260, 1e-10);
}

BOOST_AUTO_TEST_CASE(BasinFitRangeLaliberteBinary)
{
    // MgCl2 has no Mao and Duan viscosity, so its Laliberte binary at the ionic strength
    // limits the range: from 15 C and to w 0.386, which is I = 3 w / ((1 - w) M). Its upper
    // temperature is the 200 C of the Magomedov and Alkhasov check, not Laliberte's 70 C.
    Opm::SaltArray<double, Opm::SaltMolality> molality;
    molality[SaltIndex::NA] = 1.0;
    molality[SaltIndex::MG] = 0.5;
    molality[SaltIndex::CL] = 2.0;
    const auto pvt
        = pvtFor(molality, Pvt::DensityModel::MaoDuanYoung, Pvt::ViscosityModel::MaoDuan);

    const double mass = Opm::SaltElectrolytes<double>::molarMass(Pvt::Electrolyte::MgCl2);
    const double maxIonic = 3.0 * 0.386 / ((1.0 - 0.386) * mass);
    const auto& viscosity = pvt.fitRange(0, Pvt::FitProperty::Viscosity);
    BOOST_CHECK_CLOSE(viscosity.minTemperature, 15.0 + kelvin, 1e-10);
    BOOST_CHECK_CLOSE(viscosity.maxTemperature, 200.0 + kelvin, 1e-10);
    BOOST_CHECK_CLOSE(viscosity.maxPressure, 100e6, 1e-10);
    BOOST_CHECK_CLOSE(viscosity.maxIonicStrength, std::min(6.0, maxIonic), 1e-10);

    // MgCl2 density is Mao and Duan: 543 K, 30 MPa, 3 mol/kg (I = 9)
    const auto& density = pvt.fitRange(0, Pvt::FitProperty::Density);
    BOOST_CHECK_CLOSE(density.maxTemperature, 543.0, 1e-10);
    BOOST_CHECK_CLOSE(density.maxPressure, 30e6, 1e-10);
    BOOST_CHECK_CLOSE(density.maxIonicStrength, 6.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(SalarMgCl2KeepsLaliberteRange)
{
    // The extended MgCl2 range belongs to the basin viscosity only
    Opm::SaltArray<double, Opm::SaltMolality> molality;
    molality[SaltIndex::NA] = 1.0;
    molality[SaltIndex::MG] = 0.5;
    molality[SaltIndex::CL] = 2.0;
    const auto pvt
        = pvtFor(molality, Pvt::DensityModel::LaliberteCooper, Pvt::ViscosityModel::Laliberte);
    BOOST_CHECK_CLOSE(
        pvt.fitRange(0, Pvt::FitProperty::Viscosity).maxTemperature, 70.0 + kelvin, 1e-10);
}

BOOST_AUTO_TEST_CASE(OutOfRange)
{
    const auto pvt = directPvt(Pvt::DensityModel::MaoDuanYoung,
                               Pvt::ViscosityModel::MaoDuan,
                               Pvt::EnthalpyModel::Driesner);
    const double c = concentrationAt(pvt, 3.0);
    constexpr unsigned both = Pvt::densityOutOfRange | Pvt::viscosityOutOfRange;

    BOOST_CHECK_EQUAL(pvt.outOfRange(0, tRes, pRes, c), 0U);
    // Above CaCl2 density's 523 K, below the viscosity limit of 575 K
    BOOST_CHECK_EQUAL(pvt.outOfRange(0, 540.0, pRes, c), Pvt::densityOutOfRange);
    BOOST_CHECK_EQUAL(pvt.outOfRange(0, 600.0, pRes, c), both);
    BOOST_CHECK_EQUAL(pvt.outOfRange(0, 270.0, pRes, c), both);
    BOOST_CHECK_EQUAL(pvt.outOfRange(0, tRes, 70e6, c), both);
    BOOST_CHECK_EQUAL(pvt.outOfRange(0, tRes, pRes, concentrationAt(pvt, 6.5)), both);
    BOOST_CHECK_EQUAL(pvt.outOfRange(0, tRes, pRes, concentrationAt(pvt, 5.9)), 0U);

    // The cell path takes the cell's ionic strength
    const auto cell = pvt.cellFractions(pvt.salinity(0));
    BOOST_CHECK_EQUAL(pvt.outOfRange(0, 540.0, pRes, c, &cell), Pvt::densityOutOfRange);
}

BOOST_AUTO_TEST_CASE(ClampToFitRange)
{
    using Eval = Opm::DenseAd::Evaluation<double, 1>;
    auto pvt = directPvt(Pvt::DensityModel::MaoDuanYoung,
                         Pvt::ViscosityModel::MaoDuan,
                         Pvt::EnthalpyModel::Driesner);
    const double c = concentrationAt(pvt, 3.0);
    const double bIn = pvt.inverseFormationVolumeFactor(0, tRes, pRes, 0.0, c);
    const double b523 = pvt.inverseFormationVolumeFactor(0, 523.0, pRes, 0.0, c);
    const double b540 = pvt.inverseFormationVolumeFactor(0, 540.0, pRes, 0.0, c);
    const double mu540 = pvt.viscosity(0, 540.0, pRes, 0.0, c);
    const double b70 = pvt.inverseFormationVolumeFactor(0, tRes, 70e6, 0.0, c);
    const double b60 = pvt.inverseFormationVolumeFactor(0, tRes, 60e6, 0.0, c);
    BOOST_CHECK(std::abs(b540 - b523) > 1e-4 * b523);

    pvt.setClampToFitRange(true);
    BOOST_CHECK_EQUAL(pvt.inverseFormationVolumeFactor(0, tRes, pRes, 0.0, c), bIn);
    BOOST_CHECK_CLOSE(pvt.inverseFormationVolumeFactor(0, 540.0, pRes, 0.0, c), b523, 1e-12);
    BOOST_CHECK_CLOSE(pvt.inverseFormationVolumeFactor(0, tRes, 70e6, 0.0, c), b60, 1e-12);
    BOOST_CHECK(std::abs(b70 - b60) > 1e-6 * b60);

    // Viscosity has its own range, so 540 K is not clamped there
    BOOST_CHECK_EQUAL(pvt.viscosity(0, 540.0, pRes, 0.0, c), mu540);

    // Beyond the range the clamped variable has no derivative
    const Eval T = Eval::createVariable(540.0, 0);
    const Eval b = pvt.inverseFormationVolumeFactor(0, T, Eval(pRes), Eval(0.0), Eval(c));
    BOOST_CHECK_EQUAL(b.derivative(0), 0.0);

    // Salinity is not clamped
    const double cHigh = concentrationAt(pvt, 6.5);
    pvt.setClampToFitRange(false);
    const double bHigh = pvt.inverseFormationVolumeFactor(0, tRes, pRes, 0.0, cHigh);
    pvt.setClampToFitRange(true);
    BOOST_CHECK_EQUAL(pvt.inverseFormationVolumeFactor(0, tRes, pRes, 0.0, cHigh), bHigh);
}

BOOST_AUTO_TEST_CASE(RangePolicyFromDeck)
{
    {
        const Setup setup(deckWith("1* 1* 1* 1* 1* 1* CLAMP", salinitc));
        Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
        BOOST_CHECK(brinePvt(multiplexer, setup).clampToFitRange());
    }
    {
        const Setup setup(deckWith("", salinitc));
        Opm::WaterPvtMultiplexer<double, true, true> multiplexer;
        BOOST_CHECK(!brinePvt(multiplexer, setup).clampToFitRange());
    }
}
