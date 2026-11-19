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

#include <config.h>
#include <opm/material/fluidsystems/blackoilpvt/BrineMulticompPvt.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <string>

#include <stdexcept>

namespace Opm
{

template <class Scalar>
void
BrineMulticompPvt<Scalar>::initFromState(const EclipseState& eclState, const Schedule&)
{
    const auto& saltmc = eclState.runspec().multiCompSalt();
    const auto density = saltmc.densityModel() == Saltmc::DensityModel::MaoDuanYoung
        ? DensityModel::MaoDuanYoung
        : DensityModel::LaliberteCooper;
    const auto viscosity = saltmc.viscosityModel() == Saltmc::ViscosityModel::MaoDuan
        ? ViscosityModel::MaoDuan
        : ViscosityModel::Laliberte;
    EnthalpyModel enthalpy = EnthalpyModel::Driesner;
    switch (saltmc.enthalpyModel()) {
    case Saltmc::EnthalpyModel::Driesner:
        enthalpy = EnthalpyModel::Driesner;
        break;
    case Saltmc::EnthalpyModel::Laliberte:
        enthalpy = EnthalpyModel::Laliberte;
        break;
    case Saltmc::EnthalpyModel::Table:
        enthalpy = EnthalpyModel::Table;
        break;
    case Saltmc::EnthalpyModel::Thco2mix:
        throw std::logic_error("The THCO2MIX enthalpy belongs to the CO2STORE water model");
    }
    setModels(density, viscosity, enthalpy, PairingMethod::Proportional);
    setLaliberteCooperWater(saltmc.LaliberteCooperDensity());
    setClampToFitRange(saltmc.rangePolicy() == Saltmc::RangePolicy::Clamp);
    setEnableSaltConcentration(eclState.runspec().phases().active(Phase::BRINE));

    const auto& tables = eclState.getTableManager();
    setReferenceConditions(tables.stCond().temperature, tables.stCond().pressure);
    if (!tables.getPvtwTable().empty() || !tables.getPvtwSaltTables().empty()) {
        OpmLog::warning("SALTMC: water properties come from the brine models; "
                        "PVTW and PVTWSALT are ignored.");
    }

    const std::size_t regions = eclState.runspec().tabdims().getNumPVTTables();
    setNumRegions(regions);
    const auto& config = eclState.getCo2StoreConfig();
    for (std::size_t regionIdx = 0; regionIdx < regions; ++regionIdx) {
        const auto& salinity = config.numSaltRegions() > regionIdx
            ? config.saltComponents(regionIdx)
            : config.saltComponents();
        setSalinity(regionIdx, SaltArray<Scalar, SaltMassFraction>(salinity));
    }

    const auto& simulation = eclState.getSimulationConfig();
    if (enthalpy == EnthalpyModel::Table && (simulation.isThermal() || simulation.isTemp())) {
        if (!tables.hasTables("SPECHEAT")) {
            throw std::invalid_argument("SALTMC ENTHALPY_MODEL TABLE needs SPECHEAT");
        }
        // Integrate the tabulated heat capacity, as for thermal water
        for (std::size_t regionIdx = 0; regionIdx < regions; ++regionIdx) {
            const auto& table = tables.getSpecheatTables()[regionIdx];
            const auto& temperature = table.getColumn("TEMPERATURE");
            const auto& cv = table.getColumn("CV_WATER");
            std::vector<Scalar> u(temperature.size());
            u[0] = temperature[0] * cv[0];
            for (std::size_t i = 1; i < temperature.size(); ++i) {
                u[i] = u[i - 1] + 0.5 * (cv[i - 1] + cv[i]) * (temperature[i] - temperature[i - 1]);
            }
            TabulatedOneDFunction curve;
            curve.setXYContainers(temperature.vectorCopy(), u);
            setInternalEnergyCurve(regionIdx, curve);
        }
    }

    initEnd();

    const auto describe = [](const FitRange& r) {
        std::string text = r.maxTemperature < FitRange::open
            ? fmt::format(
                  "{:.0f} to {:.0f} C", r.minTemperature - 273.15, r.maxTemperature - 273.15)
            : fmt::format("above {:.0f} C", r.minTemperature - 273.15);
        if (r.maxPressure < FitRange::open) {
            text += fmt::format(", {:g} MPa", r.maxPressure * 1e-6);
        }
        if (r.maxIonicStrength < FitRange::open) {
            text += fmt::format(", ionic strength {:.3g} mol/kg", r.maxIonicStrength);
        }
        if (r.maxSaltMassFraction < FitRange::open) {
            text += fmt::format(", salt mass fraction {:.3g}", r.maxSaltMassFraction);
        }
        return text;
    };
    for (std::size_t regionIdx = 0; regionIdx < regions; ++regionIdx) {
        OpmLog::info(fmt::format("SALTMC region {}: brine surface density {:.2f} kg/m3, "
                                 "SALINITC salt concentration {:.3f} kg/sm3",
                                 regionIdx + 1,
                                 waterReferenceDensity(regionIdx),
                                 saltMassFraction_[regionIdx] * waterReferenceDensity(regionIdx)));
        OpmLog::info(fmt::format("SALTMC region {}: density fitted for {}; viscosity for {}",
                                 regionIdx + 1,
                                 describe(fitRange(regionIdx, FitProperty::Density)),
                                 describe(fitRange(regionIdx, FitProperty::Viscosity))));
    }
}

template <class Scalar>
void
BrineMulticompPvt<Scalar>::setNumRegions(std::size_t numRegions)
{
    unitElectrolytes_.resize(numRegions);
    salinity_.resize(numRegions);
    saltMassFraction_.resize(numRegions, 0.0);
    referenceDensity_.resize(numRegions, 0.0);
    densityFactor_.resize(numRegions, 1.0);
    viscosityFactor_.resize(numRegions, 1.0);
    fitRanges_.resize(numRegions);
}

template <class Scalar>
void
BrineMulticompPvt<Scalar>::setSalinity(unsigned regionIdx,
                                       const SaltArray<Scalar, SaltMassFraction>& salinity)
{
    salinity_[regionIdx] = salinity;
}

template <class Scalar>
void
BrineMulticompPvt<Scalar>::setModels(DensityModel density,
                                     ViscosityModel viscosity,
                                     EnthalpyModel enthalpy,
                                     PairingMethod pairing)
{
    densityModel_ = density;
    viscosityModel_ = viscosity;
    enthalpyModel_ = enthalpy;
    pairing_ = pairing;
}

template <class Scalar>
void
BrineMulticompPvt<Scalar>::setReferenceConditions(Scalar temperature, Scalar pressure)
{
    referenceTemperature_ = temperature;
    referencePressure_ = pressure;
}

template <class Scalar>
void
BrineMulticompPvt<Scalar>::setCalibrationFactors(unsigned regionIdx,
                                                 Scalar density,
                                                 Scalar viscosity)
{
    if (!(density > 0.0) || !(viscosity > 0.0)) {
        throw std::invalid_argument("Brine calibration factors must be positive");
    }
    densityFactor_[regionIdx] = density;
    viscosityFactor_[regionIdx] = viscosity;
}

template <class Scalar>
void
BrineMulticompPvt<Scalar>::setInternalEnergyCurve(unsigned regionIdx,
                                                  const TabulatedOneDFunction& curve)
{
    internalEnergyCurves_.resize(salinity_.size());
    internalEnergyCurves_[regionIdx] = curve;
}

template <class Scalar>
void
BrineMulticompPvt<Scalar>::initEnd()
{
    using Electrolytes = SaltElectrolytes<Scalar>;
    for (std::size_t regionIdx = 0; regionIdx < salinity_.size(); ++regionIdx) {
        const Scalar total = salinity_[regionIdx].sum();
        saltMassFraction_[regionIdx] = total;

        ElectrolyteList<Scalar> unit;
        if (total > 0.0) {
            for (const auto& [electrolyte, w] :
                 Electrolytes::fromSaltComponents(salinity_[regionIdx], pairing_)) {
                unit.push_back({electrolyte, w / total});
            }
        } else {
            // No composition given: salt from BRINE counts as NaCl
            unit.push_back({Electrolyte::NaCl, 1.0});
        }
        unitElectrolytes_[regionIdx] = unit;
        computeFitRanges_(regionIdx);

        ElectrolyteList<Scalar> salts;
        for (const auto& [electrolyte, fraction] : unit) {
            salts.push_back({electrolyte, fraction * total});
        }
        referenceDensity_[regionIdx]
            = density_(regionIdx, referenceTemperature_, referencePressure_, salts);
    }
}

// Mao and Duan evaluate each binary at the brine's ionic strength, so their molality limit
// is an ionic strength limit. Laliberte density and viscosity take the total salt mass
// fraction, except inside the Mao and Duan viscosity, where a binary is taken at the ionic
// strength.
template <class Scalar>
void
BrineMulticompPvt<Scalar>::computeFitRanges_(unsigned regionIdx)
{
    using MaoDuan = MaoDuanBrine<Scalar>;
    auto& density = fitRanges_[regionIdx][static_cast<std::size_t>(FitProperty::Density)];
    auto& viscosity = fitRanges_[regionIdx][static_cast<std::size_t>(FitProperty::Viscosity)];
    density = FitRange {};
    viscosity = FitRange {};
    const auto narrowTemperature = [](FitRange& range, Scalar minT, Scalar maxT) {
        range.minTemperature = std::max(range.minTemperature, minT);
        range.maxTemperature = std::min(range.maxTemperature, maxT);
    };

    for (const auto& [electrolyte, share] : unitElectrolytes_[regionIdx]) {
        if (share < minRangeShare) {
            continue;
        }
        const Scalar ionic = Brine::ionicStrengthPerMole(electrolyte);

        if (densityModel_ == DensityModel::MaoDuanYoung && MaoDuan::covers(electrolyte)) {
            const auto r = MaoDuan::range(electrolyte);
            narrowTemperature(density, r.minT, r.maxT);
            density.maxPressure = std::min(density.maxPressure, Scalar(r.maxP * 1e6));
            density.maxIonicStrength = std::min(density.maxIonicStrength, r.maxMolality * ionic);
        } else {
            const auto r = Brine::laliberteDensityRange(electrolyte);
            narrowTemperature(density, r.minT, r.maxT);
            density.maxSaltMassFraction = std::min(density.maxSaltMassFraction, r.maxMassFraction);
        }

        if (viscosityModel_ == ViscosityModel::Laliberte) {
            const auto r = Brine::laliberteViscosityRange(electrolyte);
            narrowTemperature(viscosity, r.minT, r.maxT);
            viscosity.maxSaltMassFraction
                = std::min(viscosity.maxSaltMassFraction, r.maxMassFraction);
        } else if (MaoDuan::coversViscosity(electrolyte)) {
            const auto r = MaoDuan::viscosityRange(electrolyte);
            narrowTemperature(viscosity, r.minT, r.maxT);
            viscosity.maxPressure = std::min(viscosity.maxPressure, Scalar(r.maxP * 1e6));
            viscosity.maxIonicStrength
                = std::min(viscosity.maxIonicStrength, r.maxMolality * ionic);
        } else if (ionic > 0.0) {
            auto r = Brine::laliberteViscosityRange(electrolyte);
            if (electrolyte == Electrolyte::MgCl2) {
                // Checked against Magomedov and Alkhasov (2007), 75 to 200 C and to 100 MPa,
                // with the Mao and Duan water viscosity carrying pressure
                r.maxT = 200.0 + 273.15;
                viscosity.maxPressure = std::min(viscosity.maxPressure, Scalar {100e6});
            }
            const Scalar w = r.maxMassFraction;
            const Scalar maxMolality
                = w / ((1.0 - w) * SaltElectrolytes<Scalar>::molarMass(electrolyte));
            narrowTemperature(viscosity, r.minT, r.maxT);
            viscosity.maxIonicStrength = std::min(viscosity.maxIonicStrength, maxMolality * ionic);
        }
    }
}

template <class Scalar>
ElectrolyteFractions
BrineMulticompPvt<Scalar>::cellFractions(const SaltArray<Scalar, SaltMassFraction>& ions) const
{
    ElectrolyteFractions result {};
    const Scalar total = ions.sum();
    if (!(total > 0.0)) {
        return result;
    }
    for (const auto& [electrolyte, w] :
         SaltElectrolytes<Scalar>::fromSaltComponents(ions, pairing_)) {
        result[static_cast<std::size_t>(electrolyte)] = static_cast<float>(w / total);
    }
    return result;
}

template <class Scalar>
bool
BrineMulticompPvt<Scalar>::operator==(const BrineMulticompPvt& data) const
{
    return unitElectrolytes_ == data.unitElectrolytes_ && salinity_ == data.salinity_
        && saltMassFraction_ == data.saltMassFraction_
        && referenceDensity_ == data.referenceDensity_ && densityFactor_ == data.densityFactor_
        && viscosityFactor_ == data.viscosityFactor_
        && internalEnergyCurves_ == data.internalEnergyCurves_
        && densityModel_ == data.densityModel_ && viscosityModel_ == data.viscosityModel_
        && enthalpyModel_ == data.enthalpyModel_ && pairing_ == data.pairing_
        && laliberteCooperWater_ == data.laliberteCooperWater_
        && enableSaltConcentration_ == data.enableSaltConcentration_
        && referenceTemperature_ == data.referenceTemperature_
        && referencePressure_ == data.referencePressure_ && fitRanges_ == data.fitRanges_
        && clampToFitRange_ == data.clampToFitRange_;
}

template class BrineMulticompPvt<double>;
template class BrineMulticompPvt<float>;

} // namespace Opm
