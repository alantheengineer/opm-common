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

    for (std::size_t regionIdx = 0; regionIdx < regions; ++regionIdx) {
        OpmLog::info(fmt::format("SALTMC region {}: brine surface density {:.2f} kg/m3, "
                                 "SALINITC salt concentration {:.3f} kg/sm3",
                                 regionIdx + 1,
                                 waterReferenceDensity(regionIdx),
                                 saltMassFraction_[regionIdx] * waterReferenceDensity(regionIdx)));
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

        ElectrolyteList<Scalar> salts;
        for (const auto& [electrolyte, fraction] : unit) {
            salts.push_back({electrolyte, fraction * total});
        }
        referenceDensity_[regionIdx] = density_(referenceTemperature_, referencePressure_, salts);
    }
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
        && referencePressure_ == data.referencePressure_;
}

template class BrineMulticompPvt<double>;
template class BrineMulticompPvt<float>;

} // namespace Opm
