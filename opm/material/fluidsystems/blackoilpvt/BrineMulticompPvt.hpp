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
/*!
 * \file
 * \copydoc Opm::BrineMulticompPvt
 */
#ifndef OPM_BRINE_MULTICOMP_PVT_HPP
#define OPM_BRINE_MULTICOMP_PVT_HPP

#include <opm/common/utility/SaltArray.hpp>
#include <opm/common/utility/SaltElectrolytes.hpp>

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/components/BrineDynamic.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/fluidsystems/BlackOilFunctions.hpp>

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Opm
{

class EclipseState;
class Schedule;

/*!
 * \brief Water phase PVT for a brine of several salts, set by SALTMC and SALINITC.
 *
 * Each PVT region has fixed ion ratios. The total salt comes from the BRINE salt concentration
 * when that is active, otherwise from SALINITC. The concentration is salt mass per surface
 * volume of brine, so the salt mass fraction is the concentration over the region's reference
 * density, the density of the SALINITC brine at standard conditions.
 *
 * Ions are paired into electrolytes once per region. Pairing scales with the total salt, so
 * each call only scales the stored electrolyte fractions, with no heap allocation.
 */
template <class Scalar>
class BrineMulticompPvt
{
public:
    using H2O = ::Opm::H2O<Scalar>;
    using Brine = BrineDynamic<Scalar, H2O>;
    using PairingMethod = typename Brine::PairingMethod;
    using Electrolyte = typename Brine::Electrolyte;
    template <class Evaluation>
    using ElectrolyteList = typename Brine::template ElectrolyteList<Evaluation>;
    using TabulatedOneDFunction = Tabulated1DFunction<Scalar>;

    enum class DensityModel { MaoDuanYoung, LaliberteCooper };
    enum class ViscosityModel { MaoDuan, Laliberte };
    enum class EnthalpyModel { Driesner, Laliberte, Table };

    BrineMulticompPvt() = default;

    void initFromState(const EclipseState& eclState, const Schedule& schedule);

    void setNumRegions(std::size_t numRegions);

    /*!
     * \brief Set the ion composition of a region.
     *
     * \param salinity Ion mass fractions in the brine [kg ion/kg brine]
     */
    void setSalinity(unsigned regionIdx, const SaltArray<Scalar, SaltMassFraction>& salinity);

    void setModels(DensityModel density,
                   ViscosityModel viscosity,
                   EnthalpyModel enthalpy,
                   PairingMethod pairing = PairingMethod::Proportional);

    /// Use the Laliberte and Cooper pure water density in that density model
    void setLaliberteCooperWater(bool yesno)
    {
        laliberteCooperWater_ = yesno;
    }

    /// Take the total salt from the BRINE salt concentration
    void setEnableSaltConcentration(bool yesno)
    {
        enableSaltConcentration_ = yesno;
    }

    void setReferenceConditions(Scalar temperature, Scalar pressure);

    /// Factors on the model density and viscosity, from lab calibration
    void setCalibrationFactors(unsigned regionIdx, Scalar density, Scalar viscosity);

    /// Internal energy against temperature, for the TABLE enthalpy model
    void setInternalEnergyCurve(unsigned regionIdx, const TabulatedOneDFunction& curve);

    /// Finish setup: pair the ions and compute the reference densities
    void initEnd();

    unsigned numRegions() const
    {
        return static_cast<unsigned>(referenceDensity_.size());
    }

    void setVapPars(const Scalar, const Scalar)
    {
    }

    Scalar hVap(unsigned) const
    {
        return 0.0;
    }

    /// The internal energy is evaluated directly, not by a mixing model
    static constexpr bool mixingEnergy()
    {
        return false;
    }

    Scalar waterReferenceDensity(unsigned regionIdx) const
    {
        return densityFactor_[regionIdx] * referenceDensity_[regionIdx];
    }

    /// Salt mass fraction [kg salt/kg brine] from the BRINE salt concentration [kg/sm3]
    template <class Evaluation>
    Evaluation saltMassFraction(unsigned regionIdx, const Evaluation& saltConcentration) const
    {
        if (!enableSaltConcentration_) {
            return Evaluation(saltMassFraction_[regionIdx]);
        }
        return saltConcentration / waterReferenceDensity(regionIdx);
    }

    /// Electrolytes and their mass fractions at the given salt concentration
    template <class Evaluation>
    ElectrolyteList<Evaluation> electrolytes(unsigned regionIdx,
                                             const Evaluation& saltConcentration) const
    {
        const Evaluation w = saltMassFraction(regionIdx, saltConcentration);
        ElectrolyteList<Evaluation> result;
        for (const auto& [electrolyte, fraction] : unitElectrolytes_[regionIdx]) {
            result.push_back({electrolyte, fraction * w});
        }
        return result;
    }

    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure,
                              const Evaluation& /*Rsw*/,
                              const Evaluation& saltConcentration) const
    {
        if (enthalpyModel_ == EnthalpyModel::Table) {
            if (internalEnergyCurves_.empty()) {
                throw std::runtime_error("SALTMC ENTHALPY_MODEL TABLE needs SPECHEAT");
            }
            return internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true);
        }
        const auto salts = electrolytes(regionIdx, saltConcentration);
        const Evaluation h = enthalpyModel_ == EnthalpyModel::Driesner
            ? Brine::liquidEnthalpyLaliberteDriesner(temperature, pressure, salts)
            : Brine::liquidEnthalpyLaliberte(temperature, pressure, salts);
        return h - pressure / density_(temperature, pressure, salts);
    }

    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& /*Rsw*/,
                         const Evaluation& saltConcentration) const
    {
        return saturatedViscosity(regionIdx, temperature, pressure, saltConcentration);
    }

    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure,
                                  const Evaluation& saltConcentration) const
    {
        const auto salts = electrolytes(regionIdx, saltConcentration);
        return viscosityFactor_[regionIdx] * viscosity_(temperature, pressure, salts);
    }

    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& /*Rsw*/,
                                            const Evaluation& saltConcentration) const
    {
        return saturatedInverseFormationVolumeFactor(
            regionIdx, temperature, pressure, saltConcentration);
    }

    // The density factor scales both densities, so it cancels here
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure,
                                                     const Evaluation& saltConcentration) const
    {
        const auto salts = electrolytes(regionIdx, saltConcentration);
        return density_(temperature, pressure, salts) / referenceDensity_[regionIdx];
    }

    template <class FluidState, class LhsEval = typename FluidState::ValueType>
    std::pair<LhsEval, LhsEval>
    inverseFormationVolumeFactorAndViscosity(const FluidState& fluidState, unsigned regionIdx) const
    {
        const LhsEval& T = decay<LhsEval>(fluidState.temperature(FluidState::waterPhaseIdx));
        const LhsEval& p = decay<LhsEval>(fluidState.pressure(FluidState::waterPhaseIdx));
        const LhsEval& saltConcentration
            = BlackOil::template getSaltConcentration_<FluidState, LhsEval>(fluidState, regionIdx);
        const auto salts = electrolytes(regionIdx, saltConcentration);
        return {density_(T, p, salts) / referenceDensity_[regionIdx],
                viscosityFactor_[regionIdx] * viscosity_(T, p, salts)};
    }

    template <class Evaluation>
    Evaluation saturationPressure(unsigned /*regionIdx*/,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& /*Rs*/,
                                  const Evaluation& /*saltConcentration*/) const
    {
        return 0.0; /* dead water */
    }

    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned /*regionIdx*/,
                                             const Evaluation& /*temperature*/,
                                             const Evaluation& /*pressure*/,
                                             const Evaluation& /*saltConcentration*/) const
    {
        return 0.0; /* dead water */
    }

    template <class Evaluation>
    Evaluation diffusionCoefficient(const Evaluation& /*temperature*/,
                                    const Evaluation& /*pressure*/,
                                    unsigned /*compIdx*/,
                                    unsigned /*regionIdx*/ = 0) const
    {
        throw std::runtime_error("Not implemented: The PVT model does not provide "
                                 "a diffusionCoefficient()");
    }

    DensityModel densityModel() const
    {
        return densityModel_;
    }

    ViscosityModel viscosityModel() const
    {
        return viscosityModel_;
    }

    EnthalpyModel enthalpyModel() const
    {
        return enthalpyModel_;
    }

    PairingMethod pairing() const
    {
        return pairing_;
    }

    bool operator==(const BrineMulticompPvt& data) const;

private:
    template <class Evaluation>
    Evaluation density_(const Evaluation& temperature,
                        const Evaluation& pressure,
                        const ElectrolyteList<Evaluation>& salts) const
    {
        if (densityModel_ == DensityModel::MaoDuanYoung) {
            return Brine::liquidDensityMaoDuan(temperature, pressure, salts);
        }
        const Evaluation rhoW = laliberteCooperWater_
            ? Brine::liquidDensityPureWaterLaliberteCooper(temperature)
            : H2O::liquidDensity(temperature, pressure, /*extrapolate=*/true);
        if (salts.empty()) {
            return rhoW;
        }
        return Brine::liquidDensityLaliberteCooper(temperature, salts, rhoW);
    }

    template <class Evaluation>
    Evaluation viscosity_(const Evaluation& temperature,
                          const Evaluation& pressure,
                          const ElectrolyteList<Evaluation>& salts) const
    {
        if (viscosityModel_ == ViscosityModel::MaoDuan) {
            return Brine::liquidViscosityMaoDuan(temperature, pressure, salts);
        }
        return Brine::liquidViscosityLaliberte(temperature, pressure, salts);
    }

    // Per region: electrolyte mass fractions per unit salt mass fraction
    std::vector<ElectrolyteList<Scalar>> unitElectrolytes_ {};
    std::vector<SaltArray<Scalar, SaltMassFraction>> salinity_ {};
    std::vector<Scalar> saltMassFraction_ {};
    std::vector<Scalar> referenceDensity_ {};
    std::vector<Scalar> densityFactor_ {};
    std::vector<Scalar> viscosityFactor_ {};
    std::vector<TabulatedOneDFunction> internalEnergyCurves_ {};

    DensityModel densityModel_ {DensityModel::MaoDuanYoung};
    ViscosityModel viscosityModel_ {ViscosityModel::MaoDuan};
    EnthalpyModel enthalpyModel_ {EnthalpyModel::Driesner};
    PairingMethod pairing_ {PairingMethod::Proportional};
    bool laliberteCooperWater_ {false};
    bool enableSaltConcentration_ {false};
    Scalar referenceTemperature_ {288.71};
    Scalar referencePressure_ {1.01325e5};
};

} // namespace Opm

#endif
