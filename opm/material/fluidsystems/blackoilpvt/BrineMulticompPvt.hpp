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

#include <opm/common/utility/ElectrolyteFractions.hpp>
#include <opm/common/utility/SaltArray.hpp>
#include <opm/common/utility/SaltElectrolytes.hpp>

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/components/BrineDynamic.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/fluidsystems/BlackOilFunctions.hpp>

#include <array>
#include <cstddef>
#include <limits>
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
 *
 * In TRACER mode a cell can carry its own electrolyte fractions, built by cellFractions()
 * from the ion tracers. Calls given those fractions use them; the rest use the region's.
 *
 * Each region has a fitted range per property: the overlap of the ranges of the models that
 * evaluate its main electrolytes. outOfRange() reports the properties a state falls outside
 * of. With clamping on, temperature and pressure are held inside the range before evaluation;
 * salinity is not clamped.
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

    enum class FitProperty { Density, Viscosity };

    /// Bits returned by outOfRange()
    static constexpr unsigned densityOutOfRange = 1;
    static constexpr unsigned viscosityOutOfRange = 2;

    /// Electrolytes with a smaller share of the salt mass are left out of the fitted range
    static constexpr Scalar minRangeShare = 0.05;

    /// Fitted range of one property in one region
    struct FitRange
    {
        static constexpr Scalar open = std::numeric_limits<Scalar>::max();

        Scalar minTemperature {0.0}; // K
        Scalar maxTemperature {open}; // K
        Scalar maxPressure {open}; // Pa
        Scalar maxIonicStrength {open}; // mol/kg water
        Scalar maxSaltMassFraction {open};

        bool operator==(const FitRange&) const = default;
    };

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

    /// Hold temperature and pressure inside the fitted range when evaluating
    void setClampToFitRange(bool yesno)
    {
        clampToFitRange_ = yesno;
    }

    bool clampToFitRange() const
    {
        return clampToFitRange_;
    }

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

    /// Calibration factor on density in use for a region
    Scalar densityFactor(unsigned regionIdx) const
    {
        return densityFactor_[regionIdx];
    }

    /// Calibration factor on viscosity in use for a region
    Scalar viscosityFactor(unsigned regionIdx) const
    {
        return viscosityFactor_[regionIdx];
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

    /// Ion mass fractions from SALINITC for a region [kg ion/kg brine]
    const SaltArray<Scalar, SaltMassFraction>& salinity(unsigned regionIdx) const
    {
        return salinity_[regionIdx];
    }

    /*!
     * \brief Electrolyte fractions per unit salt for a cell's ion composition.
     *
     * \param ions Ion mass fractions of the cell's brine [kg ion/kg brine]
     */
    ElectrolyteFractions cellFractions(const SaltArray<Scalar, SaltMassFraction>& ions) const;

    /*!
     * \brief Electrolytes and their mass fractions at the given salt concentration.
     *
     * \param cell The cell's electrolyte fractions, or nullptr for the region's
     */
    template <class Evaluation>
    ElectrolyteList<Evaluation> electrolytes(unsigned regionIdx,
                                             const Evaluation& saltConcentration,
                                             const ElectrolyteFractions* cell = nullptr) const
    {
        const Evaluation w = saltMassFraction(regionIdx, saltConcentration);
        ElectrolyteList<Evaluation> result;
        if (cell) {
            for (std::size_t i = 0; i < cell->size(); ++i) {
                if ((*cell)[i] > 0.0f) {
                    result.push_back({static_cast<Electrolyte>(i), Scalar((*cell)[i]) * w});
                }
            }
            return result;
        }
        for (const auto& [electrolyte, fraction] : unitElectrolytes_[regionIdx]) {
            result.push_back({electrolyte, fraction * w});
        }
        return result;
    }

    /// Ionic strength [mol/kg water]; cell fractions as for electrolytes()
    template <class Evaluation>
    Evaluation ionicStrength(unsigned regionIdx,
                             const Evaluation& saltConcentration,
                             const ElectrolyteFractions* cell = nullptr) const
    {
        return Brine::ionicStrength(electrolytes(regionIdx, saltConcentration, cell));
    }

    /// Mean molar mass of the dissolved salt [kg/mol]; cell fractions as for electrolytes()
    template <class Evaluation>
    Evaluation meanSaltMolarMass(unsigned regionIdx,
                                 const Evaluation& saltConcentration,
                                 const ElectrolyteFractions* cell = nullptr) const
    {
        return Brine::meanSaltMolarMass(electrolytes(regionIdx, saltConcentration, cell));
    }

    const FitRange& fitRange(unsigned regionIdx, FitProperty property) const
    {
        return fitRanges_[regionIdx][static_cast<std::size_t>(property)];
    }

    /// The properties a state is outside the fitted range of, as densityOutOfRange and
    /// viscosityOutOfRange bits; cell fractions as for electrolytes()
    unsigned outOfRange(unsigned regionIdx,
                        Scalar temperature,
                        Scalar pressure,
                        Scalar saltConcentration,
                        const ElectrolyteFractions* cell = nullptr) const
    {
        const Scalar w = saltMassFraction(regionIdx, saltConcentration);
        const Scalar ionic = ionicStrength(regionIdx, saltConcentration, cell);
        const auto outside = [&](const FitRange& r) {
            return temperature < r.minTemperature || temperature > r.maxTemperature
                || pressure > r.maxPressure || ionic > r.maxIonicStrength
                || w > r.maxSaltMassFraction;
        };
        unsigned flags = 0;
        if (outside(fitRange(regionIdx, FitProperty::Density))) {
            flags |= densityOutOfRange;
        }
        if (outside(fitRange(regionIdx, FitProperty::Viscosity))) {
            flags |= viscosityOutOfRange;
        }
        return flags;
    }

    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure,
                              const Evaluation& Rsw,
                              const Evaluation& saltConcentration) const
    {
        return internalEnergy(regionIdx, temperature, pressure, Rsw, saltConcentration, nullptr);
    }

    /// Internal energy with the cell's electrolyte fractions, or the region's for nullptr
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure,
                              const Evaluation& /*Rsw*/,
                              const Evaluation& saltConcentration,
                              const ElectrolyteFractions* cell) const
    {
        if (enthalpyModel_ == EnthalpyModel::Table) {
            if (internalEnergyCurves_.empty()) {
                throw std::runtime_error("SALTMC ENTHALPY_MODEL TABLE needs SPECHEAT");
            }
            return internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true);
        }
        const auto salts = electrolytes(regionIdx, saltConcentration, cell);
        const Evaluation h = enthalpyModel_ == EnthalpyModel::Driesner
            ? Brine::liquidEnthalpyLaliberteDriesner(temperature, pressure, salts)
            : Brine::liquidEnthalpyLaliberte(temperature, pressure, salts);
        return h - pressure / density_(regionIdx, temperature, pressure, salts);
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
        return viscosityFactor_[regionIdx] * viscosity_(regionIdx, temperature, pressure, salts);
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
        return density_(regionIdx, temperature, pressure, salts) / referenceDensity_[regionIdx];
    }

    template <class FluidState, class LhsEval = typename FluidState::ValueType>
    std::pair<LhsEval, LhsEval>
    inverseFormationVolumeFactorAndViscosity(const FluidState& fluidState, unsigned regionIdx) const
    {
        const LhsEval& T = decay<LhsEval>(fluidState.temperature(FluidState::waterPhaseIdx));
        const LhsEval& p = decay<LhsEval>(fluidState.pressure(FluidState::waterPhaseIdx));
        const LhsEval& saltConcentration
            = BlackOil::template getSaltConcentration_<FluidState, LhsEval>(fluidState, regionIdx);
        const auto salts = electrolytes(
            regionIdx, saltConcentration, BlackOil::getBrineComposition_(fluidState));
        return {density_(regionIdx, T, p, salts) / referenceDensity_[regionIdx],
                viscosityFactor_[regionIdx] * viscosity_(regionIdx, T, p, salts)};
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
    void computeFitRanges_(unsigned regionIdx);

    template <class Evaluation>
    static Evaluation clamp_(const Evaluation& x, Scalar lo, Scalar hi)
    {
        if (x < lo) {
            return Evaluation(lo);
        }
        if (x > hi) {
            return Evaluation(hi);
        }
        return x;
    }

    template <class Evaluation>
    Evaluation density_(unsigned regionIdx,
                        const Evaluation& temperature,
                        const Evaluation& pressure,
                        const ElectrolyteList<Evaluation>& salts) const
    {
        if (clampToFitRange_) {
            const auto& r = fitRange(regionIdx, FitProperty::Density);
            return densityAt_(clamp_(temperature, r.minTemperature, r.maxTemperature),
                              clamp_(pressure, Scalar {0.0}, r.maxPressure),
                              salts);
        }
        return densityAt_(temperature, pressure, salts);
    }

    template <class Evaluation>
    Evaluation viscosity_(unsigned regionIdx,
                          const Evaluation& temperature,
                          const Evaluation& pressure,
                          const ElectrolyteList<Evaluation>& salts) const
    {
        if (clampToFitRange_) {
            const auto& r = fitRange(regionIdx, FitProperty::Viscosity);
            return viscosityAt_(clamp_(temperature, r.minTemperature, r.maxTemperature),
                                clamp_(pressure, Scalar {0.0}, r.maxPressure),
                                salts);
        }
        return viscosityAt_(temperature, pressure, salts);
    }

    template <class Evaluation>
    Evaluation densityAt_(const Evaluation& temperature,
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
    Evaluation viscosityAt_(const Evaluation& temperature,
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
    std::vector<std::array<FitRange, 2>> fitRanges_ {};

    DensityModel densityModel_ {DensityModel::MaoDuanYoung};
    ViscosityModel viscosityModel_ {ViscosityModel::MaoDuan};
    EnthalpyModel enthalpyModel_ {EnthalpyModel::Driesner};
    PairingMethod pairing_ {PairingMethod::Proportional};
    bool laliberteCooperWater_ {false};
    bool enableSaltConcentration_ {false};
    bool clampToFitRange_ {false};
    Scalar referenceTemperature_ {288.71};
    Scalar referencePressure_ {1.01325e5};
};

} // namespace Opm

#endif
