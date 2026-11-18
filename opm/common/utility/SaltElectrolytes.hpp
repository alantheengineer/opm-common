/*
  Copyright 2026 NORCE.
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
#ifndef OPM_SALTELECTROLYTES_HPP
#define OPM_SALTELECTROLYTES_HPP

#include <opm/common/utility/ElectrolyteFractions.hpp>
#include <opm/common/utility/FixedCapacityList.hpp>
#include <opm/common/utility/SaltArray.hpp>

#include <opm/material/common/MathToolbox.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>

namespace Opm
{

/*!
 * Speciates a SaltArray of dissolved ions into the electrolytes (NaCl, CaCl2, ...)
 * needed by the Laliberte (viscosity) and Laliberte & Cooper (density) correlations
 * for multicomponent-salt brine.
 *
 * @tparam Scalar Value type used for molar masses
 */
template <class Scalar>
class SaltElectrolytes
{
public:
    /// Supported electrolytes for Laliberte viscosity and Laliberte-Cooper density
    // clang-format off
    enum class Electrolyte
    {
        CaCl2, KCl, K2SO4, MgCl2, MgSO4, NaCl, Na2SO4,
        LiCl, Li2SO4, SrCl2, NaBr, KBr, NaHCO3, Na2CO3, KHCO3, K2CO3,
        BOH3, // boric acid, a neutral solute
        NUM // Must be last: the number of electrolytes
    };
    // clang-format on

    static constexpr std::size_t numElectrolytes = static_cast<std::size_t>(Electrolyte::NUM);
    static_assert(numElectrolytes == numBrineElectrolytes);

    /// How ions are paired into electrolytes
    enum class PairingMethod {
        Proportional, ///< every cation with every anion, by share of charge
        Greedy ///< anions in turn take cations in a fixed order; the original method
    };

    /// Electrolytes and their mass fractions, held without heap allocation. Each
    /// cation-anion pair is visited once, so no electrolyte appears twice.
    template <class Evaluation>
    using ElectrolyteList = FixedCapacityList<std::pair<Electrolyte, Evaluation>, numElectrolytes>;

    /// An electrolyte and the moles of cation and anion in one mole of it. A neutral
    /// solute lists its own component as the cation and has no anion.
    struct Pairing
    {
        Electrolyte electrolyte;
        SaltIndex cation;
        Scalar cationsPer;
        SaltIndex anion;
        Scalar anionsPer;
    };

    /// Every supported electrolyte, in Electrolyte order
    // clang-format off
    static constexpr std::array<Pairing, numElectrolytes> pairings = {{
        {Electrolyte::CaCl2, SaltIndex::CA, 1, SaltIndex::CL, 2},
        {Electrolyte::KCl, SaltIndex::K, 1, SaltIndex::CL, 1},
        {Electrolyte::K2SO4, SaltIndex::K, 2, SaltIndex::SO4, 1},
        {Electrolyte::MgCl2, SaltIndex::MG, 1, SaltIndex::CL, 2},
        {Electrolyte::MgSO4, SaltIndex::MG, 1, SaltIndex::SO4, 1},
        {Electrolyte::NaCl, SaltIndex::NA, 1, SaltIndex::CL, 1},
        {Electrolyte::Na2SO4, SaltIndex::NA, 2, SaltIndex::SO4, 1},
        {Electrolyte::LiCl, SaltIndex::LI, 1, SaltIndex::CL, 1},
        {Electrolyte::Li2SO4, SaltIndex::LI, 2, SaltIndex::SO4, 1},
        {Electrolyte::SrCl2, SaltIndex::SR, 1, SaltIndex::CL, 2},
        {Electrolyte::NaBr, SaltIndex::NA, 1, SaltIndex::BR, 1},
        {Electrolyte::KBr, SaltIndex::K, 1, SaltIndex::BR, 1},
        {Electrolyte::NaHCO3, SaltIndex::NA, 1, SaltIndex::HCO3, 1},
        {Electrolyte::Na2CO3, SaltIndex::NA, 2, SaltIndex::CO3, 1},
        {Electrolyte::KHCO3, SaltIndex::K, 1, SaltIndex::HCO3, 1},
        {Electrolyte::K2CO3, SaltIndex::K, 2, SaltIndex::CO3, 1},
        {Electrolyte::BOH3, SaltIndex::B, 1, SaltIndex::B, 0},
    }};
    // clang-format on

    static constexpr bool isNeutral(const Pairing& p)
    {
        return p.anionsPer == 0;
    }

    static constexpr bool pairingsInElectrolyteOrder()
    {
        for (std::size_t i = 0; i < numElectrolytes; ++i) {
            if (pairings[i].electrolyte != static_cast<Electrolyte>(i)) {
                return false;
            }
        }
        return true;
    }

    /*!
     * Molar mass of electrolyte
     *
     * @param e Electrolyte
     * @return Molar mass [kg/mol]
     */
    static Scalar molarMass(const Electrolyte e)
    {
        static_assert(pairingsInElectrolyteOrder(), "pairings must follow Electrolyte order");
        const auto& p = pairings[static_cast<std::size_t>(e)];
        return p.cationsPer * saltMolarMass<Scalar>(p.cation)
            + p.anionsPer * saltMolarMass<Scalar>(p.anion);
    }

    /*!
     * Generate electrolytes from salt ions
     *
     * @param salinity Array of salt ions
     * @param method How ions are paired
     * @return Electrolytes and their mass fractions
     */
    template <class Evaluation>
    static ElectrolyteList<Evaluation>
    fromSaltComponents(const SaltArray<Evaluation, SaltMassFraction>& salinity,
                       const PairingMethod method = PairingMethod::Proportional)
    {
        auto molality = salinity.template convert_to<SaltMolality>();
        auto electrolytes = fromIonMolality_(molality, method);

        // Convert salt electrolytes from molality to mass fraction
        molalToMassFrac_(electrolytes);

        return electrolytes;
    }

    /*!
     * Molality of the ions left unpaired by fromSaltComponents. Those ions take no part in
     * the brine properties, so input checks warn when this is above zero.
     *
     * @param salinity Array of salt ions
     * @param method How ions are paired
     * @return Sum of unpaired ion molalities [mol/kg water]
     */
    template <class Evaluation>
    static Evaluation unpairedMolality(const SaltArray<Evaluation, SaltMassFraction>& salinity,
                                       const PairingMethod method = PairingMethod::Proportional)
    {
        auto molality = salinity.template convert_to<SaltMolality>();
        fromIonMolality_(molality, method);
        return molality.sum();
    }

    static constexpr Scalar molalTolerance = 1e-6; ///< tolerance for molality calculations

private:
    static constexpr std::size_t numCations_ = numSaltCations;
    static constexpr std::size_t numAnions_ = numSaltAnions;
    static constexpr int maxFitIterations_ = 200;
    static constexpr Scalar fitTolerance_ = 100 * std::numeric_limits<Scalar>::epsilon();

    /*!
     * Generate electrolytes (in molality) from salt ion molalities
     *
     * @param molalSalinity Molality of each salt ion; left holding the unpaired molalities
     * @param method How ions are paired
     * @return Electrolytes and their molality [mol electrolyte/kg water]
     */
    template <class Evaluation>
    static ElectrolyteList<Evaluation>
    fromIonMolality_(SaltArray<Evaluation, SaltMolality>& molalSalinity, const PairingMethod method)
    {
        auto electrolytes = method == PairingMethod::Greedy ? pairGreedy_(molalSalinity)
                                                            : pairProportional_(molalSalinity);

        // Neutral solutes need no partner
        for (const auto& p : pairings) {
            auto& molal = molalSalinity[p.cation];
            if (isNeutral(p) && molal >= molalTolerance) {
                electrolytes.push_back({p.electrolyte, molal});
                molal = 0.0;
            }
        }

        return electrolytes;
    }

    /*!
     * Pair ions anion by anion, strongest charge first, each taking cations in a fixed order
     *
     * @param molalSalinity Molality of each salt ion; left holding the unpaired molalities
     * @return Electrolytes and their molality [mol electrolyte/kg water]
     */
    template <class Evaluation>
    static ElectrolyteList<Evaluation>
    pairGreedy_(SaltArray<Evaluation, SaltMolality>& molalSalinity)
    {
        // Generate (valid) electrolytes from salt components
        const auto [cations, anions] = molalSalinity.cations_and_anions();
        ElectrolyteList<Evaluation> electrolytes;
        for (const auto& anionIdx : anions) {
            auto& molalAnion = molalSalinity[anionIdx];
            for (std::size_t i = 0; i < cations.size() && molalAnion >= molalTolerance; ++i) {
                const auto& cationIdx = cations[i];
                auto& molalCation = molalSalinity[cationIdx];
                if (molalCation < molalTolerance) {
                    continue;
                }

                // Generate salt electrolyte (in molal), and subtract used cation and anion
                // molalities. Note: Electrolytes are made from cation and anion pairs in
                // decreasing ion strength, and incompatible (or not implemented) pairs are skipped
                auto electrolyte =
                    electrolyteAndRemainingMolal_(cationIdx,
                                                  molalCation,
                                                  anionIdx,
                                                  molalAnion);
                if (electrolyte) {
                    electrolytes.push_back(*electrolyte);
                }
            }
        }

        return electrolytes;
    }

    /*!
     * Pair every cation with every anion in proportion to their charge, which is Young's rule
     * with equivalent fractions (Reilly and Wood, J. Phys. Chem. 73, 1969). Pairs without
     * Laliberte coefficients are left out, and the others are rescaled by iterative
     * proportional fitting so that each ion's charge is used once. Charge left over from an
     * unbalanced brine, or on ions with no supported partner, stays unpaired.
     *
     * @param molalSalinity Molality of each salt ion; left holding the unpaired molalities
     * @return Electrolytes and their molality [mol electrolyte/kg water]
     */
    template <class Evaluation>
    static ElectrolyteList<Evaluation>
    pairProportional_(SaltArray<Evaluation, SaltMolality>& molalSalinity)
    {
        const auto [cations, anions] = molalSalinity.cations_and_anions();
        const std::size_t nc = cations.size();
        const std::size_t na = anions.size();

        // Charge equivalents of each ion, and the supported pairs
        std::array<Evaluation, numCations_> cationEq {};
        std::array<Evaluation, numAnions_> anionEq {};
        std::array<std::array<int, numAnions_>, numCations_> pairIdx {};
        for (std::size_t i = 0; i < nc; ++i) {
            const auto& m = molalSalinity[cations[i]];
            const Scalar z = ionStrength(cations[i]);
            cationEq[i] = m >= molalTolerance ? m * z : Evaluation(0.0);
        }
        for (std::size_t j = 0; j < na; ++j) {
            const auto& m = molalSalinity[anions[j]];
            const Scalar z = -ionStrength(anions[j]);
            anionEq[j] = m >= molalTolerance ? m * z : Evaluation(0.0);
        }

        // Ions with no supported partner take no part
        std::array<bool, numCations_> cationPaired {};
        std::array<bool, numAnions_> anionPaired {};
        for (std::size_t i = 0; i < nc; ++i) {
            for (std::size_t j = 0; j < na; ++j) {
                const bool present = cationEq[i] > 0.0 && anionEq[j] > 0.0;
                pairIdx[i][j] = present ? pairingIndex_(cations[i], anions[j]) : -1;
                if (pairIdx[i][j] >= 0) {
                    cationPaired[i] = true;
                    anionPaired[j] = true;
                }
            }
        }
        Evaluation totalCationEq = 0.0;
        Evaluation totalAnionEq = 0.0;
        for (std::size_t i = 0; i < nc; ++i) {
            if (!cationPaired[i]) {
                cationEq[i] = 0.0;
            }
            totalCationEq += cationEq[i];
        }
        for (std::size_t j = 0; j < na; ++j) {
            if (!anionPaired[j]) {
                anionEq[j] = 0.0;
            }
            totalAnionEq += anionEq[j];
        }

        ElectrolyteList<Evaluation> electrolytes;
        if (totalCationEq < molalTolerance || totalAnionEq < molalTolerance) {
            return electrolytes;
        }

        // The larger side gives up its excess in proportion
        const Evaluation total = min(totalCationEq, totalAnionEq);
        std::array<Evaluation, numCations_> rowTarget {};
        std::array<Evaluation, numAnions_> colTarget {};
        for (std::size_t i = 0; i < nc; ++i) {
            rowTarget[i] = cationEq[i] * (total / totalCationEq);
        }
        for (std::size_t j = 0; j < na; ++j) {
            colTarget[j] = anionEq[j] * (total / totalAnionEq);
        }

        // Equivalents in each pair, starting from the proportional share
        std::array<std::array<Evaluation, numAnions_>, numCations_> eq {};
        for (std::size_t i = 0; i < nc; ++i) {
            for (std::size_t j = 0; j < na; ++j) {
                eq[i][j]
                    = pairIdx[i][j] >= 0 ? rowTarget[i] * colTarget[j] / total : Evaluation(0.0);
            }
        }
        fitToTargets_(eq, rowTarget, colTarget, nc, na);

        for (std::size_t i = 0; i < nc; ++i) {
            for (std::size_t j = 0; j < na; ++j) {
                if (pairIdx[i][j] < 0 || eq[i][j] <= 0.0) {
                    continue;
                }
                const auto& p = pairings[static_cast<std::size_t>(pairIdx[i][j])];
                const Evaluation m = eq[i][j] / (p.cationsPer * ionStrength(p.cation));
                molalSalinity[p.cation] -= m * p.cationsPer;
                molalSalinity[p.anion] -= m * p.anionsPer;
                electrolytes.push_back({p.electrolyte, m});
            }
        }
        return electrolytes;
    }

    /*!
     * Scale rows and columns in turn until the row and column sums meet their targets. Rows
     * and columns that are all zero are left alone. If the targets cannot all be met, the
     * columns are capped at their targets so no ion is used beyond what is present.
     */
    template <class Evaluation>
    static void fitToTargets_(std::array<std::array<Evaluation, numAnions_>, numCations_>& eq,
                              const std::array<Evaluation, numCations_>& rowTarget,
                              const std::array<Evaluation, numAnions_>& colTarget,
                              const std::size_t nc,
                              const std::size_t na)
    {
        const auto scaleColumns = [&](bool capOnly) {
            for (std::size_t j = 0; j < na; ++j) {
                Evaluation sum = 0.0;
                for (std::size_t i = 0; i < nc; ++i) {
                    sum += eq[i][j];
                }
                if (sum > 0.0 && (!capOnly || sum > colTarget[j])) {
                    const Evaluation f = colTarget[j] / sum;
                    for (std::size_t i = 0; i < nc; ++i) {
                        eq[i][j] *= f;
                    }
                }
            }
        };

        for (int iter = 0; iter < maxFitIterations_; ++iter) {
            Scalar worst = 0.0;
            for (std::size_t i = 0; i < nc; ++i) {
                Evaluation sum = 0.0;
                for (std::size_t j = 0; j < na; ++j) {
                    sum += eq[i][j];
                }
                if (sum > 0.0) {
                    const Evaluation f = rowTarget[i] / sum;
                    worst = std::max(worst, static_cast<Scalar>(std::abs(getValue(f) - 1.0)));
                    for (std::size_t j = 0; j < na; ++j) {
                        eq[i][j] *= f;
                    }
                }
            }
            if (worst < fitTolerance_) {
                return;
            }
            scaleColumns(false);
        }
        scaleColumns(true);
    }

    /// Index into pairings of the electrolyte made of a cation and an anion, or -1
    static constexpr int pairingIndex_(const SaltIndex cation, const SaltIndex anion)
    {
        for (std::size_t k = 0; k < numElectrolytes; ++k) {
            const auto& p = pairings[k];
            if (!isNeutral(p) && p.cation == cation && p.anion == anion) {
                return static_cast<int>(k);
            }
        }
        return -1;
    }

    /*!
     * @brief Calculate electrolyte from cation and anion, and adjust molality of cation and
     * anion after generating the electrolyte
     *
     * The electrolyte is made from the dissociation (or chemical equilibrium) equation
     * A_pB_q = pA+ + qB-. Hence, to make 1 molal of A_pB_q (electrolyte), p molal of A+ (cation)
     * and q molal of B- (anion) are needed. In this function we make as many molal electrolyte
     * as the maximum amount of cations or anions allow, based on the dissociation equation.
     *
     * @param cationIdx Salt index of cation
     * @param cation Molality of cation [mol cation/kg water]
     * @param anionIdx Salt index of anion
     * @param anion Molality of anion [mol anion/kg water]
     * @return Electrolyte and its molality [mol electrolyte/kg water], or std::nullopt if the
     *         cation/anion pair does not form a supported electrolyte (cation and anion
     *         molalities are left untouched in that case)
     */
    template <class Evaluation>
    static std::optional<std::pair<Electrolyte, Evaluation> >
    electrolyteAndRemainingMolal_(const SaltIndex cationIdx,
                                  Evaluation& cation,
                                  const SaltIndex anionIdx,
                                  Evaluation& anion)
    {
        for (const auto& p : pairings) {
            if (isNeutral(p) || p.cation != cationIdx || p.anion != anionIdx) {
                continue;
            }
            const Evaluation electrolyteMolal = min(cation / p.cationsPer, anion / p.anionsPer);
            cation -= electrolyteMolal * p.cationsPer;
            anion -= electrolyteMolal * p.anionsPer;
            return std::pair<Electrolyte, Evaluation> {p.electrolyte, electrolyteMolal};
        }
        return std::nullopt;
    }

    /*!
     * Convert electrolyte molality to mass fraction
     *
     * @param electrolyte Electrolytes and their molality
     */
    template <class Evaluation>
    static void molalToMassFrac_(ElectrolyteList<Evaluation>& electrolyte)
    {
        Evaluation sum = 1.0;
        for (auto& salt : electrolyte) {
            auto mMsalt = molarMass(salt.first);
            salt.second *= mMsalt;
            sum += salt.second;
        }
        for (auto& salt : electrolyte) {
            salt.second /= sum;
        }
    }
};

} // namespace Opm

#endif // OPM_SALTELECTROLYTES_HPP
