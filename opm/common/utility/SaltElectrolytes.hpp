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

#include <opm/common/utility/FixedCapacityList.hpp>
#include <opm/common/utility/SaltArray.hpp>

#include <opm/material/common/MathToolbox.hpp>

#include <array>
#include <cstddef>
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
    enum class Electrolyte
    {
        CaCl2, KCl, K2SO4, MgCl2, MgSO4, NaCl, Na2SO4,
        LiCl, Li2SO4, SrCl2, NaBr, KBr, NaHCO3, Na2CO3, KHCO3, K2CO3,
        NUM // Must be last: the number of electrolytes
    };

    static constexpr std::size_t numElectrolytes = static_cast<std::size_t>(Electrolyte::NUM);

    /// Electrolytes and their mass fractions, held without heap allocation. Each
    /// cation-anion pair is visited once, so no electrolyte appears twice.
    template <class Evaluation>
    using ElectrolyteList = FixedCapacityList<std::pair<Electrolyte, Evaluation>, numElectrolytes>;

    /// An electrolyte and the moles of cation and anion in one mole of it
    struct Pairing
    {
        Electrolyte electrolyte;
        SaltIndex cation;
        Scalar cationsPer;
        SaltIndex anion;
        Scalar anionsPer;
    };

    /// Every supported cation-anion pair, in Electrolyte order
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
    }};

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
     * @return Electrolytes and their mass fractions
     */
    template <class Evaluation>
    static ElectrolyteList<Evaluation>
    fromSaltComponents(const SaltArray<Evaluation, SaltMassFraction>& salinity)
    {
        auto molality = salinity.template convert_to<SaltMolality>();
        auto electrolytes = fromIonMolality_(molality);

        // Convert salt electrolytes from molality to mass fraction
        molalToMassFrac_(electrolytes);

        return electrolytes;
    }

    /*!
     * Molality of the ions left unpaired by fromSaltComponents. Those ions take no part in
     * the brine properties, so input checks warn when this is above zero.
     *
     * @param salinity Array of salt ions
     * @return Sum of unpaired ion molalities [mol/kg water]
     */
    template <class Evaluation>
    static Evaluation unpairedMolality(const SaltArray<Evaluation, SaltMassFraction>& salinity)
    {
        auto molality = salinity.template convert_to<SaltMolality>();
        fromIonMolality_(molality);
        return molality.sum();
    }

    static constexpr Scalar molalTolerance = 1e-6; ///< tolerance for molality calculations

private:

    /*!
     * Generate electrolytes (in molality) from salt ion molalities
     *
     * @param molalSalinity Molality of each salt ion; left holding the unpaired molalities
     * @return Electrolytes and their molality [mol electrolyte/kg water]
     */
    template <class Evaluation>
    static ElectrolyteList<Evaluation>
    fromIonMolality_(SaltArray<Evaluation, SaltMolality>& molalSalinity)
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
            if (p.cation != cationIdx || p.anion != anionIdx) {
                continue;
            }
            const Evaluation electrolyteMolal = min(cation / p.cationsPer, anion / p.anionsPer);
            cation -= electrolyteMolal * p.cationsPer;
            anion -= electrolyteMolal * p.anionsPer;
            return std::pair<Electrolyte, Evaluation>{p.electrolyte, electrolyteMolal};
        }
        return std::nullopt;
    }

    /*!
     * Convert electrolyte molality to mass fraction
     *
     * @param electrolyte Electrolytes and their molality
     */
    template <class Evaluation>
    static void
    molalToMassFrac_(ElectrolyteList<Evaluation>& electrolyte)
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
