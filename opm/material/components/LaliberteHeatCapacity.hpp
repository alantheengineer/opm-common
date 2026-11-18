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
 * \copydoc Opm::LaliberteHeatCapacity
 */
#ifndef OPM_LALIBERTE_HEAT_CAPACITY_HPP
#define OPM_LALIBERTE_HEAT_CAPACITY_HPP

#include <opm/common/utility/SaltElectrolytes.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

namespace Opm
{

/*!
 * \brief Solute and water heat capacities from Laliberte, J. Chem. Eng. Data 54, 2009.
 *
 * The solution heat capacity is (1 - W) Cpw + sum(w_i Cp_i) (eq 10), with
 *     Cp_i = a1 exp(alpha) + a5 W^a6,  alpha = a2 t + a3 exp(0.01 t) + a4 W  (eqs 11 and 12),
 * where W is the total solute mass fraction and t is in degrees C.
 */
template <class Scalar>
class LaliberteHeatCapacity
{
    using Electrolyte = typename SaltElectrolytes<Scalar>::Electrolyte;
    static constexpr std::size_t numElectrolytes = SaltElectrolytes<Scalar>::numElectrolytes;

public:
    /*!
     * Heat capacity of water from Table 1 of the paper, with its quadratic interpolation
     * (eq 9) through the node at or below t and the next two. Values come from IAPWS-95 at
     * 0.1 MPa to 95 C and at saturation above.
     *
     * @param temperature Temperature [K]
     * @return Heat capacity [J/(kg K)]
     */
    template <class Evaluation>
    static Evaluation water(const Evaluation& temperature)
    {
        const Evaluation tempC = temperature - 273.15;
        const Scalar pos = (scalarValue(tempC) - waterFirstC_) / waterStepC_;
        const int last = static_cast<int>(cpWater_.size()) - 3;
        const int k = std::clamp(static_cast<int>(std::floor(pos)), 0, last);

        const Evaluation x = (tempC - (waterFirstC_ + k * waterStepC_)) / waterStepC_;
        const Scalar c1 = cpWater_[k];
        const Scalar c2 = cpWater_[k + 1];
        const Scalar c3 = cpWater_[k + 2];
        return 1e3 * (c1 + (c2 - c1) * x + 0.5 * (c3 - 2 * c2 + c1) * x * (x - 1));
    }

    /*!
     * Apparent heat capacity of one solute (eqs 11 and 12)
     *
     * @param electrolyte Electrolyte, not boric acid
     * @param tempC Temperature [C]
     * @param sumW Total solute mass fraction [-]
     * @return Heat capacity [J/(kg K)]
     */
    template <class Evaluation>
    static Evaluation
    solute(const Electrolyte electrolyte, const Evaluation& tempC, const Evaluation& sumW)
    {
        const auto& a = coefficients_[static_cast<std::size_t>(electrolyte)];
        return 1e3
            * (a[0] * exp(a[1] * tempC + a[2] * exp(0.01 * tempC) + a[3] * sumW)
               + a[4] * pow(sumW, a[5]));
    }

    /*!
     * Integral of the solute heat capacity from 0 C to t, at fixed composition
     *
     * @param electrolyte Electrolyte, not boric acid
     * @param tempC Temperature [C]
     * @param sumW Total solute mass fraction [-]
     * @return Enthalpy relative to 0 C [J/kg]
     */
    template <class Evaluation>
    static Evaluation
    soluteEnthalpy(const Electrolyte electrolyte, const Evaluation& tempC, const Evaluation& sumW)
    {
        const std::size_t e = static_cast<std::size_t>(electrolyte);
        const auto& a = coefficients_[e];
        return 1e3
            * (a[0] * exp(a[3] * sumW) * expIntegral_(e, tempC) + a[4] * pow(sumW, a[5]) * tempC);
    }

private:
    // Water heat capacity, Table 1 [kJ/(kg K)], from -15 C in steps of 5 C
    static constexpr Scalar waterFirstC_ = -15.0;
    static constexpr Scalar waterStepC_ = 5.0;
    // clang-format off
    static constexpr std::array<Scalar, 32> cpWater_ = {
        4.294027456829579, 4.256878595310196, 4.233579005755254, 4.21944, 4.20495, 4.19545,
        4.1891, 4.1848, 4.1819, 4.18002, 4.17895, 4.178856, 4.17877, 4.17956, 4.18089, 4.18277,
        4.18517, 4.1881, 4.19155, 4.19552, 4.20001, 4.20502, 4.21057, 4.21664, 4.22323, 4.23036,
        4.23807, 4.24637, 4.25528, 4.26484, 4.27508, 4.28604,
    };
    // clang-format on

    static constexpr Scalar unused_ = std::numeric_limits<Scalar>::quiet_NaN();

    // Coefficients a1 to a6 from the 2009 supplement, in Electrolyte order
    // clang-format off
    static constexpr std::array<std::array<Scalar, 6>, numElectrolytes> coefficients_ = {{
        // CaCl2
        {-1.3892271378463972, -0.014249134161856351, 0.5782474297490663, -0.7853394719779166,
         4.398953416292243, 1.1268559362341135},
        // KCl
        {-1.6395291777159238, -0.007794617046554878, 0.32842265724913666, -3.768724412355345,
         0.34049805034804986, -0.1479847834204727},
        // K2SO4
        {-0.24004643791671176, -0.05282385665347616, 2.3697421382867554, -9.727161136383726,
         -20.064280098564193, 2.1901955033828595},
        // MgCl2
        {-1.3460686403323334, -0.018554256655827786, 0.8465222110093862, -3.54621764242775,
         6.642171961158046, 3.3808284305138807},
        // MgSO4
        {-7.629919298283411, 0.01084735529840594, -1.4861671624644168, -1.8978254620975095,
         1.7200716872247217, 0.4567747446135715},
        // NaCl
        {-0.06935596689933224, -0.0782134167486952, 3.847984794086355, -11.276210924707238,
         8.731876985426716, 1.8124593047275528},
        // Na2SO4
        {-0.1289400598506014, -0.07910078502732856, 3.2220490671730038, -15.69342507488851,
         5.805447171835829, 1.1836256711409199},
        // LiCl
        {-0.11383628170471641, -0.05840807773194118, 2.7078733038320513, -6.803845972193638,
         -0.15302743092275495, -0.23409509237124265},
        // Li2SO4
        {-50.99777193898149, -3.7483598249095302, 13.195983918030464, 21.393023454955074,
         -0.042807714268368095, -0.8711840587005555},
        // SrCl2
        {-0.6775629429624657, -0.024565183831394628, 1.136773269847769, -1.9732236131755434,
         1.6619925635594976, 1.0240418255796286},
        // NaBr
        {-0.11601939402742224, -0.0583057782055604, 2.795111329461777, -8.569275106840628,
         3.0663306694850507, 2.216703450020728},
        // KBr
        {-0.4701916974368641, -0.03861897620442665, 1.263118612440269, -0.3167742239720588,
         4.526718491147356, 1.2915367831693083},
        // NaHCO3
        {-0.030688580780269314, -0.02716693132217614, 1.2634502296851688, -0.3524152511572335,
         4.158479413300207, 0.6822033111940248},
        // Na2CO3
        {-0.5368152124011092, -0.02716693132217614, 1.2634502296851688, -9.170500724149688,
         5.655502384934204, 0.8280134584293479},
        // KHCO3
        {-0.7940179046250507, -0.02716693132217614, 1.2634502296851688, -2.429610557068276,
         1.8105200307805134, 0.010342092637275445},
        // K2CO3
        {-2.2145294382324208, -0.03077332449455719, 0.337197867325901, -6.97645705984794,
         -0.23847045479667833, 0.46122724906023577},
        // B(OH)3, not used
        {unused_, unused_, unused_, unused_, unused_, unused_},
    }};

    // Six-point Gauss-Legendre rule on [-1, 1]
    static constexpr std::array<Scalar, 6> glNodes_ = {
        -0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
         0.2386191860831969,  0.6612093864662645,  0.9324695142031521,
    };
    static constexpr std::array<Scalar, 6> glWeights_ = {
        0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
        0.4679139345726910, 0.3607615730481386, 0.1713244923791704,
    };
    // clang-format on

    // Nodes of the integral table, every 1 C
    static constexpr int tableFirstC_ = -15;
    static constexpr int tableLastC_ = 250;
    static constexpr std::size_t tableSize_ = tableLastC_ - tableFirstC_ + 1;
    using Table = std::array<std::array<Scalar, tableSize_>, numElectrolytes>;

    // exp(alpha) without the composition term
    template <class Evaluation>
    static Evaluation expAlpha_(const std::array<Scalar, 6>& a, const Evaluation& tempC)
    {
        return exp(a[1] * tempC + a[2] * exp(0.01 * tempC));
    }

    template <class Evaluation>
    static Evaluation
    integrate_(const std::array<Scalar, 6>& a, const Evaluation& from, const Evaluation& to)
    {
        const Evaluation half = 0.5 * (to - from);
        const Evaluation mid = 0.5 * (to + from);
        Evaluation sum = 0.0;
        for (std::size_t k = 0; k < glNodes_.size(); ++k) {
            sum += glWeights_[k] * expAlpha_(a, mid + half * glNodes_[k]);
        }
        return half * sum;
    }

    /*!
     * Integral of exp(a2 s + a3 exp(0.01 s)) from 0 to t. The integrand can change by a
     * factor of 40 per degree, so values at 1 C steps are tabulated once and the step
     * holding t is integrated directly. Outside the table the end step is extended.
     */
    template <class Evaluation>
    static Evaluation expIntegral_(const std::size_t e, const Evaluation& tempC)
    {
        const auto& table = integralTable_();
        const Scalar pos = scalarValue(tempC) - tableFirstC_;
        const int k
            = std::clamp(static_cast<int>(std::floor(pos)), 0, static_cast<int>(tableSize_) - 1);
        const Scalar nodeC = static_cast<Scalar>(tableFirstC_ + k);
        return table[e][k] + integrate_(coefficients_[e], Evaluation(nodeC), tempC);
    }

    static const Table& integralTable_()
    {
        static const Table table = buildIntegralTable_();
        return table;
    }

    static Table buildIntegralTable_()
    {
        Table table {};
        constexpr int zero = -tableFirstC_;
        for (std::size_t e = 0; e < numElectrolytes; ++e) {
            if (SaltElectrolytes<Scalar>::isNeutral(SaltElectrolytes<Scalar>::pairings[e])) {
                continue;
            }
            const auto& a = coefficients_[e];
            for (int k = zero + 1; k < static_cast<int>(tableSize_); ++k) {
                const Scalar t = static_cast<Scalar>(tableFirstC_ + k);
                table[e][k] = table[e][k - 1] + integrate_(a, t - 1, t);
            }
            for (int k = zero - 1; k >= 0; --k) {
                const Scalar t = static_cast<Scalar>(tableFirstC_ + k);
                table[e][k] = table[e][k + 1] - integrate_(a, t, t + 1);
            }
        }
        return table;
    }
};

} // namespace Opm

#endif // OPM_LALIBERTE_HEAT_CAPACITY_HPP
