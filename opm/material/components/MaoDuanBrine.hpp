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
 * \copydoc Opm::MaoDuanBrine
 */
#ifndef OPM_MAO_DUAN_BRINE_HPP
#define OPM_MAO_DUAN_BRINE_HPP

#include <opm/common/utility/SaltElectrolytes.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/components/DebyeHuckelSlope.hpp>
#include <opm/material/components/H2O.hpp>

#include <array>
#include <cstddef>
#include <optional>

namespace Opm
{

/*!
 * \brief Apparent molar volume of binary chloride solutions from Mao and Duan, J. Chem.
 * Thermodyn. 40, 2008.
 *
 * Pitzer form, eq 2, referenced to the volume at a fixed molality, eqs 8 and 10 to 12, so that
 *     V_phi(m) = V(m_r)/m_r - 1000/(rho_w m_r) + nu |z+ z-| A_V [h(I_m) - h(I_r)]
 *                + 2 nu+ nu- R T [B_V (m - m_r) + nu+ z+ C_V (m^2 - m_r^2)],
 * with h(I) = ln(1 + b I^0.5)/(2b), T in K, P in MPa, volumes in cm3 and rho_w the IAPWS-97
 * water density in g/cm3.
 */
template <class Scalar>
class MaoDuanBrine
{
    using Electrolyte = typename SaltElectrolytes<Scalar>::Electrolyte;

public:
    /// Ranges of Table 3: T [K], P [MPa] and molality [mol/kg]
    struct Range
    {
        Scalar maxT;
        Scalar maxP;
        Scalar maxMolality;
    };

    /// Whether Mao and Duan give a model for the electrolyte
    static constexpr bool covers(const Electrolyte electrolyte)
    {
        return index_(electrolyte).has_value();
    }

    static Range range(const Electrolyte electrolyte)
    {
        const auto& s = salts_[*index_(electrolyte)];
        return {s.maxT, s.maxP, s.maxMolality};
    }

    /// Ionic strength of one mole of the electrolyte in 1 kg of water
    static constexpr Scalar ionicStrengthPerMole(const Electrolyte electrolyte)
    {
        const auto& s = salts_[*index_(electrolyte)];
        return 0.5 * (s.nuCation * s.zCation * s.zCation + s.nuAnion * s.zAnion * s.zAnion);
    }

    /*!
     * IAPWS-97 water density, extended smoothly below the vapour pressure
     *
     * @param temperature Temperature [K]
     * @param pressure Pressure [Pa]
     * @return Density [kg/m3]
     */
    template <class Evaluation>
    static Evaluation waterDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        return H2O<Scalar>::liquidDensity(temperature, pressure, /*extrapolate=*/true);
    }

    /*!
     * Apparent molar volume of a chloride in its binary solution
     *
     * @param electrolyte LiCl, NaCl, KCl, MgCl2, CaCl2 or SrCl2
     * @param temperature Temperature [K]
     * @param pressure Pressure [Pa]
     * @param molality Molality of the electrolyte [mol/kg water]
     * @param rhoW Water density [kg/m3]
     * @return Apparent molar volume [m3/mol]
     */
    template <class Evaluation>
    static Evaluation apparentMolarVolume(const Electrolyte electrolyte,
                                          const Evaluation& temperature,
                                          const Evaluation& pressure,
                                          const Evaluation& molality,
                                          const Evaluation& rhoW)
    {
        const auto& s = salts_[*index_(electrolyte)];
        const auto& c = s.c;
        const Evaluation& T = temperature;
        const Evaluation P = pressure * 1e-6; // MPa
        const Evaluation T2 = T * T;
        const Evaluation T3 = T2 * T;
        const Evaluation low = 1.0 / (T - 227.0);
        const Evaluation high = 1.0 / (647.0 - T);

        // Eqs 10 to 12
        const Evaluation volumeRef = c[0] + c[1] * T + c[2] * T2 + c[3] * T3
            + P * (c[4] + c[5] * T + c[6] * T2 + c[7] * T3);
        const Evaluation bV = c[8] * low + c[9] + c[10] * T + c[11] * T2 + c[12] * high
            + P * (c[13] * low + c[14] + c[15] * T + c[16] * T2 + c[17] * high);
        const Evaluation cV = c[18] * low + c[19] + c[20] * T + c[21] * T2 + c[22] * high;

        const Scalar mr = s.refMolality;
        const Scalar nu = s.nuCation + s.nuAnion;
        const Scalar zz = s.zCation * s.zAnion;
        const Scalar ionic = ionicStrengthPerMole(electrolyte);
        const Evaluation aV = DebyeHuckelSlope<Scalar>::volumeSlope(temperature, pressure);
        const Evaluation& m = molality;

        const Evaluation vPhi = volumeRef / mr - 1e6 / (rhoW * mr)
            + nu * zz * aV * (h_(ionic * m) - h_(Evaluation(ionic * mr)))
            + 2.0 * s.nuCation * s.nuAnion * gasConstant_ * T
                * (bV * (m - mr) + s.nuCation * s.zCation * cV * (m * m - mr * mr));
        return vPhi * 1e-6;
    }

private:
    static constexpr Scalar gasConstant_ = 8.314472; // cm3 MPa/(mol K)
    static constexpr Scalar b_ = 1.2; // kg^0.5 mol^-0.5

    template <class Evaluation>
    static Evaluation h_(const Evaluation& ionicStrength)
    {
        return log(1.0 + b_ * sqrt(ionicStrength)) / (2.0 * b_);
    }

    struct Salt
    {
        Electrolyte electrolyte;
        Scalar nuCation;
        Scalar zCation;
        Scalar nuAnion;
        Scalar zAnion; // magnitude
        Scalar refMolality;
        Scalar maxT;
        Scalar maxP;
        Scalar maxMolality;
        std::array<Scalar, 23> c;
    };

    static constexpr std::optional<std::size_t> index_(const Electrolyte electrolyte)
    {
        for (std::size_t k = 0; k < salts_.size(); ++k) {
            if (salts_[k].electrolyte == electrolyte) {
                return k;
            }
        }
        return std::nullopt;
    }

    // Tables 2 and 3. BaCl2 is left out, as it is not one of our electrolytes.
    // clang-format off
    static constexpr std::array<Salt, 6> salts_ = {{
        {Electrolyte::LiCl, 1, 1, 1, 1, 10.0, 564.0, 40.0, 10.0,
         {1.17271480E+03, 1.40527916E-01, -5.53962649E-04, 1.72402126E-06,
          -1.57184556E+00, 8.89959461E-03, -1.52090064E-05, 0,
          7.33553879E-03, 4.06701494E-04, -2.38873863E-06, 3.94900863E-09, -3.56131664E-02,
          -6.18472877E-05, 3.00484214E-06, -1.02229075E-08, 0, 2.35592818E-04,
          -2.68117086E-04, -2.17228726E-05, 1.19732313E-07, -1.51104808E-10, 3.83403994E-04}},
        {Electrolyte::NaCl, 1, 1, 1, 1, 6.0, 573.0, 100.0, 6.0,
         {1.06607098E+03, -8.39622456E-03, 5.35429127E-04, 7.55373789E-07,
          -4.19512335E-01, 1.45082899E-03, -3.47807732E-06, 0,
          1.10913788E-02, 1.14498252E-03, -5.51181270E-06, 7.05483955E-09, -5.05734723E-02,
          -1.32747828E-04, 4.77261581E-06, -1.76888377E-08, 0, 6.40541237E-04,
          3.07698827E-04, -1.64042763E-04, 7.06784935E-07, -6.50338372E-10, -4.50906014E-04}},
        {Electrolyte::KCl, 1, 1, 1, 1, 6.0, 543.0, 50.0, 4.5,
         {2.90812061E+02, 6.54111195E+00, -1.61831978E-02, 1.46280384E-05,
          1.41397987E+01, -1.07266230E-01, 2.64506021E-04, -2.19789708E-07,
          3.02182158E-02, -2.15621394E-03, 9.24163206E-06, -1.10089434E-08, 2.87018859E-02,
          -6.73119697E-04, 1.68332473E-04, -7.99645640E-07, 1.11881560E-09, -6.59292385E-03,
          -2.02369103E-03, -1.70609099E-04, 1.00510108E-06, -1.86624642E-09, 1.91919166E-02}},
        {Electrolyte::MgCl2, 1, 2, 2, 1, 2.0, 543.0, 30.0, 3.0,
         {1.18880927E+03, -1.43194546E+00, 3.87973220E-03, -2.20330377E-06,
          6.38745038E+00, -5.51728055E-02, 1.50231562E-04, -1.35757912E-07,
          8.43627549E-03, 5.25365072E-03, -1.87204100E-05, 4.20263897E-08, -1.18062548E+00,
          6.07424747E-04, -1.20268210E-04, 5.23784551E-07, -8.23940319E-10, 9.75167613E-03,
          -4.92959181E-04, -2.73642775E-04, 5.42602386E-07, -1.95602825E-09, 1.00921935E-01}},
        {Electrolyte::CaCl2, 1, 2, 2, 1, 5.0, 523.0, 60.0, 6.0,
         {1.12080057E+03, -2.61669538E-01, 1.52042960E-03, -6.89131095E-07,
          -5.11802652E-01, 2.22234857E-03, -5.66464544E-06, 2.92950266E-09,
          2.43934633E-02, -1.42746873E-03, 7.35840529E-06, -9.43615480E-09, -5.18606814E-02,
          -6.16536928E-05, -1.04523561E-05, 4.52637296E-08, -1.05076158E-10, 2.31544709E-03,
          -1.09663211E-03, 1.90836111E-04, -9.25997994E-07, 1.54388261E-09, -1.29354832E-02}},
        {Electrolyte::SrCl2, 1, 2, 2, 1, 2.0, 473.0, 2.0, 2.0,
         {1.11894213E+03, -7.37321458E-01, 1.77908655E-03, 0,
          0, 0, 0, 0,
          2.21225680E-02, 6.62517291E-04, -2.37296050E-06, 0, 0,
          0, 0, 0, 0, 0,
          -4.21300430E-03, 0, 9.46738388E-08, 0, 0}},
    }};
    // clang-format on
};

} // namespace Opm

#endif // OPM_MAO_DUAN_BRINE_HPP
