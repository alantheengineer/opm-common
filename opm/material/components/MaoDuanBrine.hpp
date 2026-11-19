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
 * \brief Binary chloride solutions after Mao and Duan: apparent molar volume from J. Chem.
 * Thermodyn. 40, 2008, and viscosity from Int. J. Thermophys. 30, 2009.
 *
 * Volume:
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
    /// Fitted range: T [K], P [MPa] and molality [mol/kg]
    struct Range
    {
        Scalar minT;
        Scalar maxT;
        Scalar maxP;
        Scalar maxMolality;
    };

    /// Whether Mao and Duan give a model for the electrolyte
    static constexpr bool covers(const Electrolyte electrolyte)
    {
        return index_(electrolyte).has_value();
    }

    /// Fitted range of the volume model, from the paper's abstract and Section 3
    static Range range(const Electrolyte electrolyte)
    {
        const auto& s = salts_[*index_(electrolyte)];
        return {s.minT, s.maxT, s.maxP, s.maxMolality};
    }

    /// Fitted range of the viscosity model
    static Range viscosityRange(const Electrolyte electrolyte)
    {
        return viscosity_[*viscosityIndex_(electrolyte)].range;
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

    /// Whether a relative viscosity model is given for the electrolyte
    static constexpr bool coversViscosity(const Electrolyte electrolyte)
    {
        return viscosityIndex_(electrolyte).has_value();
    }

    /*!
     * Viscosity of water, Mao and Duan (2009), eq 6 and Table 2, with IAPWS-97 density
     *
     * @param temperature Temperature [K]
     * @param pressure Pressure [Pa]
     * @return Viscosity [Pa s]
     */
    template <class Evaluation>
    static Evaluation waterViscosity(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation rho = waterDensity(temperature, pressure) * 1e-3; // g/cm3
        const Evaluation& T = temperature;
        const Evaluation lnMu = d_[0] / (T * T) + d_[1] / T + d_[2] + d_[3] * T + d_[4] * T * T
            + rho * (d_[5] / (T * T) + d_[6] / T + d_[7] + d_[8] * T + d_[9] * T * T);
        return exp(lnMu);
    }

    /*!
     * Logarithm of the relative viscosity of a binary solution. LiCl, NaCl and KCl follow
     * Mao and Duan (2009), eqs 1 to 5 and Table 3. CaCl2 uses the same form, fitted for this
     * code to Abdulagatov and Azizov, Fluid Phase Equilib. 240, 2006 (293 to 575 K, 0.1 to
     * 60 MPa, to 2 mol/kg) and to the data compiled by Laliberte, J. Chem. Eng. Data 54,
     * 2009 (0 to 55 C, to 6 mol/kg); it reproduces them to 0.6% and 1.6% on average.
     *
     * @param electrolyte LiCl, NaCl, KCl or CaCl2
     * @param temperature Temperature [K]
     * @param molality Molality of the electrolyte [mol/kg water]
     * @return ln(viscosity / water viscosity) [-]
     */
    template <class Evaluation>
    static Evaluation lnRelativeViscosity(const Electrolyte electrolyte,
                                          const Evaluation& temperature,
                                          const Evaluation& molality)
    {
        const auto& v = viscosity_[*viscosityIndex_(electrolyte)];
        const Evaluation t = temperature / v.tempScale;
        const Evaluation a = v.a[0] + t * (v.a[1] + t * (v.a[2] + t * v.a[3]));
        const Evaluation b = v.b[0] + t * (v.b[1] + t * v.b[2]);
        const Evaluation c = v.c[0] + t * v.c[1];
        const Evaluation& m = molality;
        return m * (a + m * (b + m * c));
    }

private:
    static constexpr Scalar gasConstant_ = 8.314472; // cm3 MPa/(mol K)

    // Mao and Duan (2009), Table 2: ln(mu) = sum d_i T^(i-3) + rho sum d_i T^(i-8)
    static constexpr std::array<Scalar, 10> d_ = {0.28853170e7,
                                                  -0.11072577e5,
                                                  -0.90834095e1,
                                                  0.30925651e-1,
                                                  -0.27407100e-4,
                                                  -0.19283851e7,
                                                  0.56216046e4,
                                                  0.13827250e2,
                                                  -0.47609523e-1,
                                                  0.35545041e-4};

    // ln(eta_r) = A m + B m^2 + C m^3, with A, B and C polynomials in T / tempScale
    struct ViscositySet
    {
        Electrolyte electrolyte;
        Range range;
        Scalar tempScale; // K
        std::array<Scalar, 4> a;
        std::array<Scalar, 3> b;
        std::array<Scalar, 2> c;
    };

    static constexpr std::optional<std::size_t> viscosityIndex_(const Electrolyte electrolyte)
    {
        for (std::size_t k = 0; k < viscosity_.size(); ++k) {
            if (viscosity_[k].electrolyte == electrolyte) {
                return k;
            }
        }
        return std::nullopt;
    }

    // clang-format off
    static constexpr std::array<ViscositySet, 4> viscosity_ = {{
        // Mao and Duan (2009), Table 3, T in K; valid ranges from its Section 4
        {Electrolyte::LiCl, {273.15, 623.0, 100.0, 16.7}, 1.0,
         {0.62204136e-2, 0.54436974e-3, -0.40443190e-6, 0},
         {0.14987325e-1, -0.66617390e-4, 0.52113332e-7},
         {0.12101624e-5, 0.17772678e-6}},
        {Electrolyte::NaCl, {273.15, 623.0, 100.0, 6.0}, 1.0,
         {-0.21319213, 0.13651589e-2, -0.12191756e-5, 0},
         {0.69161945e-1, -0.27292263e-3, 0.20852448e-6},
         {-0.25988855e-2, 0.77989227e-5}},
        {Electrolyte::KCl, {273.15, 623.0, 100.0, 4.5}, 1.0,
         {-0.42122934, 0.18286059e-2, -0.13603098e-5, 0},
         {0.11380205e-1, 0.47541391e-5, -0.99280575e-7},
         {0, 0}},
        // Fitted for this code, T / 300 K; range of the data it was fitted to
        {Electrolyte::CaCl2, {273.15, 575.0, 60.0, 6.0}, 300.0,
         {-1.3438892472e+00, 3.0573535084e+00, -1.7606044540e+00, 3.3193259481e-01},
         {2.4475740832e-01, -3.5265818370e-01, 1.1381394527e-01},
         {1.2693964694e-03, 0}},
    }};
    // clang-format on

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
        Scalar minT;
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
        {Electrolyte::LiCl, 1, 1, 1, 1, 10.0, 273.15, 564.0, 40.0, 10.0,
         {1.17271480E+03, 1.40527916E-01, -5.53962649E-04, 1.72402126E-06,
          -1.57184556E+00, 8.89959461E-03, -1.52090064E-05, 0,
          7.33553879E-03, 4.06701494E-04, -2.38873863E-06, 3.94900863E-09, -3.56131664E-02,
          -6.18472877E-05, 3.00484214E-06, -1.02229075E-08, 0, 2.35592818E-04,
          -2.68117086E-04, -2.17228726E-05, 1.19732313E-07, -1.51104808E-10, 3.83403994E-04}},
        {Electrolyte::NaCl, 1, 1, 1, 1, 6.0, 273.15, 573.0, 100.0, 6.0,
         {1.06607098E+03, -8.39622456E-03, 5.35429127E-04, 7.55373789E-07,
          -4.19512335E-01, 1.45082899E-03, -3.47807732E-06, 0,
          1.10913788E-02, 1.14498252E-03, -5.51181270E-06, 7.05483955E-09, -5.05734723E-02,
          -1.32747828E-04, 4.77261581E-06, -1.76888377E-08, 0, 6.40541237E-04,
          3.07698827E-04, -1.64042763E-04, 7.06784935E-07, -6.50338372E-10, -4.50906014E-04}},
        {Electrolyte::KCl, 1, 1, 1, 1, 6.0, 273.15, 543.0, 50.0, 4.5,
         {2.90812061E+02, 6.54111195E+00, -1.61831978E-02, 1.46280384E-05,
          1.41397987E+01, -1.07266230E-01, 2.64506021E-04, -2.19789708E-07,
          3.02182158E-02, -2.15621394E-03, 9.24163206E-06, -1.10089434E-08, 2.87018859E-02,
          -6.73119697E-04, 1.68332473E-04, -7.99645640E-07, 1.11881560E-09, -6.59292385E-03,
          -2.02369103E-03, -1.70609099E-04, 1.00510108E-06, -1.86624642E-09, 1.91919166E-02}},
        {Electrolyte::MgCl2, 1, 2, 2, 1, 2.0, 273.15, 543.0, 30.0, 3.0,
         {1.18880927E+03, -1.43194546E+00, 3.87973220E-03, -2.20330377E-06,
          6.38745038E+00, -5.51728055E-02, 1.50231562E-04, -1.35757912E-07,
          8.43627549E-03, 5.25365072E-03, -1.87204100E-05, 4.20263897E-08, -1.18062548E+00,
          6.07424747E-04, -1.20268210E-04, 5.23784551E-07, -8.23940319E-10, 9.75167613E-03,
          -4.92959181E-04, -2.73642775E-04, 5.42602386E-07, -1.95602825E-09, 1.00921935E-01}},
        {Electrolyte::CaCl2, 1, 2, 2, 1, 5.0, 273.15, 523.0, 60.0, 6.0,
         {1.12080057E+03, -2.61669538E-01, 1.52042960E-03, -6.89131095E-07,
          -5.11802652E-01, 2.22234857E-03, -5.66464544E-06, 2.92950266E-09,
          2.43934633E-02, -1.42746873E-03, 7.35840529E-06, -9.43615480E-09, -5.18606814E-02,
          -6.16536928E-05, -1.04523561E-05, 4.52637296E-08, -1.05076158E-10, 2.31544709E-03,
          -1.09663211E-03, 1.90836111E-04, -9.25997994E-07, 1.54388261E-09, -1.29354832E-02}},
        {Electrolyte::SrCl2, 1, 2, 2, 1, 2.0, 298.15, 473.0, 2.0, 2.0,
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
