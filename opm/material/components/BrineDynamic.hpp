// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::BrineDynamic
 */
#ifndef OPM_BRINEDYNAMIC_HPP
#define OPM_BRINEDYNAMIC_HPP

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/common/utility/SaltArray.hpp>
#include <opm/common/utility/SaltElectrolytes.hpp>
#include <opm/material/components/Component.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <numeric>
#include <span>
#include <string_view>

namespace Opm {

/*!
 * \ingroup Components
 *
 * \brief A class for the brine fluid properties.
 *
 * \tparam Scalar The type used for scalar values
 * \tparam H2O Static polymorphism: the Brine class can access all properties of the H2O class
 */
template <class Scalar, class H2O>
class BrineDynamic : public Component<Scalar, BrineDynamic<Scalar, H2O> >
{
public:
    using PairingMethod = typename SaltElectrolytes<Scalar>::PairingMethod;

    /*!
     * \copydoc Component::name
     */
    static std::string_view name()
    { return "Brine"; }

    /*!
     * \copydoc H2O::gasIsIdeal
     */
    OPM_HOST_DEVICE static bool gasIsIdeal()
    { return H2O::gasIsIdeal(); }

    /*!
     * \copydoc H2O::gasIsCompressible
     */
    OPM_HOST_DEVICE static bool gasIsCompressible()
    { return H2O::gasIsCompressible(); }

    /*!
     * \copydoc H2O::liquidIsCompressible
     */
    OPM_HOST_DEVICE static bool liquidIsCompressible()
    { return H2O::liquidIsCompressible(); }

    /*!
     * \copydoc Component::molarMass
     *
     * This assumes that the salt is pure NaCl.
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation molarMass(const Evaluation& salinity)
    {
        const Scalar M1 = H2O::molarMass();
        const Evaluation X2 = salinity; // mass fraction of salt in brine
        return M1*mM_salt()/(mM_salt() + X2*(M1 - mM_salt()));
    }

    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation
    invAvgMolarMassFromMassFrac(const SaltArray<Evaluation, SaltMassFraction>& salinity)
    {
        const Scalar mH2O = H2O::molarMass();
        Evaluation s = 1.0 / mH2O;
        for (std::size_t i = 0; i < salinity.size(); ++i) {
            auto sIdx = static_cast<SaltIndex>(i);
            auto mIon = saltMolarMass<Scalar>(sIdx);
            s += salinity[sIdx] * ((mH2O - mIon) / (mH2O * mIon));
        }
        return s;
    }

    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation
    avgMolarMassFromMoleFrac(const SaltArray<Evaluation, SaltMoleFraction>& salinity)
    {
        const Scalar mH2O = H2O::molarMass();
        Evaluation s = mH2O;
        for (std::size_t i = 0; i < salinity.size(); ++i) {
            auto sIdx = static_cast<SaltIndex>(i);
            auto mIon = saltMolarMass<Scalar>(sIdx);
            s += salinity[sIdx] * (mIon - mH2O);
        }
        return s;
    }

    /*!
     * \copydoc H2O::criticalTemperature
     */
    OPM_HOST_DEVICE static Scalar criticalTemperature()
    { return H2O::criticalTemperature(); /* [K] */ }

    /*!
     * \copydoc H2O::criticalPressure
     */
    OPM_HOST_DEVICE static Scalar criticalPressure()
    { return H2O::criticalPressure(); /* [N/m^2] */ }

    /*!
     * \copydoc H2O::criticalVolume
     */
    OPM_HOST_DEVICE static Scalar criticalVolume()
    { return H2O::criticalVolume(); /* [m3/kmol] */ }

    /*!
     * \copydoc H20::acentricFactor
     */
    OPM_HOST_DEVICE static Scalar acentricFactor()
    { return H2O::acentricFactor(); }

    /*!
     * \copydoc H2O::tripleTemperature
     */
    OPM_HOST_DEVICE static Scalar tripleTemperature()
    { return H2O::tripleTemperature(); /* [K] */ }

    /*!
     * \copydoc H2O::triplePressure
     */
    OPM_HOST_DEVICE static Scalar triplePressure()
    { return H2O::triplePressure(); /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure water
     *        at a given temperature.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param T Absolute temperature of the system in \f$\mathrm{[K]}\f$
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation vaporPressure(const Evaluation& T)
    { return H2O::vaporPressure(T); /* [N/m^2] */ }

    /*!
     * \copydoc Component::gasEnthalpy
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation gasEnthalpy(const Evaluation& temperature,
                                  const Evaluation& pressure)
    { return H2O::gasEnthalpy(temperature, pressure); /* [J/kg] */ }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component in liquid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param salinity Salinity of component
     *
     *
     * Equations given in:
     * - Palliser & McKibbin 1997
     * - Michaelides 1981
     * - Daubert & Danner 1989
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation liquidEnthalpy(const Evaluation& temperature,
                                     const Evaluation& pressure,
                                     const Evaluation& salinity)
    {
        // Numerical coefficents from Palliser and McKibbin
        static constexpr Scalar f[] = {
            2.63500e-1, 7.48368e-6, 1.44611e-6, -3.80860e-10
        };

        // Numerical coefficents from Michaelides for the enthalpy of brine
        static constexpr Scalar a[4][3] = {
            { -9633.6, -4080.0, +286.49 },
            { +166.58, +68.577, -4.6856 },
            { -0.90963, -0.36524, +0.249667e-1 },
            { +0.17965e-2, +0.71924e-3, -0.4900e-4 }
        };

        const Evaluation theta = temperature - 273.15;

        Evaluation S = salinity;
        const Evaluation S_lSAT =
            f[0]
            + f[1]*theta
            + f[2]*pow(theta, 2)
            + f[3]*pow(theta, 3);

        // Regularization
        if (S > S_lSAT)
            S = S_lSAT;

        const Evaluation hw = H2O::liquidEnthalpy(temperature, pressure)/1e3; // [kJ/kg]

        // From Daubert and Danner
        const Evaluation h_NaCl =
            (3.6710e4*temperature
             + (6.2770e1/2)*temperature*temperature
             - (6.6670e-2/3)*temperature*temperature*temperature
             + (2.8000e-5/4)*pow(temperature, 4.0))/58.44e3
            - 2.045698e+02; // [kJ/kg]

        const Evaluation m = S/(1-S)/58.44e-3;

        Evaluation d_h = 0;
        for (int i = 0; i<=3; ++i) {
            for (int j = 0; j <= 2; ++j) {
                d_h += a[i][j] * pow(theta, i) * pow(m, j);
            }
        }

        const Evaluation delta_h = 4.184/(1e3 + (58.44 * m))*d_h;

        // Enthalpy of brine
        const Evaluation h_ls = (1-S)*hw + S*h_NaCl + S*delta_h; // [kJ/kg]
        return h_ls*1e3; // convert to [J/kg]
    }


    /*!
     * \copydoc H2O::liquidHeatCapacity
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation liquidHeatCapacity(const Evaluation& temperature,
                                         const Evaluation& pressure)
    {
        Scalar eps = scalarValue(temperature)*1e-8;
        return (liquidEnthalpy(temperature + eps, pressure) - liquidEnthalpy(temperature, pressure))/eps;
    }

    /*!
     * \copydoc H2O::gasHeatCapacity
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation gasHeatCapacity(const Evaluation& temperature,
                                      const Evaluation& pressure)
    { return H2O::gasHeatCapacity(temperature, pressure); }

    /*!
     * \copydoc H2O::gasInternalEnergy
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation gasInternalEnergy(const Evaluation& temperature,
                                        const Evaluation& pressure)
    {
        return
            gasEnthalpy(temperature, pressure) -
            pressure/gasDensity(temperature, pressure);
    }

    /*!
     * \copydoc H2O::liquidInternalEnergy
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation liquidInternalEnergy(const Evaluation& temperature,
                                           const Evaluation& pressure)
    {
        return
            liquidEnthalpy(temperature, pressure) -
            pressure/liquidDensity(temperature, pressure);
    }

    /*!
     * \copydoc H2O::gasDensity
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation gasDensity(const Evaluation& temperature, const Evaluation& pressure)
    { return H2O::gasDensity(temperature, pressure); }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the liquid component at a given pressure in \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param salinity Salinity of component
     * \param extrapolate True to use extrapolation
     *
     * Equations given in:
     * - Batzle & Wang (1992)
     * - cited by: Adams & Bachu in Geofluids (2002) 2, 257-271
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation liquidDensity(const Evaluation& temperature, const Evaluation& pressure,
                                                    const Evaluation& salinity, bool extrapolate = false)
    {
        const Evaluation rhow = H2O::liquidDensity(temperature, pressure, extrapolate);
        return liquidDensity(temperature, pressure, salinity, rhow);
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the liquid component at
     *        a given pressure in \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param salinity Concentration of dissolved salts
     * \param rhow Density of liquid
     *
     * Equations given in:
     * - Batzle & Wang (1992)
     * - cited by: Adams & Bachu in Geofluids (2002) 2, 257-271
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation liquidDensity(const Evaluation& temperature,
                                                    const Evaluation& pressure,
                                                    const Evaluation& salinity,
                                                    const Evaluation& rhow)
    {
        Evaluation tempC = temperature - 273.15;
        Evaluation pMPa = pressure/1.0E6;
        return
            rhow +
            1000*salinity*(
                0.668 +
                0.44*salinity +
                1.0E-6*(
                    300*pMPa -
                    2400*pMPa*salinity +
                    tempC*(
                        80.0 +
                        3*tempC -
                        3300*salinity -
                        13*pMPa +
                        47*pMPa*salinity)));
    }

    /*!
     * Density of liquid brine using multicomponent salts from Laliberte & Cooper,
     * Journal of Chemical and Engineering Data, Vol. 49, No. 5, 2004.
     *
     * @param temperature Temperature [K]
     * @param pressure Pressure [Pa], only used if useH2ODensity is true
     * @param salinity Array of salt-ion mass fractions in the total brine [kg ion/kg brine]
     * @param useH2ODensity True to compute the pure water density from H2O::liquidDensity
     *        (pressure dependent) instead of the Laliberte & Cooper correlation
     *        (pressure independent)
     * @param extrapolate True to use extrapolation, only used if useH2ODensity is true
     * @param pairing How ions are paired into electrolytes
     * @return Brine density [kg/m3]
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation
    liquidDensityMulticompSalt(const Evaluation& temperature,
                               const Evaluation& pressure,
                               const SaltArray<Evaluation, SaltMassFraction>& salinity,
                               bool useH2ODensity = false,
                               bool extrapolate = false,
                               PairingMethod pairing = PairingMethod::Proportional)
    {
        // Density of pure water
        const Evaluation rhow = useH2ODensity
            ? H2O::liquidDensity(temperature, pressure, extrapolate)
            : liquidDensityPureWaterLaliberteCooper_(temperature);

        // Return pure water density if salinity is zero
        if (!salinity.any_nonzero()) {
            return rhow;
        }
        return liquidDensityLaliberteCooper(temperature, salinity, rhow, pairing);
    }

    /*!
     * Density of liquid brine using multicomponent salts from Laliberte & Cooper,
     * Journal of Chemical and Engineering Data, Vol. 49, No. 5, 2004 (all equations refs. are
     * from this paper).
     *
     * @param temperature Temperature [K]
     * @param salinity Array of salt ions [kg ion/kg water]
     * @param rhow Density pure water [kg/m3]
     * @param pairing How ions are paired into electrolytes
     * @return Brine density [kg/m3]
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation
    liquidDensityLaliberteCooper(const Evaluation& temperature,
                                 const SaltArray<Evaluation, SaltMassFraction>& salinity,
                                 const Evaluation& rhow,
                                 PairingMethod pairing = PairingMethod::Proportional)
    {
#if OPM_IS_INSIDE_DEVICE_FUNCTION
        // Electrolyte speciation (SaltElectrolytes::fromSaltComponents) throws on unknown
        // ions and electrolytes, and is therefore not supported on the GPU yet.
        assert(false && "liquidDensityLaliberteCooper is not supported on the GPU");
        return rhow;
#else
        // Generate (valid) electrolytes from salt components
        auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity, pairing);

        // Sum of electrolyte mass fractions
        // OBS: If there are residual salt molalities when generating electrolytes, sumSalinity
        // will not be equal to salinity.sum(); see SaltElectrolytes::fromSaltComponents()
        Evaluation sumSalinity = std::accumulate(
            electrolytes.begin(), electrolytes.end(), Evaluation(0.0),
            [](const Evaluation& acc, const auto& salt) { return acc + salt.second; });

        // Calculate brine density with Eq. (6)
        Evaluation tempC = temperature - 273.15;
        Evaluation sumAppVolTimesMf = 0.0;
        for (const auto& salt : electrolytes) {
            if (salt.first == Electrolyte::BOH3) {
                sumAppVolTimesMf += salt.second * boricAcidSpecificVolume_(temperature);
                continue;
            }

            // Get constants
            const auto c = densityCoefficients(salt.first);

            // Calculate apparent volume with Eq. (10) (where we use (1 - w_H2O) instead of w_i)
            // multiplied by mass fraction salt electrolyte to get the sum term in Eq. (6)
            Evaluation tmpExponent = tempC + c[4];
            sumAppVolTimesMf +=
                salt.second * (sumSalinity + c[2] + c[3] * tempC)
                / ((c[0] * sumSalinity + c[1]) * exp(1e-6 * tmpExponent * tmpExponent));
        }

        // Calculate the rest of Eq. (6) to get brine density
        auto wH2O = 1.0 - sumSalinity;
        Evaluation rho = 1.0 / (wH2O / rhow + sumAppVolTimesMf);

        return rho;
#endif
    }

    /*!
     * \copydoc H2O::gasPressure
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation gasPressure(const Evaluation& temperature, const Evaluation& density)
    { return H2O::gasPressure(temperature, density); }

    /*!
     * \copydoc H2O::liquidPressure
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation liquidPressure(const Evaluation& temperature, const Evaluation& density)
    {
        // We use the newton method for this. For the initial value we
        // assume the pressure to be 10% higher than the vapor
        // pressure
        Evaluation pressure = 1.1*vaporPressure(temperature);
        Scalar eps = scalarValue(pressure)*1e-7;

        Evaluation deltaP = pressure*2;
        for (int i = 0;
             i < 5
                 && std::abs(scalarValue(pressure)*1e-9) < std::abs(scalarValue(deltaP));
             ++i)
        {
            const Evaluation f = liquidDensity(temperature, pressure) - density;

            Evaluation df_dp = liquidDensity(temperature, pressure + eps);
            df_dp -= liquidDensity(temperature, pressure - eps);
            df_dp /= 2*eps;

            deltaP = - f/df_dp;

            pressure += deltaP;
        }

        return pressure;
    }

    /*!
     * \copydoc H2O::gasViscosity
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation gasViscosity(const Evaluation& temperature, const Evaluation& pressure)
    { return H2O::gasViscosity(temperature, pressure); }

    /*!
     * \brief The dynamic liquid viscosity \f$\mathrm{[Pa*s]}\f$ of the pure component.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param salinity Salinity of component
     *
     * Equation given in:
     * - Batzle & Wang (1992)
     * - cited by: Bachu & Adams (2002)
     *   "Equations of State for basin geofluids"
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation liquidViscosity(const Evaluation& temperature,
                                                      const Evaluation& /*pressure*/,
                                                      const Evaluation& salinity)
    {
        Evaluation T_C = temperature - 273.15;
        if(temperature <= 275.) // regularization
            T_C = 275.0;

        Evaluation A = (0.42*Opm::pow((Opm::pow(salinity, 0.8)-0.17), 2) + 0.045)*pow(T_C, 0.8);
        Evaluation mu_brine = 0.1 + 0.333*salinity + (1.65+91.9*salinity*salinity*salinity)*exp(-A);

        return mu_brine/1000.0; // convert to [Pa s] (todo: check if correct cP->Pa s is times 10...)
    }

    /*!
     * Viscosity of liquid brine using multicomponent salts from Laliberte, Journal of Chemical
     * and Engineering Data, Vol. 52, No. 2, 2007 and corrections paper Laliberte, Journal of
     * Chemical and Engineering Data, Vol. 52, No. 4, 2007.
     *
     * @param temperature Temperature [K]
     * @param salinity Array of salt-ion mass fractions in the total brine [kg ion/kg brine]
     * @param pairing How ions are paired into electrolytes
     * @return Viscosity of brine [Pa*s]
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation
    liquidViscosityMulticompSalt(const Evaluation& temperature,
                                 const SaltArray<Evaluation, SaltMassFraction>& salinity,
                                 PairingMethod pairing = PairingMethod::Proportional)
    {
        // Pure water viscosity
        const Evaluation muW = liquidViscosityPureWaterLaliberte_(temperature);

        // Return pure water viscosity for zero salinity
        if (!salinity.any_nonzero()) {
            return muW * 1e-3;
        }

        return liquidViscosityLaliberte(temperature, salinity, muW, pairing);
    }

    /*!
     * Viscosity of liquid brine using multicomponent salts from Laliberte, Journal of Chemical
     * and Engineering Data, Vol. 52, No. 2, 2007 and corrections paper Laliberte, Journal of
     * Chemical and Engineering Data, Vol. 52, No. 4, 2007 (all equations refs. are from these
     * papers).
     *
     * @param temperature Temperature [K]
     * @param salinity Array of salt ions [kg ion/kg water]
     * @param muW Viscosity of pure water [mPa*s]
     * @param pairing How ions are paired into electrolytes
     * @return Viscosity of brine [Pa*s]
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation
    liquidViscosityLaliberte(const Evaluation& temperature,
                             const SaltArray<Evaluation, SaltMassFraction>& salinity,
                             const Evaluation& muW,
                             PairingMethod pairing = PairingMethod::Proportional)
    {
#if OPM_IS_INSIDE_DEVICE_FUNCTION
        // Electrolyte speciation (SaltElectrolytes::fromSaltComponents) throws on unknown
        // ions and electrolytes, and is therefore not supported on the GPU yet.
        assert(false && "liquidViscosityLaliberte is not supported on the GPU");
        return muW * 1e-3;
#else
        // Generate (valid) electrolytes from salt components
        auto electrolytes = SaltElectrolytes<Scalar>::fromSaltComponents(salinity, pairing);

        // Sum of electrolyte mass fractions
        // OBS: If there are residual salt molalities when generating electrolytes, sumSalinity
        // will not be equal to salinity.sum(); see SaltElectrolytes::fromSaltComponents()
        Evaluation sumSalinity = std::accumulate(
            electrolytes.begin(), electrolytes.end(), Evaluation(0.0),
            [](const Evaluation& acc, const auto& salt) { return acc + salt.second; });

        // ln-viscosity of mix from Eq. (8)
        Evaluation tempC = temperature - 273.15;
        Evaluation lnMuMix = 0.0;
        for (const auto& salt : electrolytes) {
            // Laliberte gives no viscosity set for boric acid, so it counts as water
            if (salt.first == Electrolyte::BOH3) {
                lnMuMix += salt.second * log(muW);
                continue;
            }

            // Get constants
            const auto v = viscosityCoefficients(salt.first);

            // Calculate ln of Eq. (12).
            // OBS: Different from Eq. (13) in correction paper (which seems to be incorrect)
            lnMuMix += salt.second
                * ((v[0] * pow(sumSalinity, v[1]) + v[2]) / (v[3] * tempC + 1.0)
                - log(v[4] * pow(sumSalinity, v[5]) + 1.0));
        }

        // Calculate the rest of Eq. (8)
        lnMuMix += (1.0 - sumSalinity) * log(muW);

        // OBS: mPa*s to Pa*s
        return exp(lnMuMix) * 1e-3;
#endif
    }

private:
    using Electrolyte = typename SaltElectrolytes<Scalar>::Electrolyte;
    static constexpr std::size_t numElectrolytes = SaltElectrolytes<Scalar>::numElectrolytes;

    //Molar mass salt (assumes pure NaCl) [kg/mol]
    static constexpr Scalar mM_salt()
    {
        return 58.44e-3;
    }

    /*!
     * Get coefficients for calculating density with electrolyte
     *
     * @param electrolyte Electrolyte index
     * @return Coefficients for electrolyte
     */
    static std::span<const Scalar, 5>
    densityCoefficients(const Electrolyte electrolyte)
    {
        return cDens[static_cast<std::size_t>(electrolyte)];
    }

    /*!
     * Get coefficients for calculating viscosity with electrolyte
     *
     * @param electrolyte Electrolyte index
     * @return Coefficients for electrolyte
     */
    static std::span<const Scalar, 6>
    viscosityCoefficients(const Electrolyte electrolyte)
    {
        return cVisc[static_cast<std::size_t>(electrolyte)];
    }

    /*!
     * Density of pure water from Laliberte & Cooper, Journal of Chemical and Engineering Data,
     * Vol. 49, No. 5, 2004.
     *
     * @param temperature Temperature [K]
     * @return Density [kg/m3]
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation
    liquidDensityPureWaterLaliberteCooper_(const Evaluation& temperature)
    {
        const auto tempC = temperature - 273.15;
        return (((((-2.8054253e-10 * tempC + 1.0556302e-7) * tempC - 4.6170461e-5) * tempC
                - 0.0079870401) * tempC + 16.945176) * tempC + 999.83952)
            / (1.0 + 0.01687985 * tempC);
    }

    /*!
     * Viscosity of pure water from Laliberte, Journal of Chemical and Engineering Data, Vol. 52,
     * No. 2, 2007.
     *
     * @param temperature Temperature [K]
     * @return Viscosity [mPa*s]
     */
    template <class Evaluation>
    OPM_HOST_DEVICE static Evaluation
    liquidViscosityPureWaterLaliberte_(const Evaluation& temperature)
    {
        const auto tempC = temperature - 273.15;
        return (tempC + 246.0) / ((0.05594 * tempC + 5.2842) * tempC + 137.37);
    }

    /*!
     * Apparent specific volume of dissolved boric acid. The apparent molar volume is a line
     * in T fitted to the infinite-dilution values of Abdulagatov and Azizov, J. Solution Chem.
     * 33, 2004, Table VII, from 296 to 523 K. They found pressure and concentration effects
     * small below 500 K, so both are ignored.
     *
     * @param temperature Temperature [K]
     * @return Apparent specific volume [m3/kg]
     */
    template <class Evaluation>
    static Evaluation boricAcidSpecificVolume_(const Evaluation& temperature)
    {
        const Evaluation molarVolume = (30.882 + 0.031814 * temperature) * 1e-6; // m3/mol
        return molarVolume / saltMolarMass<Scalar>(SaltIndex::B);
    }

    static constexpr Scalar unused_ = std::numeric_limits<Scalar>::quiet_NaN();

    // Density coefficients c0 to c4 in Electrolyte order. The first seven are from Laliberte
    // and Cooper, J. Chem. Eng. Data 49, 2004, Table 3; the rest from the refit in Laliberte,
    // J. Chem. Eng. Data 54, 2009, Supporting Information.
    // clang-format off
    static constexpr std::array<std::array<Scalar, 5>, numElectrolytes> cDens = {{
        // CaCl2
        {-0.63254, 0.93995, 4.2785, 0.048319, 3180.9},
        // KCl
        {-0.46782, 4.30800, 2.3780, 0.022044, 2714.0},
        // K2SO4
        {-2.6681e-5, 3.0412e-5, 0.97118, 0.019816, 4366.1},
        // MgCl2
        {-0.00051500, 0.0013444, 0.58358, 0.0085832, 3843.6},
        // MgSO4
        {3.9412e-7, 1.4425e-6, -0.05372, 0.002062, 4563.3},
        // NaCl
        {-0.00433, 0.06471, 1.01660, 0.014624, 3315.6},
        // Na2SO4
        {-1.2095e-7, 4.3474e-7, 0.15364, 0.0072514, 4731.5},
        // LiCl
        {1777.71168869463, 208.095675885873, 0.0924032897423372, -9.65138194640463e-05,
         -303.212122198705},
        // Li2SO4
        {0.00711493013455434, 0.00869246384846451, 0.151649349676449, 0.00360862647418256,
         3458.25184736606},
        // SrCl2
        {1.35341758029363e-06, -7.48767819877433e-07, -1.93563339996439, 0.010704329282332,
         -4882.14470009011},
        // NaBr
        {23.765937876294, 210.721942848473, 2.01023472140861, 0.0163941235513688, 1941.82908042744},
        // KBr
        {-13.151388032057, 190.74540854419, 5.39738604657711, 0.0402662461365671, 2153.22549041443},
        // NaHCO3
        {0.0377378060615898, 0.00710029652458366, 0.0735733388535116, 0.00168369765511221,
         3291.38382252025},
        // Na2CO3
        {5.45486543455015e-07, 1.72154929405038e-07, -0.0825850207323492, 0.00168822949783412,
         4705.39822185746},
        // KHCO3
        {0.0458479253730966, 0.0749127808255568, 0.505766905030703, 0.00643951863958203,
         3163.65175558093},
        // K2CO3
        {2290.10867876218, 1038.77363895405, 0.109769327294915, 0.000497783340935823,
         -74.0786862142127},
        // B(OH)3, see boricAcidSpecificVolume_
        {unused_, unused_, unused_, unused_, unused_},
    }};
    // clang-format on

    // Viscosity coefficients v1 to v6 in Electrolyte order. The first seven are from
    // Laliberte, J. Chem. Eng. Data 52 (2), 2007 and its correction, J. Chem. Eng. Data 52 (4),
    // 2007; the rest from Laliberte, J. Chem. Eng. Data 54, 2009, Supporting Information.
    // clang-format off
    static constexpr std::array<std::array<Scalar, 6>, numElectrolytes> cVisc = {{
        // CaCl2
        {32.028, 0.78792, -1.1495, 0.0026995, 780860, 5.8442},
        // KCl
        {6.4883, 1.3175, -0.77785, 0.092722, -1.3, 2.0811},
        // K2SO4
        {-983.76, 0.00019331, 984.52, 0.0038473, -9.5001, 2.1916},
        // MgCl2
        {24.032, 2.2694, 3.7108, 0.021853, -1.1236, 0.14474},
        // MgSO4
        {72.269, 2.2238, 6.6037, 0.0079004, 3340.1, 6.1304},
        // NaCl
        {16.222, 1.3229, 1.4849, 0.0074691, 30.78, 2.0583},
        // Na2SO4
        {26.519, 1.5746, 3.4966, 0.010388, 106.23, 2.9738},
        // LiCl
        {18.6178234588751, 0.773036318483134, 2.15660166137342, 0.00435445218019151,
         1023.45333257758, 2.38089336779273},
        // Li2SO4
        {42.4909603154037, 2.21333353374412, 6.06402999044027, 0.00814290129030605,
         24.1674397087677, 2.96004938607259},
        // SrCl2
        {8.03309316438313, 1.85314812521892, 1.99358226403085, 0.00848650077967969,
         -4.98502701890093, 2.88704932687605},
        // NaBr
        {13.02905395428, 1.74783698289335, 0.604127008790876, 0.0108036956912202, 17.6807414205322,
         2.3830591098608},
        // KBr
        {348.31972605609, -0.000302529305665952, -349.153198551268, -0.00425795862818544,
         -1.10436123629246, 0.763153789481199},
        // NaHCO3
        {31.5819207092739, 0.0794532561718167, -16.9096255040761, 0.0013751751230389,
         44188.6551718116, 1.88902235502575},
        // Na2CO3
        {27.6831996816743, 1.56324040665499, 6.80529471973068, 0.0150777050822261, 22452.3440720574,
         28.2362397327323},
        // KHCO3
        {0.23545396225936, -1.07084725909778, 2.35937892744089, 0.0639877898997282,
         -1.00034737120549, 0.501215910081986},
        // K2CO3
        {41.3874564880263, 3.67791349310373, 3.26807400039543, 0.00876702412403997,
         8684.04398141182, 12.7276295460364},
        // B(OH)3, not used
        {unused_, unused_, unused_, unused_, unused_, unused_},
    }};
    // clang-format on
};

} // namespace Opm

#endif
