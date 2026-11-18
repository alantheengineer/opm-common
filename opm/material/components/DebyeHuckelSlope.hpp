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
 * \copydoc Opm::DebyeHuckelSlope
 */
#ifndef OPM_DEBYE_HUCKEL_SLOPE_HPP
#define OPM_DEBYE_HUCKEL_SLOPE_HPP

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/components/H2O.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace Opm
{

/*!
 * \brief Debye-Huckel slopes of water from the dielectric constant of Bradley and Pitzer,
 * J. Phys. Chem. 83, 1979, and IAPWS-97 water density.
 *
 * The volume slope needs the pressure derivative of the osmotic slope, which involves the
 * compressibility of water. Inside a Newton loop that would mean a second derivative, so the
 * volume slope is tabulated once in T and p and read through a bicubic Hermite patch, which
 * gives smooth values and first derivatives.
 */
template <class Scalar>
class DebyeHuckelSlope
{
    using Water = H2O<double>;

public:
    /*!
     * Dielectric constant of water, Bradley and Pitzer (1979), eq 1 and Table I
     *
     * @param temperature Temperature [K]
     * @param pressure Pressure [Pa]
     */
    static double dielectricConstant(const double temperature, const double pressure)
    {
        const auto [d1000, c, b] = dielectricTerms_(temperature);
        return d1000 + c * std::log((b + pressure * 1e-5) / (b + 1000.0));
    }

    /*!
     * Debye-Huckel slope for the osmotic coefficient, A_phi
     *
     * @param temperature Temperature [K]
     * @param pressure Pressure [Pa]
     * @return A_phi [kg^0.5 mol^-0.5]
     */
    static double osmoticSlope(const double temperature, const double pressure)
    {
        const double rhoW = waterDensity_(temperature, pressure);
        const double d = dielectricConstant(temperature, pressure);
        const double lengthScale = electronCharge_ * electronCharge_
            / (4.0 * pi_ * vacuumPermittivity_ * d * boltzmann_ * temperature);
        return std::sqrt(2.0 * pi_ * avogadro_ * rhoW) * std::pow(lengthScale, 1.5) / 3.0;
    }

    /*!
     * Volumetric Debye-Huckel slope A_V = -4RT (dA_phi/dP)_T, computed directly
     *
     * @param temperature Temperature [K]
     * @param pressure Pressure [Pa]
     * @return A_V [cm3 kg^0.5 mol^-1.5]
     */
    static double volumeSlopeDirect(const double temperature, const double pressure)
    {
        // Pressure derivatives per MPa
        const double dp = 1e4; // Pa
        const double dLnRhoW = (std::log(waterDensity_(temperature, pressure + dp))
                                - std::log(waterDensity_(temperature, pressure - dp)))
            / (2.0 * dp * 1e-6);
        const auto [d1000, c, b] = dielectricTerms_(temperature);
        const double pBar = pressure * 1e-5;
        const double d = d1000 + c * std::log((b + pBar) / (b + 1000.0));
        const double dLnD = 10.0 * c / ((b + pBar) * d);

        const double aPhi = osmoticSlope(temperature, pressure);
        return -4.0 * gasConstant_ * temperature * aPhi * (0.5 * dLnRhoW - 1.5 * dLnD);
    }

    /*!
     * Volumetric Debye-Huckel slope from the table, with derivatives
     *
     * @param temperature Temperature [K]
     * @param pressure Pressure [Pa]
     * @return A_V [cm3 kg^0.5 mol^-1.5]
     */
    template <class Evaluation>
    static Evaluation volumeSlope(const Evaluation& temperature, const Evaluation& pressure)
    {
        const auto& table = table_();
        const Evaluation x = (temperature - tableFirstT_) / tableStepT_;
        const Evaluation y = (pressure - tableFirstP_) / tableStepP_;
        const int i = std::clamp(static_cast<int>(std::floor(scalarValue(x))), 0, numT_ - 2);
        const int j = std::clamp(static_cast<int>(std::floor(scalarValue(y))), 0, numP_ - 2);
        const Evaluation u = x - i;
        const Evaluation v = y - j;

        // Hermite basis: value and slope at each end of the unit interval
        const auto basis = [](const Evaluation& s) {
            const Evaluation s2 = s * s;
            const Evaluation s3 = s2 * s;
            return std::array<Evaluation, 4> {
                2 * s3 - 3 * s2 + 1, s3 - 2 * s2 + s, -2 * s3 + 3 * s2, s3 - s2};
        };
        const auto bu = basis(u);
        const auto bv = basis(v);

        Evaluation result = 0.0;
        for (int a = 0; a < 2; ++a) {
            for (int b = 0; b < 2; ++b) {
                const Node& n = table[i + a][j + b];
                const Evaluation& hu = bu[2 * a];
                const Evaluation& gu = bu[2 * a + 1];
                const Evaluation& hv = bv[2 * b];
                const Evaluation& gv = bv[2 * b + 1];
                result += hu * hv * n.value + gu * hv * n.dx + hu * gv * n.dy + gu * gv * n.dxy;
            }
        }
        return result;
    }

private:
    static constexpr double pi_ = 3.14159265358979323846;
    static constexpr double avogadro_ = 6.0221415e23; // 1/mol
    static constexpr double electronCharge_ = 1.60217733e-19; // C
    static constexpr double boltzmann_ = 1.3806505e-23; // J/K
    static constexpr double vacuumPermittivity_ = 8.854187817e-12; // F/m
    static constexpr double gasConstant_ = 8.314472; // cm3 MPa/(mol K)

    // Bradley and Pitzer (1979), Table I
    static constexpr std::array<double, 9> u_ = {
        3.4279e2, -5.0866e-3, 9.4690e-7, -2.0525, 3.1159e3, -1.8289e2, -8.0325e3, 4.2142e6, 2.1417};

    struct DielectricTerms
    {
        double d1000;
        double c;
        double b;
    };

    static DielectricTerms dielectricTerms_(const double temperature)
    {
        const double t = temperature;
        return {u_[0] * std::exp(u_[1] * t + u_[2] * t * t),
                u_[3] + u_[4] / (u_[5] + t),
                u_[6] + u_[7] / t + u_[8] * t};
    }

    // IAPWS-97 liquid density [kg/m3]; below the vapour pressure OPM extends it smoothly
    static double waterDensity_(const double temperature, const double pressure)
    {
        return Water::liquidDensity(temperature, pressure, /*extrapolate=*/true);
    }

    // Table from 0 to 350 C every 2.5 C, and from 0.1 to 150.1 MPa every 2.5 MPa
    static constexpr double tableFirstT_ = 273.15;
    static constexpr double tableStepT_ = 2.5;
    static constexpr double tableFirstP_ = 1.0e5;
    static constexpr double tableStepP_ = 2.5e6;
    static constexpr int numT_ = 141;
    static constexpr int numP_ = 61;

    // Value and scaled derivatives, per unit step in T and in p
    struct Node
    {
        Scalar value;
        Scalar dx;
        Scalar dy;
        Scalar dxy;
    };
    using Table = std::array<std::array<Node, numP_>, numT_>;

    static const Table& table_()
    {
        static const Table table = buildTable_();
        return table;
    }

    static Table buildTable_()
    {
        std::array<std::array<double, numP_>, numT_> f {};
        for (int i = 0; i < numT_; ++i) {
            for (int j = 0; j < numP_; ++j) {
                const double p = tableFirstP_ + j * tableStepP_;
                f[i][j] = volumeSlopeDirect(tableFirstT_ + i * tableStepT_, p);
            }
        }

        // Derivatives by central differences, one-sided at the edges
        const auto diff = [](int k, int n) {
            const int lo = std::max(k - 1, 0);
            const int hi = std::min(k + 1, n - 1);
            return std::array<int, 3> {lo, hi, hi - lo};
        };
        Table table {};
        for (int i = 0; i < numT_; ++i) {
            const auto [il, ih, di] = diff(i, numT_);
            for (int j = 0; j < numP_; ++j) {
                const auto [jl, jh, dj] = diff(j, numP_);
                Node& n = table[i][j];
                n.value = f[i][j];
                n.dx = (f[ih][j] - f[il][j]) / di;
                n.dy = (f[i][jh] - f[i][jl]) / dj;
                n.dxy = (f[ih][jh] - f[ih][jl] - f[il][jh] + f[il][jl]) / (di * dj);
            }
        }
        return table;
    }
};

} // namespace Opm

#endif // OPM_DEBYE_HUCKEL_SLOPE_HPP
