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
 * \copydoc Opm::CO3Ion
 */
#ifndef OPM_CO3ION_HPP
#define OPM_CO3ION_HPP

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/material/components/Component.hpp>

namespace Opm
{

/*!
 * \ingroup Components
 *
 * \brief Properties of CO3-- (carbonate) ion
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class CO3Ion : public Component<Scalar, CO3Ion<Scalar> >
{
public:
    /*!
     * \copydoc Component::name
     */
    static std::string_view name()
    {
        return "CO3--";
    }

    /*!
     * \copydoc Component::molarMass
     */
    OPM_HOST_DEVICE static Scalar molarMass()
    {
        return 60.008e-3;
    }
};

}

#endif //OPM_CO3ION_HPP
