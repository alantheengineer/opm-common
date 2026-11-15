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
 * \copydoc Opm::BOH3
 */
#ifndef OPM_BOH3_HPP
#define OPM_BOH3_HPP

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/material/components/Component.hpp>

namespace Opm
{

/*!
 * \ingroup Components
 *
 * \brief Properties of B(OH)3 (boric acid), the form dissolved boron takes at reservoir pH
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class BOH3 : public Component<Scalar, BOH3<Scalar> >
{
public:
    /*!
     * \copydoc Component::name
     */
    static std::string_view name()
    {
        return "B(OH)3";
    }

    /*!
     * \copydoc Component::molarMass
     */
    OPM_HOST_DEVICE static Scalar molarMass()
    {
        return 61.83e-3;
    }
};

}

#endif //OPM_BOH3_HPP
