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
#ifndef OPM_SALTARRAY_HPP
#define OPM_SALTARRAY_HPP

#include <opm/input/eclipse/Deck/DeckRecord.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/FixedCapacityList.hpp>
#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/material/components/BOH3.hpp>
#include <opm/material/components/BrIon.hpp>
#include <opm/material/components/CO3Ion.hpp>
#include <opm/material/components/CaIon.hpp>
#include <opm/material/components/ClIon.hpp>
#include <opm/material/components/HCO3Ion.hpp>
#include <opm/material/components/KIon.hpp>
#include <opm/material/components/LiIon.hpp>
#include <opm/material/components/MgIon.hpp>
#include <opm/material/components/NaIon.hpp>
#include <opm/material/components/SO4Ion.hpp>
#include <opm/material/components/SrIon.hpp>

#include <array>
#include <cassert>
#include <cstddef>
#include <utility>

namespace Opm
{

// SaltUnit tags for SaltArray
struct SaltMolality
{
};

struct SaltMassFraction
{
};

struct SaltMoleFraction
{
};

// Helper conversion struct
template <class FromSaltUnit, class ToSaltUnit>
struct SaltUnitConverter;

/// Salt ion indices for SaltArray. New ions go at the end so that serialized arrays keep
/// their meaning. B is total dissolved boron, carried as B(OH)3.
// clang-format off
enum class SaltIndex : std::size_t
{
    NA, K, CA, MG, CL, SO4,
    LI, SR, BR, HCO3, CO3, B,
    NCOMP // DO NOT CHANGE! Must be the last enum to set size of SaltArray!
};
// clang-format on

/// Number of salt ions in SaltIndex
inline constexpr std::size_t numSaltIons = static_cast<std::size_t>(SaltIndex::NCOMP);

/// Largest ion charge in SaltIndex, by absolute value
inline constexpr short maxIonCharge = 2;

/// Salt ion indices held without heap allocation
using SaltIndexList = FixedCapacityList<SaltIndex, numSaltIons>;

/// Enum for classifying cations, anions and neutral solutes
enum class IonType { Cation, Anion, Neutral };

/*!
 * Converts ion string name to index
 *
 * @param s Ion name
 * @return Index of salt ion
 */
inline SaltIndex
saltIndexFromString(const std::string& s)
{
    if (s == "NA") {
        return SaltIndex::NA;
    }
    if (s == "K") {
        return SaltIndex::K;
    }
    if (s == "CA") {
        return SaltIndex::CA;
    }
    if (s == "MG") {
        return SaltIndex::MG;
    }
    if (s == "CL") {
        return SaltIndex::CL;
    }
    if (s == "SO4") {
        return SaltIndex::SO4;
    }
    if (s == "LI") {
        return SaltIndex::LI;
    }
    if (s == "SR") {
        return SaltIndex::SR;
    }
    if (s == "BR") {
        return SaltIndex::BR;
    }
    if (s == "HCO3") {
        return SaltIndex::HCO3;
    }
    if (s == "CO3") {
        return SaltIndex::CO3;
    }
    if (s == "B") {
        return SaltIndex::B;
    }
    throw std::invalid_argument("SaltIndex not valid!");
}

/*!
 * Categorize ion in cation or anion
 *
 * @param s Index of salt ion
 * @return Cation or anion category
 */
inline IonType
categorizeIon(const SaltIndex s)
{
    switch (s) {
    case SaltIndex::CA:
    case SaltIndex::NA:
    case SaltIndex::MG:
    case SaltIndex::K:
    case SaltIndex::LI:
    case SaltIndex::SR:
        return IonType::Cation;
    case SaltIndex::CL:
    case SaltIndex::SO4:
    case SaltIndex::BR:
    case SaltIndex::HCO3:
    case SaltIndex::CO3:
        return IonType::Anion;
    case SaltIndex::B:
        return IonType::Neutral;
    default:
        throw std::runtime_error("Unknown SaltIndex");
    }
}

/*!
 * Returns the ion strength or salt ion. Positive for cations; negative for anions.
 *
 * @param s Index of salt ion
 * @return Ion strength
 */
inline constexpr short
ionStrength(const SaltIndex s)
{
    switch (s) {
    case SaltIndex::CA:
        return 2;
    case SaltIndex::NA:
        return 1;
    case SaltIndex::MG:
        return 2;
    case SaltIndex::K:
        return 1;
    case SaltIndex::CL:
        return -1;
    case SaltIndex::SO4:
        return -2;
    case SaltIndex::LI:
        return 1;
    case SaltIndex::SR:
        return 2;
    case SaltIndex::BR:
        return -1;
    case SaltIndex::HCO3:
        return -1;
    case SaltIndex::CO3:
        return -2;
    case SaltIndex::B:
        return 0;
    default:
        throw std::runtime_error("Unknown SaltIndex");
    }
}

namespace detail
{
    constexpr std::size_t countSaltIons(const short sign)
    {
        std::size_t n = 0;
        for (std::size_t i = 0; i < numSaltIons; ++i) {
            if (ionStrength(static_cast<SaltIndex>(i)) * sign > 0) {
                ++n;
            }
        }
        return n;
    }
} // namespace detail

/// Number of cation and anion species in SaltIndex
inline constexpr std::size_t numSaltCations = detail::countSaltIons(1);
inline constexpr std::size_t numSaltAnions = detail::countSaltIons(-1);

/*!
 * Get molar mass of salt ion
 *
 * @param ind Index of salt ion
 * @return Molar mass
 */
template <class Scalar>
OPM_HOST_DEVICE Scalar
saltMolarMass(const SaltIndex ind)
{
    switch (ind) {
    case SaltIndex::NA:
        return NaIon<Scalar>::molarMass();
    case SaltIndex::K:
        return KIon<Scalar>::molarMass();
    case SaltIndex::CA:
        return CaIon<Scalar>::molarMass();
    case SaltIndex::CL:
        return ClIon<Scalar>::molarMass();
    case SaltIndex::MG:
        return MgIon<Scalar>::molarMass();
    case SaltIndex::SO4:
        return SO4Ion<Scalar>::molarMass();
    case SaltIndex::LI:
        return LiIon<Scalar>::molarMass();
    case SaltIndex::SR:
        return SrIon<Scalar>::molarMass();
    case SaltIndex::BR:
        return BrIon<Scalar>::molarMass();
    case SaltIndex::HCO3:
        return HCO3Ion<Scalar>::molarMass();
    case SaltIndex::CO3:
        return CO3Ion<Scalar>::molarMass();
    case SaltIndex::B:
        return BOH3<Scalar>::molarMass();
    default:
#if OPM_IS_INSIDE_DEVICE_FUNCTION
        assert(false && "SaltIndex not valid!");
        return Scalar{};
#else
        OPM_THROW(std::invalid_argument, "SaltIndex not valid!");
#endif
    }
}

/*!
 * Wrapper around std::array for salt components with some convenience functions and conversion
 * to SaltArrays with different units
 *
 * @tparam T Value type (typically Scalar or Evaluation)
 * @tparam SaltUnit Unit of salt array elements
 */
template <class T, class SaltUnit>
class SaltArray
{
public:
    using SaltIndexPair = std::pair<SaltIndexList, SaltIndexList>;

    static constexpr std::size_t ncomp = numSaltIons;

    /*!
     * Default constructor
     */
    SaltArray() = default;

    /*!
     * Default (not converting) copy constructor.
     */
    OPM_HOST_DEVICE SaltArray(const SaltArray&) = default;

    /*!
     * Conversion constructor
     *
     * @tparam U Value type different from T
     * @param other Salt array to convert from
     */
    template <class U>
    OPM_HOST_DEVICE explicit SaltArray(const SaltArray<U, SaltUnit>& other)
    {
        for (std::size_t i = 0; i < size(); ++i) {
            saltData_[i] = static_cast<T>(other[static_cast<SaltIndex>(i)]);
        }
    }

    /*!
     * Assign array elements from a deck array
     *
     * @param record Deck record array
     */
    void assignFromDeckRecord(const DeckRecord& record)
    {
        assert(record.size() <= saltData_.size());
        for (const auto& item : record) {
            auto recItem = static_cast<T>(item.get<double>(0));
            auto index = saltIndexFromString(item.name());
            (*this)[index] = recItem;
        }
    }

    /*!
     * Return read-only view of array element at a salt ion index
     *
     * @param ind Index of salt ion
     * @return Array element
     */
    OPM_HOST_DEVICE const T& operator[](const SaltIndex ind) const
    {
        return saltData_[static_cast<std::size_t>(ind)];
    }

    /*!
     * Return array element at a salt ion index
     *
     * @param ind Index of salt ion
     * @return Array element
     */
    OPM_HOST_DEVICE T& operator[](const SaltIndex ind)
    {
        return saltData_[static_cast<std::size_t>(ind)];
    }

    /*!
     * Check for equal SaltArray
     *
     * @param other Comparison array
     * @return True if arrays are equal
     */
    OPM_HOST_DEVICE bool operator==(const SaltArray& other) const
    {
        return saltData_ == other.saltData_;
    }

    /*!
     * Default copy assignment
     *
     * @param other Array to copy
     * @return Copy of other array
     */
    OPM_HOST_DEVICE SaltArray& operator=(const SaltArray& other) = default;

    /*!
     * Conversion copy assignment for simple assignment to new value type (that supports conversion)
     *
     * @tparam U Value type of other array
     * @param other Array to copy
     * @return Copy of array
     */
    template <class U>
    OPM_HOST_DEVICE SaltArray& operator=(const SaltArray<U, SaltUnit>& other)
    {
        for (std::size_t i = 0; i < size(); ++i) {
            saltData_[i] = static_cast<T>(other[static_cast<SaltIndex>(i)]);
        }
        return *this;
    }

    /*!
     * Array begin iterator
     *
     * @return Iterator
     */
    OPM_HOST_DEVICE auto begin()
    {
        return saltData_.begin();
    }

    /*!
     * Array end iterator
     *
     * @return Iterator
     */
    OPM_HOST_DEVICE auto end()
    {
        return saltData_.end();
    }

    /*!
     * Read-only view of array begin iterator
     *
     * @return Iterator
     */
    [[nodiscard]] OPM_HOST_DEVICE auto begin() const
    {
        return saltData_.begin();
    }

    /*!
     * Read-only view of array end iterator
     *
     * @return Iterator
     */
    [[nodiscard]] OPM_HOST_DEVICE auto end() const
    {
        return saltData_.end();
    }

    /*!
     * Size of array, which should be equal to number of supported ions in @ref SaltIndex
     *
     * @return Size of array
     */
    [[nodiscard]] OPM_HOST_DEVICE constexpr std::size_t size() const
    {
        return saltData_.size();
    }

    /*!
     * Sum of array elements.
     * @note For mass and mole fraction arrays, 1 - sum() is equal to mass/mole fraction of H2O
     *
     * @return Sum
     */
    OPM_HOST_DEVICE T sum() const
    {
        T s{};
        for (std::size_t i = 0; i < size(); ++i) {
            s += saltData_[i];
        }
        return s;
    }

    /*!
     * Fill all elements with zero
     */
    void clear() noexcept
    {
        saltData_.fill(T{});
    }

    /*!
     * Check if any elements are non-zero
     *
     * @return True if any element non-zero
     */
    [[nodiscard]] OPM_HOST_DEVICE bool any_nonzero() const noexcept
    {
        for (std::size_t i = 0; i < size(); ++i) {
            if (saltData_[i] != T{}) {
                return true;
            }
        }
        return false;
    }

    /*!
     * Get cations and anions of non-zero ions, cations by falling charge and anions by
     * rising charge. Ions of equal charge keep their SaltIndex order. Neutral solutes are
     * in neither list.
     *
     * @return Pair of cation and anion lists
     */
    [[nodiscard]] SaltIndexPair cations_and_anions() const
    {
        SaltIndexPair ions;
        for (short z = maxIonCharge; z > 0; --z) {
            for (std::size_t i = 0; i < size(); ++i) {
                const auto ind = static_cast<SaltIndex>(i);
                if (saltData_[i] == T {}) {
                    continue;
                }
                const auto charge = ionStrength(ind);
                if (charge == z) {
                    ions.first.push_back(ind);
                } else if (charge == -z) {
                    ions.second.push_back(ind);
                }
            }
        }
        return ions;
    }

    /*!
     * Return an array with elements converted to new unit
     *
     * @tparam NewSaltUnit New unit for returning array
     * @return Converted array
     */
    template <class NewSaltUnit>
    OPM_HOST_DEVICE SaltArray<T, NewSaltUnit> convert_to() const
    {
        return SaltUnitConverter<SaltUnit, NewSaltUnit>::template convert<T>(*this);
    }

    /*!
     * The ion count goes first, so data written before ions were appended to SaltIndex still
     * loads, with the ions it lacks set to zero.
     */
    template <class Serializer>
    void serializeOp(Serializer& serializer)
    {
        std::size_t count = ncomp;
        serializer(count);
        if (count > ncomp) {
            OPM_THROW(std::runtime_error, "Serialized SaltArray has more ions than SaltIndex");
        }
        for (std::size_t i = 0; i < count; ++i) {
            serializer(saltData_[i]);
        }
        for (std::size_t i = count; i < ncomp; ++i) {
            saltData_[i] = T {};
        }
    }

private:
    std::array<T, ncomp> saltData_{}; ///< Array of salt ion values
}; // class SaltArray

// ------------------------------------------
// Conversion between different salt units
// ------------------------------------------
template <>
struct SaltUnitConverter<SaltMassFraction, SaltMolality>
{
    /*!
     * Convert mass fraction to molality:
     * https://en.wikipedia.org/wiki/Molality#Mass_fraction
     *
     * @tparam T Value type
     * @param saltArray Mass fraction array
     * @return Molality array
     */
    template <class T>
    OPM_HOST_DEVICE static SaltArray<T, SaltMolality> convert(const SaltArray<T, SaltMassFraction>& saltArray)
    {
        SaltArray<T, SaltMolality> molalityArray;
        for (std::size_t i = 0; i < saltArray.size(); ++i) {
            auto sIdx = static_cast<SaltIndex>(i);
            molalityArray[sIdx] = saltArray[sIdx] / saltMolarMass<T>(sIdx);
        }

        T s = 1.0 - saltArray.sum();
        for (auto& elem : molalityArray) {
            elem /= s;
        }

        return molalityArray;
    }
};

template <>
struct SaltUnitConverter<SaltMassFraction, SaltMoleFraction>
{
    /*!
     * Convert mass fraction to mole fraction:
     * https://en.wikipedia.org/wiki/Mass_fraction_(chemistry)#Mole_fraction
     *
     * @tparam T Value type
     * @param saltArray Mass fraction array
     * @return Mole fraction array
     */
    template <class T>
    OPM_HOST_DEVICE static SaltArray<T, SaltMoleFraction> convert(const SaltArray<T, SaltMassFraction>& saltArray)
    {
        SaltArray<T, SaltMoleFraction> moleFracArray;
        T s = 1.0;
        for (std::size_t i = 0; i < saltArray.size(); ++i) {
            auto sIdx = static_cast<SaltIndex>(i);
            moleFracArray[sIdx] = saltArray[sIdx] * 18.01518e-3 / saltMolarMass<T>(sIdx);
            s += moleFracArray[sIdx] - saltArray[sIdx];
        }

        for (auto& elem : moleFracArray) {
            elem /= s;
        }

        return moleFracArray;
    }
};

template <>
struct SaltUnitConverter<SaltMoleFraction, SaltMassFraction>
{
    /*!
     * Convert mole fraction to mass fraction:
     * https://en.wikipedia.org/wiki/Mole_fraction#Mass_fraction
     *
     * @tparam T Value type
     * @param saltArray Mole fraction array
     * @return Mass fraction array
     */
    template <class T>
    OPM_HOST_DEVICE static SaltArray<T, SaltMassFraction> convert(const SaltArray<T, SaltMoleFraction>& saltArray)
    {
        SaltArray<T, SaltMassFraction> massFracArray;
        T s = 18.01518e-3;
        for (std::size_t i = 0; i < saltArray.size(); ++i) {
            auto sIdx = static_cast<SaltIndex>(i);
            massFracArray[sIdx] = saltArray[sIdx] * saltMolarMass<T>(sIdx);
            s += massFracArray[sIdx] - saltArray[sIdx] * 18.01518e-3;
        }

        for (auto& elem : massFracArray) {
            elem /= s;
        }

        return massFracArray;
    }
};

template <>
struct SaltUnitConverter<SaltMoleFraction, SaltMolality>
{
    /*!
     * Convert mole fraction to molality:
     * https://en.wikipedia.org/wiki/Molality#Mole_fraction
     *
     * @tparam T Value type
     * @param saltArray Mole fraction array
     * @return Molality array
     */
    template <class T>
    OPM_HOST_DEVICE static SaltArray<T, SaltMolality> convert(const SaltArray<T, SaltMoleFraction>& saltArray)
    {
        SaltArray<T, SaltMolality> molalityArray;
        for (std::size_t i = 0; i < saltArray.size(); ++i) {
            auto sIdx = static_cast<SaltIndex>(i);
            molalityArray[sIdx] = saltArray[sIdx] / 18.01518e-3;
        }

        T s = 1.0 - saltArray.sum();
        for (auto& elem : molalityArray) {
            elem /= s;
        }

        return molalityArray;
    }
};

template <>
struct SaltUnitConverter<SaltMolality, SaltMoleFraction>
{
    /*!
     * Convert molality to mole fraction:
     * https://en.wikipedia.org/wiki/Molality#Mole_fraction
     *
     * @tparam T Value type
     * @param saltArray Molality array
     * @return Mole fraction array
     */
    template <class T>
    OPM_HOST_DEVICE static SaltArray<T, SaltMoleFraction> convert(const SaltArray<T, SaltMolality>& saltArray)
    {
        SaltArray<T, SaltMoleFraction> moleFracArray;
        T s = 1.0;
        for (std::size_t i = 0; i < saltArray.size(); ++i) {
            auto sIdx = static_cast<SaltIndex>(i);
            moleFracArray[sIdx] = saltArray[sIdx] * 18.01518e-3;
            s += moleFracArray[sIdx];
        }

        for (auto& elem : moleFracArray) {
            elem /= s;
        }

        return moleFracArray;
    }
};

template <>
struct SaltUnitConverter<SaltMolality, SaltMassFraction>
{
    /*!
     * Convert molality to mass fraction:
     * https://en.wikipedia.org/wiki/Molality#Mass_fraction
     *
     * @tparam T Value type
     * @param saltArray Molality array
     * @return Mass fraction array
     */
    template <class T>
    OPM_HOST_DEVICE static SaltArray<T, SaltMassFraction> convert(const SaltArray<T, SaltMolality>& saltArray)
    {
        SaltArray<T, SaltMassFraction> massFracArray;
        T s = 1.0;
        for (std::size_t i = 0; i < saltArray.size(); ++i) {
            auto sIdx = static_cast<SaltIndex>(i);
            massFracArray[sIdx] = saltArray[sIdx] * saltMolarMass<T>(sIdx);
            s += massFracArray[sIdx];
        }

        for (auto& elem : massFracArray) {
            elem /= s;
        }

        return massFracArray;
    }
};

}

#endif //OPM_SALTARRAY_HPP
