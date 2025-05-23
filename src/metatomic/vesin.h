/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2024 The METATOMIC-PLUMED team
(see the PEOPLE-METATOMIC file at the root of this folder for a list of names)

See https://docs.metatensor.org/metatomic/ for more information about the
metatomic package that this module allows you to call from PLUMED.

This file is part of METATOMIC-PLUMED module.

The METATOMIC-PLUMED module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The METATOMIC-PLUMED module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the METATOMIC-PLUMED module. If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_metatomic_vesin_h
#define __PLUMED_metatomic_vesin_h
/*INDENT-OFF*/


#include <cstddef>
#include <cstdint>

#if defined(VESIN_SHARED)
    #if defined(VESIN_EXPORTS)
        #if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
            #define VESIN_API __attribute__((visibility("default")))
        #elif defined(_MSC_VER)
            #define VESIN_API __declspec(dllexport)
        #else
            #define VESIN_API
        #endif
    #else
        #if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
            #define VESIN_API __attribute__((visibility("default")))
        #elif defined(_MSC_VER)
            #define VESIN_API __declspec(dllimport)
        #else
            #define VESIN_API
        #endif
    #endif
#else
    #define VESIN_API
#endif

#ifdef __cplusplus
namespace PLMD {
namespace metatomic {
namespace vesin {
extern "C" {
#endif

/// Options for a neighbor list calculation
struct VesinOptions {
    /// Spherical cutoff, only pairs below this cutoff will be included
    double cutoff;
    /// Should the returned neighbor list be a full list (include both `i -> j`
    /// and `j -> i` pairs) or a half list (include only `i -> j`)?
    bool full;
    /// Should the neighbor list be sorted? If yes, the returned pairs will be
    /// sorted using lexicographic order.
    bool sorted;

    /// Should the returned `VesinNeighborList` contain `shifts`?
    bool return_shifts;
    /// Should the returned `VesinNeighborList` contain `distances`?
    bool return_distances;
    /// Should the returned `VesinNeighborList` contain `vector`?
    bool return_vectors;
};

/// Device on which the data can be
enum VesinDevice {
    /// Unknown device, used for default initialization and to indicate no
    /// allocated data.
    VesinUnknownDevice = 0,
    /// CPU device
    VesinCPU = 1,
};


/// The actual neighbor list
///
/// This is organized as a list of pairs, where each pair can contain the
/// following data:
///
/// - indices of the points in the pair;
/// - distance between points in the pair, accounting for periodic boundary
///   conditions;
/// - vector between points in the pair, accounting for periodic boundary
///   conditions;
/// - periodic shift that created the pair. This is only relevant when using
///   periodic boundary conditions, and contains the number of bounding box we
///   need to cross to create the pair. If the positions of the points are `r_i`
///   and `r_j`, the bounding box is described by a matrix of three vectors `H`,
///   and the periodic shift is `S`, the distance vector for a given pair will
///   be given by `r_ij = r_j - r_i + S @ H`.
///
/// Under periodic boundary conditions, two atoms can be part of multiple pairs,
/// each pair having a different periodic shift.
struct VESIN_API VesinNeighborList {
#ifdef __cplusplus
    VesinNeighborList():
        length(0),
        device(VesinUnknownDevice),
        pairs(nullptr),
        shifts(nullptr),
        distances(nullptr),
        vectors(nullptr)
    {}
#endif

    /// Number of pairs in this neighbor list
    size_t length;
    /// Device used for the data allocations
    VesinDevice device;
    /// Array of pairs (storing the indices of the first and second point in the
    /// pair), containing `length` elements.
    size_t (*pairs)[2];
    /// Array of box shifts, one for each `pair`. This is only set if
    /// `options.return_pairs` was `true` during the calculation.
    int32_t (*shifts)[3];
    /// Array of pair distance (i.e. distance between the two points), one for
    /// each pair. This is only set if `options.return_distances` was `true`
    /// during the calculation.
    double *distances;
    /// Array of pair vector (i.e. vector between the two points), one for
    /// each pair. This is only set if `options.return_vector` was `true`
    /// during the calculation.
    double (*vectors)[3];

    // TODO: custom memory allocators?
};

/// Free all allocated memory inside a `VesinNeighborList`, according the it's
/// `device`.
void VESIN_API vesin_free(struct VesinNeighborList* neighbors);

/// Compute a neighbor list.
///
/// The data is returned in a `VesinNeighborList`. For an initial call, the
/// `VesinNeighborList` should be zero-initialized (or default-initalized in
/// C++). The `VesinNeighborList` can be re-used across calls to this functions
/// to re-use memory allocations, and once it is no longer needed, users should
/// call `vesin_free` to release the corresponding memory.
///
/// @param points positions of all points in the system;
/// @param n_points number of elements in the `points` array
/// @param box bounding box for the system. If the system is non-periodic,
///     this is ignored. This should contain the three vectors of the bounding
///     box, one vector per row of the matrix.
/// @param periodic is the system using periodic boundary conditions?
/// @param device device where the `points` and `box` data is allocated.
/// @param options options for the calculation
/// @param neighbors non-NULL pointer to `VesinNeighborList` that will be used
///     to store the computed list of neighbors.
/// @param error_message Pointer to a `char*` that wil be set to the error
///     message if this function fails. This does not need to be freed when no
///     longer needed.
int VESIN_API vesin_neighbors(
    const double (*points)[3],
    size_t n_points,
    const double box[3][3],
    bool periodic,
    VesinDevice device,
    struct VesinOptions options,
    struct VesinNeighborList* neighbors,
    const char** error_message
);


#ifdef __cplusplus

} // extern "C"
} // namespace vesin
} // namespace metatomic
} // namespace PLMD

#endif

#endif
