#ifndef VESIN_H
#define VESIN_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/// Options for a neighbors list calculation
typedef struct VesinOptions {
    /// Spherical cutoff, only pairs below this cutoff will be included
    double cutoff;
    /// Should the returned neighbors list be a full list (include both `i -> j`
    /// and `j -> i` pairs) or a half list (include only `i -> j`).
    bool full;
    // TODO: sort option?

    /// Should the returned `VesinNeighborsList` contain `shifts`?
    bool return_shifts;
    /// Should the returned `VesinNeighborsList` contain `distances`?
    bool return_distances;
    /// Should the returned `VesinNeighborsList` contain `vector`?
    bool return_vectors;
} VesinOptions;

/// Device on which the data can be
enum VesinDevice {
    /// Unknown device, used for default initialization and to indicate no
    /// allocated data.
    VesinUnknownDevice = 0,
    /// CPU device
    VesinCPU = 1,
};


/// The actual neighbors list
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
typedef struct VesinNeighborsList {
#ifdef __cplusplus
    VesinNeighborsList():
        length(0),
        device(VesinUnknownDevice),
        pairs(nullptr),
        shifts(nullptr),
        distances(nullptr),
        vectors(nullptr)
    {}
#endif

    /// Number of pairs in this neighbors list
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
} VesinNeighborList;

/// Free all allocated memory inside a `VesinNeighborsList`, according the it's
/// `device`.
void vesin_free(VesinNeighborList* neighbors);

/// Compute a neighbors list.
///
/// The data is returned in a `VesinNeighborsList`. For an initial call, the
/// `VesinNeighborsList` should be zero-initialized (or default-initalized in
/// C++). The `VesinNeighborsList` can be re-used across calls to this functions
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
/// @param neighbors non-NULL pointer to `VesinNeighborsList` that will be used
///     to store the computed list of neighbors.
/// @param error_message Pointer to a `char*` that wil be set to the error
///     message if this function fails. This does not need to be freed when no
///     longer needed.
int vesin_neighbors(
    const double (*points)[3],
    size_t n_points,
    const double box[3][3],
    bool periodic,
    VesinDevice device,
    VesinOptions options,
    VesinNeighborList* neighbors,
    const char** error_message
);


#ifdef __cplusplus

} // extern "C"

#endif

#endif
