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
/*INDENT-OFF*/
#include "vesin.h"
// automatically generated 
// vesin version: 0.5.7

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <numeric>
#include <tuple>

#ifndef VESIN_CPU_CELL_LIST_HPP
#define VESIN_CPU_CELL_LIST_HPP

#include <vector>

#include "vesin.h"

#ifndef VESIN_TYPES_HPP
#define VESIN_TYPES_HPP

#include <array>
#include <cassert>
#include <string>

#ifndef VESIN_MATH_HPP
#define VESIN_MATH_HPP

#include <array>
#include <cmath>
#include <stdexcept>

namespace PLMD {
namespace metatomic {
namespace vesin {
struct Vector;

Vector operator*(Vector vector, double scalar);

struct Vector: public std::array<double, 3> {
    double dot(Vector other) const {
        return (*this)[0] * other[0] + (*this)[1] * other[1] + (*this)[2] * other[2];
    }

    double norm() const {
        return std::sqrt(this->dot(*this));
    }

    Vector normalize() const {
        return *this * (1.0 / this->norm());
    }

    Vector cross(Vector other) const {
        return Vector{
            (*this)[1] * other[2] - (*this)[2] * other[1],
            (*this)[2] * other[0] - (*this)[0] * other[2],
            (*this)[0] * other[1] - (*this)[1] * other[0],
        };
    }
};

inline Vector operator+(Vector u, Vector v) {
    return Vector{
        u[0] + v[0],
        u[1] + v[1],
        u[2] + v[2],
    };
}

inline Vector operator-(Vector u, Vector v) {
    return Vector{
        u[0] - v[0],
        u[1] - v[1],
        u[2] - v[2],
    };
}

inline Vector operator*(double scalar, Vector vector) {
    return Vector{
        scalar * vector[0],
        scalar * vector[1],
        scalar * vector[2],
    };
}

inline Vector operator*(Vector vector, double scalar) {
    return Vector{
        scalar * vector[0],
        scalar * vector[1],
        scalar * vector[2],
    };
}

struct Matrix: public std::array<std::array<double, 3>, 3> {
    double determinant() const {
        // clang-format off
        return (*this)[0][0] * ((*this)[1][1] * (*this)[2][2] - (*this)[2][1] * (*this)[1][2])
             - (*this)[0][1] * ((*this)[1][0] * (*this)[2][2] - (*this)[1][2] * (*this)[2][0])
             + (*this)[0][2] * ((*this)[1][0] * (*this)[2][1] - (*this)[1][1] * (*this)[2][0]);
        // clang-format on
    }

    Matrix inverse() const {
        auto det = this->determinant();

        if (std::abs(det) < 1e-30) {
            throw std::runtime_error("this matrix is not invertible");
        }

        auto inverse = Matrix();
        inverse[0][0] = ((*this)[1][1] * (*this)[2][2] - (*this)[2][1] * (*this)[1][2]) / det;
        inverse[0][1] = ((*this)[0][2] * (*this)[2][1] - (*this)[0][1] * (*this)[2][2]) / det;
        inverse[0][2] = ((*this)[0][1] * (*this)[1][2] - (*this)[0][2] * (*this)[1][1]) / det;
        inverse[1][0] = ((*this)[1][2] * (*this)[2][0] - (*this)[1][0] * (*this)[2][2]) / det;
        inverse[1][1] = ((*this)[0][0] * (*this)[2][2] - (*this)[0][2] * (*this)[2][0]) / det;
        inverse[1][2] = ((*this)[1][0] * (*this)[0][2] - (*this)[0][0] * (*this)[1][2]) / det;
        inverse[2][0] = ((*this)[1][0] * (*this)[2][1] - (*this)[2][0] * (*this)[1][1]) / det;
        inverse[2][1] = ((*this)[2][0] * (*this)[0][1] - (*this)[0][0] * (*this)[2][1]) / det;
        inverse[2][2] = ((*this)[0][0] * (*this)[1][1] - (*this)[1][0] * (*this)[0][1]) / det;
        return inverse;
    }
};

inline Vector operator*(Matrix matrix, Vector vector) {
    return Vector{
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    };
}

inline Vector operator*(Vector vector, Matrix matrix) {
    return Vector{
        vector[0] * matrix[0][0] + vector[1] * matrix[1][0] + vector[2] * matrix[2][0],
        vector[0] * matrix[0][1] + vector[1] * matrix[1][1] + vector[2] * matrix[2][1],
        vector[0] * matrix[0][2] + vector[1] * matrix[1][2] + vector[2] * matrix[2][2],
    };
}

} // namespace vesin
} // namespace metatomic
} // namespace PLMD

#endif

namespace PLMD {
namespace metatomic {
namespace vesin {

class BoundingBox {
public:
    BoundingBox(const BoundingBox&) = delete;
    BoundingBox& operator=(const BoundingBox&) = delete;

    BoundingBox(BoundingBox&&) = default;
    BoundingBox& operator=(BoundingBox&&) = default;

    BoundingBox(Matrix matrix, const bool periodic[3]):
        matrix_(matrix),
        periodic_({periodic[0], periodic[1], periodic[2]}),
        max_positions_({-1e300, -1e300, -1e300}),
        min_positions_({1e300, 1e300, 1e300}) {

        // find number of periodic directions and their indices
        int n_periodic = 0;
        int periodic_idx_1 = -1;
        int periodic_idx_2 = -1;
        for (int i = 0; i < 3; ++i) {
            if (periodic_[i]) {
                n_periodic += 1;
                if (periodic_idx_1 == -1) {
                    periodic_idx_1 = i;
                } else if (periodic_idx_2 == -1) {
                    periodic_idx_2 = i;
                }
            }
        }

        // adjust the box matrix to have a simple orthogonal dimension along
        // non-periodic directions
        if (n_periodic == 0) {
            matrix_ = Matrix{
                std::array<double, 3>{1, 0, 0},
                std::array<double, 3>{0, 1, 0},
                std::array<double, 3>{0, 0, 1},
            };
        } else if (n_periodic == 1) {
            assert(periodic_idx_1 != -1);
            // Make the two non-periodic directions orthogonal to the periodic one
            auto a = Vector{matrix_[periodic_idx_1]};
            auto b = Vector{0, 1, 0};
            if (std::abs(a.normalize().dot(b)) > 0.9) {
                b = Vector{0, 0, 1};
            }
            auto c = a.cross(b).normalize();
            b = c.cross(a).normalize();

            // Assign back to the matrix picking the "non-periodic" indices without ifs
            matrix_[(periodic_idx_1 + 1) % 3] = b;
            matrix_[(periodic_idx_1 + 2) % 3] = c;
        } else if (n_periodic == 2) {
            assert(periodic_idx_1 != -1 && periodic_idx_2 != -1);
            // Make the one non-periodic direction orthogonal to the two periodic ones
            auto a = Vector{matrix_[periodic_idx_1]};
            auto b = Vector{matrix_[periodic_idx_2]};
            auto c = a.cross(b).normalize();

            // Assign back to the matrix picking the "non-periodic" index without ifs
            matrix_[(3 - periodic_idx_1 - periodic_idx_2)] = c;
        }

        // precompute the inverse matrix
        auto det = matrix_.determinant();
        if (std::abs(det) < 1e-30) {
            throw std::runtime_error("the box matrix is not invertible");
        }

        this->inverse_ = matrix_.inverse();

        // precompute distances between faces of the bounding box
        auto a = Vector{matrix_[0]};
        auto b = Vector{matrix_[1]};
        auto c = Vector{matrix_[2]};

        // Plans normal vectors
        auto na = b.cross(c).normalize();
        auto nb = c.cross(a).normalize();
        auto nc = a.cross(b).normalize();

        distances_between_faces_ = Vector{
            periodic_[0] ? std::abs(na.dot(a)) : max_positions_[0] - min_positions_[0],
            periodic_[1] ? std::abs(nb.dot(b)) : max_positions_[1] - min_positions_[1],
            periodic_[2] ? std::abs(nc.dot(c)) : max_positions_[2] - min_positions_[2],
        };
    }

    const Matrix& matrix() const {
        return this->matrix_;
    }

    bool periodic(size_t spatial) const {
        return this->periodic_[spatial];
    }

    /// Convert a vector from cartesian coordinates to fractional coordinates
    ///
    /// For non-periodic dimensions, the fractional coordinates are not wrapped
    /// inside [0, 1], but are normalized by the corresponding box length.
    Vector cartesian_to_fractional(Vector cartesian) const {
        auto fractional = cartesian * inverse_;
        if (!periodic_[0]) {
            fractional[0] = (cartesian[0] - min_positions_[0]) / distances_between_faces_[0];
        }

        if (!periodic_[1]) {
            fractional[1] = (cartesian[1] - min_positions_[1]) / distances_between_faces_[1];
        }

        if (!periodic_[2]) {
            fractional[2] = (cartesian[2] - min_positions_[2]) / distances_between_faces_[2];
        }

        return fractional;
    }

    /// Convert a vector from fractional coordinates to cartesian coordinates
    Vector fractional_to_cartesian(Vector fractional) const {
        auto cartesian = fractional * matrix_;

        if (!periodic_[0]) {
            cartesian[0] *= distances_between_faces_[0];
            cartesian[0] += min_positions_[0];
        }

        if (!periodic_[1]) {
            cartesian[1] *= distances_between_faces_[1];
            cartesian[1] += min_positions_[1];
        }

        if (!periodic_[2]) {
            cartesian[2] *= distances_between_faces_[2];
            cartesian[2] += min_positions_[2];
        }

        return cartesian;
    }

    /// Get the three distances between faces of the bounding box
    Vector distances_between_faces() const {
        return distances_between_faces_;
    }

    void make_bounding_for(const double (*points)[3], size_t n_points) {
        // find the min and max coordinates along each axis
        for (size_t i = 0; i < n_points; i++) {
            for (size_t spatial = 0; spatial < 3; spatial++) {
                if (!std::isfinite(points[i][spatial])) {
                    throw std::runtime_error(
                        "point " + std::to_string(i) + " has non-finite coordinate " +
                        "along axis " + std::to_string(spatial) + ": " +
                        std::to_string(points[i][spatial])
                    );
                }

                if (points[i][spatial] < min_positions_[spatial]) {
                    min_positions_[spatial] = points[i][spatial];
                }
                if (points[i][spatial] > max_positions_[spatial]) {
                    max_positions_[spatial] = points[i][spatial];
                }
            }
        }

        for (int dim = 0; dim < 3; dim++) {
            // if all atoms have the same coordinate in this dimension, pretend
            // that the bounding box is at least 1 unit wide to avoid numerical issues
            if (max_positions_[dim] - min_positions_[dim] < 1e-6) {
                max_positions_[dim] = min_positions_[dim] + 1;
            }

            if (!periodic_[dim]) {
                // add a 1% margin to make sure all points are strictly inside the
                // bounding box
                distances_between_faces_[dim] = max_positions_[dim] * 1.01 - min_positions_[dim];
            }
        }
    }

private:
    Matrix matrix_;
    std::array<bool, 3> periodic_;

    Matrix inverse_;
    Vector min_positions_;
    Vector max_positions_;
    Vector distances_between_faces_;
};

/// A cell shift represents the displacement along cell axis between the actual
/// position of an atom and a periodic image of this atom.
///
/// The cell shift can be used to reconstruct the vector between two points,
/// wrapped inside the unit cell.
struct CellShift: public std::array<int32_t, 3> {
    /// Compute the shift vector in cartesian coordinates, using the given cell
    /// matrix (stored in row major order).
    Vector cartesian(const BoundingBox& box) const {
        assert(box.periodic(0) || (*this)[0] == 0);
        assert(box.periodic(1) || (*this)[1] == 0);
        assert(box.periodic(2) || (*this)[2] == 0);

        auto vector = Vector{
            static_cast<double>((*this)[0]),
            static_cast<double>((*this)[1]),
            static_cast<double>((*this)[2]),
        };
        return vector * box.matrix();
    }
};

inline CellShift operator+(CellShift a, CellShift b) {
    return CellShift{
        a[0] + b[0],
        a[1] + b[1],
        a[2] + b[2],
    };
}

inline CellShift operator-(CellShift a, CellShift b) {
    return CellShift{
        a[0] - b[0],
        a[1] - b[1],
        a[2] - b[2],
    };
}

} // namespace vesin
} // namespace metatomic
} // namespace PLMD

#endif

namespace PLMD {
namespace metatomic {
namespace vesin {
namespace cpu {

void free_neighbors(VesinNeighborList& neighbors);

void neighbors(
    const Vector* points,
    size_t n_points,
    BoundingBox box,
    VesinOptions options,
    VesinNeighborList& neighbors
);

/// The cell list is used to sort atoms inside bins/cells.
///
/// The list of potential pairs is then constructed by looking through all
/// neighboring cells (the number of cells to search depends on the cutoff and
/// the size of the cells) for each atom to create pair candidates.
class CellList {
public:
    /// Create a new `CellList` for the given bounding box and cutoff,
    /// determining all required parameters.
    CellList(BoundingBox box, double cutoff);

    /// Add a single point to the cell list at the given `position`. The point
    /// is uniquely identified by its `index`.
    void add_point(size_t index, Vector position);

    /// Iterate over all possible pairs, calling the given callback every time
    template <typename Function>
    void foreach_pair(Function callback);

private:
    /// How many cells do we need to look at when searching neighbors to include
    /// all neighbors below cutoff
    std::array<int32_t, 3> n_search_;

    /// the cells themselves are a list of points & corresponding
    /// shift to place the point inside the cell
    struct Point {
        size_t index;
        CellShift shift;
    };
    struct Cell: public std::vector<Point> {};

    // raw data for the cells
    std::vector<Cell> cells_;
    // shape of the cell array
    std::array<size_t, 3> cells_shape_;

    BoundingBox box_;

    Cell& get_cell(std::array<int32_t, 3> index);
};

/// Wrapper around `VesinNeighborList` that behaves like a std::vector,
/// automatically growing memory allocations.
class GrowableNeighborList {
public:
    VesinNeighborList& neighbors;
    size_t capacity;
    VesinOptions options;

    size_t length() const {
        return neighbors.length;
    }

    void increment_length() {
        neighbors.length += 1;
    }

    void set_pair(size_t index, size_t first, size_t second);
    void set_shift(size_t index, PLMD::metatomic::vesin::CellShift shift);
    void set_distance(size_t index, double distance);
    void set_vector(size_t index, PLMD::metatomic::vesin::Vector vector);

    // reset length to 0, and allocate/deallocate members of
    // `neighbors` according to `options`
    void reset();

    // allocate more memory & update capacity
    void grow();

    // sort the pairs currently in the neighbor list
    void sort();
};

} // namespace cpu
} // namespace vesin
} // namespace metatomic
} // namespace PLMD

#endif

using namespace PLMD::metatomic::vesin;
using namespace PLMD::metatomic::vesin::cpu;

void PLMD::metatomic::vesin::cpu::neighbors(
    const Vector* points,
    size_t n_points,
    BoundingBox box,
    VesinOptions options,
    VesinNeighborList& raw_neighbors
) {
    if (options.algorithm == VesinAutoAlgorithm || options.algorithm == VesinCellList) {
        // all good, this is the only thing we implement
    } else {
        throw std::runtime_error("only VesinAutoAlgorithm and VesinCellList are supported on CPU");
    }

    auto cell_list = CellList(std::move(box), options.cutoff);

    for (size_t i = 0; i < n_points; i++) {
        cell_list.add_point(i, points[i]);
    }

    auto cutoff2 = options.cutoff * options.cutoff;

    // the cell list creates too many pairs, we only need to keep the
    // one where the distance is actually below the cutoff
    auto neighbors = GrowableNeighborList{raw_neighbors, raw_neighbors.length, options};
    neighbors.reset();

    cell_list.foreach_pair([&](size_t first, size_t second, CellShift shift) {
        if (!options.full) {
            // filter out some pairs for half neighbor lists
            if (first > second) {
                return;
            }

            if (first == second) {
                // When creating pairs between a point and one of its periodic
                // images, the code generate multiple redundant pairs (e.g. with
                // shifts 0 1 1 and 0 -1 -1); and we want to only keep one of
                // these.
                if (shift[0] + shift[1] + shift[2] < 0) {
                    // drop shifts on the negative half-space
                    return;
                }

                if ((shift[0] + shift[1] + shift[2] == 0) && (shift[2] < 0 || (shift[2] == 0 && shift[1] < 0))) {
                    // drop shifts in the negative half plane or the negative
                    // shift[1] axis. See below for a graphical representation:
                    // we are keeping the shifts indicated with `O` and dropping
                    // the ones indicated with `X`
                    //
                    //  O O O │ O O O
                    //  O O O │ O O O
                    //  O O O │ O O O
                    // ─X─X─X─┼─O─O─O─
                    //  X X X │ X X X
                    //  X X X │ X X X
                    //  X X X │ X X X
                    return;
                }
            }
        }

        auto vector = points[second] - points[first] + shift.cartesian(box);
        auto distance2 = vector.dot(vector);

        if (distance2 < cutoff2) {
            auto index = neighbors.length();
            neighbors.set_pair(index, first, second);

            if (options.return_shifts) {
                neighbors.set_shift(index, shift);
            }

            if (options.return_distances) {
                neighbors.set_distance(index, std::sqrt(distance2));
            }

            if (options.return_vectors) {
                neighbors.set_vector(index, vector);
            }

            neighbors.increment_length();
        }
    });

    if (options.sorted) {
        neighbors.sort();
    }
}

/* ========================================================================== */

/// Maximal number of cells, we need to use this to prevent having too many
/// cells with a small bounding box and a large cutoff
#define MAX_NUMBER_OF_CELLS 1e5

/// Function to compute both quotient and remainder of the division of a by b.
/// This function follows Python convention, making sure the remainder have the
/// same sign as `b`.
static std::tuple<int32_t, int32_t> divmod(int32_t a, size_t b) {
    assert(b < (std::numeric_limits<int32_t>::max()));
    auto b_32 = static_cast<int32_t>(b);
    auto quotient = a / b_32;
    auto remainder = a % b_32;
    if (remainder < 0) {
        remainder += b_32;
        quotient -= 1;
    }
    return std::make_tuple(quotient, remainder);
}

/// Apply the `divmod` function to three components at the time
static std::tuple<std::array<int32_t, 3>, std::array<int32_t, 3>>
divmod(std::array<int32_t, 3> a, std::array<size_t, 3> b) {
    auto [qx, rx] = divmod(a[0], b[0]);
    auto [qy, ry] = divmod(a[1], b[1]);
    auto [qz, rz] = divmod(a[2], b[2]);
    return std::make_tuple(
        std::array<int32_t, 3>{qx, qy, qz},
        std::array<int32_t, 3>{rx, ry, rz}
    );
}

CellList::CellList(BoundingBox box, double cutoff):
    n_search_({0, 0, 0}),
    cells_shape_({0, 0, 0}),
    box_(std::move(box)) {
    auto distances_between_faces = box_.distances_between_faces();

    auto n_cells = Vector{
        std::clamp(std::trunc(distances_between_faces[0] / cutoff), 1.0, HUGE_VAL),
        std::clamp(std::trunc(distances_between_faces[1] / cutoff), 1.0, HUGE_VAL),
        std::clamp(std::trunc(distances_between_faces[2] / cutoff), 1.0, HUGE_VAL),
    };

    assert(std::isfinite(n_cells[0]) && std::isfinite(n_cells[1]) && std::isfinite(n_cells[2]));

    // limit memory consumption by ensuring we have less than `MAX_N_CELLS`
    // cells to look though
    auto n_cells_total = n_cells[0] * n_cells[1] * n_cells[2];
    if (n_cells_total > MAX_NUMBER_OF_CELLS) {
        // set the total number of cells close to MAX_N_CELLS, while keeping
        // roughly the ratio of cells in each direction
        auto ratio_x_y = n_cells[0] / n_cells[1];
        auto ratio_y_z = n_cells[1] / n_cells[2];

        n_cells[2] = std::trunc(std::cbrt(MAX_NUMBER_OF_CELLS / (ratio_x_y * ratio_y_z * ratio_y_z)));
        n_cells[1] = std::trunc(ratio_y_z * n_cells[2]);
        n_cells[0] = std::trunc(ratio_x_y * n_cells[1]);
    }

    // number of cells to search in each direction to make sure all possible
    // pairs below the cutoff are accounted for.
    this->n_search_ = std::array<int32_t, 3>{
        static_cast<int32_t>(std::ceil(cutoff * n_cells[0] / distances_between_faces[0])),
        static_cast<int32_t>(std::ceil(cutoff * n_cells[1] / distances_between_faces[1])),
        static_cast<int32_t>(std::ceil(cutoff * n_cells[2] / distances_between_faces[2])),
    };

    this->cells_shape_ = std::array<size_t, 3>{
        static_cast<size_t>(n_cells[0]),
        static_cast<size_t>(n_cells[1]),
        static_cast<size_t>(n_cells[2]),
    };

    for (size_t spatial = 0; spatial < 3; spatial++) {
        if (n_search_[spatial] < 1) {
            n_search_[spatial] = 1;
        }
    }

    this->cells_.resize(cells_shape_[0] * cells_shape_[1] * cells_shape_[2]);
}

void CellList::add_point(size_t index, Vector position) {
    auto fractional = box_.cartesian_to_fractional(position);

    // find the cell in which this atom should go
    auto cell_index = std::array<int32_t, 3>{
        static_cast<int32_t>(std::floor(fractional[0] * static_cast<double>(cells_shape_[0]))),
        static_cast<int32_t>(std::floor(fractional[1] * static_cast<double>(cells_shape_[1]))),
        static_cast<int32_t>(std::floor(fractional[2] * static_cast<double>(cells_shape_[2]))),
    };

    // deal with pbc by wrapping the atom inside if it was outside of the cell
    CellShift shift;
    for (size_t spatial = 0; spatial < 3; spatial++) {
        auto result = divmod(cell_index[spatial], cells_shape_[spatial]);
        shift[spatial] = std::get<0>(result);
        cell_index[spatial] = std::get<1>(result);

        assert(box_.periodic(spatial) || shift[spatial] == 0);
    }

    this->get_cell(cell_index).emplace_back(Point{index, shift});
}

// clang-format off
template <typename Function>
void CellList::foreach_pair(Function callback) {
    for (int32_t cell_i_x=0; cell_i_x<static_cast<int32_t>(cells_shape_[0]); cell_i_x++) {
    for (int32_t cell_i_y=0; cell_i_y<static_cast<int32_t>(cells_shape_[1]); cell_i_y++) {
    for (int32_t cell_i_z=0; cell_i_z<static_cast<int32_t>(cells_shape_[2]); cell_i_z++) {
        const auto& current_cell = this->get_cell({cell_i_x, cell_i_y, cell_i_z});
        // look through each neighboring cell
        for (int32_t delta_x=-n_search_[0]; delta_x<=n_search_[0]; delta_x++) {
        for (int32_t delta_y=-n_search_[1]; delta_y<=n_search_[1]; delta_y++) {
        for (int32_t delta_z=-n_search_[2]; delta_z<=n_search_[2]; delta_z++) {
            auto cell_i = std::array<int32_t, 3>{
                cell_i_x + delta_x,
                cell_i_y + delta_y,
                cell_i_z + delta_z,
            };

            // shift vector from one cell to the other and index of
            // the neighboring cell
            auto [cell_shift, neighbor_cell_i] = divmod(cell_i, cells_shape_);

            for (const auto& atom_i: current_cell) {
                for (const auto& atom_j: this->get_cell(neighbor_cell_i)) {
                    auto shift = CellShift{cell_shift} + atom_i.shift - atom_j.shift;
                    auto shift_is_zero = shift[0] == 0 && shift[1] == 0 && shift[2] == 0;

                    if ((shift[0] != 0 && !box_.periodic(0)) ||
                        (shift[1] != 0 && !box_.periodic(1)) ||
                        (shift[2] != 0 && !box_.periodic(2)))
                    {
                        // do not create pairs crossing the periodic
                        // boundaries in a non-periodic box
                        continue;
                    }

                    if (atom_i.index == atom_j.index && shift_is_zero) {
                        // only create pairs with the same atom twice if the
                        // pair spans more than one bounding box
                        continue;
                    }

                    callback(atom_i.index, atom_j.index, shift);
                }
            } // loop over atoms in current neighbor cells
        }}}
    }}} // loop over neighboring cells
}

CellList::Cell& CellList::get_cell(std::array<int32_t, 3> index) {
    size_t linear_index = (cells_shape_[0] * cells_shape_[1] * index[2])
                        + (cells_shape_[0] * index[1])
                        + index[0];
    return cells_[linear_index];
}
// clang-format on

/* ========================================================================== */

void GrowableNeighborList::set_pair(size_t index, size_t first, size_t second) {
    if (index >= this->capacity) {
        this->grow();
    }

    this->neighbors.pairs[index][0] = first;
    this->neighbors.pairs[index][1] = second;
}

void GrowableNeighborList::set_shift(size_t index, PLMD::metatomic::vesin::CellShift shift) {
    if (index >= this->capacity) {
        this->grow();
    }

    this->neighbors.shifts[index][0] = shift[0];
    this->neighbors.shifts[index][1] = shift[1];
    this->neighbors.shifts[index][2] = shift[2];
}

void GrowableNeighborList::set_distance(size_t index, double distance) {
    if (index >= this->capacity) {
        this->grow();
    }

    this->neighbors.distances[index] = distance;
}

void GrowableNeighborList::set_vector(size_t index, PLMD::metatomic::vesin::Vector vector) {
    if (index >= this->capacity) {
        this->grow();
    }

    this->neighbors.vectors[index][0] = vector[0];
    this->neighbors.vectors[index][1] = vector[1];
    this->neighbors.vectors[index][2] = vector[2];
}

template <typename scalar_t, size_t N>
using array_ptr = scalar_t (*)[N];

template <typename scalar_t, size_t N>
static array_ptr<scalar_t, N> alloc(array_ptr<scalar_t, N> ptr, size_t size, size_t new_size) {
    auto* new_ptr = reinterpret_cast<scalar_t(*)[N]>(std::realloc(ptr, new_size * sizeof(scalar_t[N])));

    if (new_ptr == nullptr) {
        return nullptr;
    }

#ifndef NDEBUG
    // initialize with a bit pattern that maps to NaN for double
    std::memset(new_ptr + size, 0b11111111, (new_size - size) * sizeof(scalar_t[N]));
#endif

    return new_ptr;
}

template <typename scalar_t>
static scalar_t* alloc(scalar_t* ptr, size_t size, size_t new_size) {
    auto* new_ptr = reinterpret_cast<scalar_t*>(std::realloc(ptr, new_size * sizeof(scalar_t)));

    if (new_ptr == nullptr) {
        return nullptr;
    }

#ifndef NDEBUG
    // initialize with a bit pattern that maps to NaN for double
    std::memset(new_ptr + size, 0b11111111, (new_size - size) * sizeof(scalar_t));
#endif

    return new_ptr;
}

void GrowableNeighborList::grow() {
    auto new_size = neighbors.length * 2;
    if (new_size == 0) {
        new_size = 1;
    }

    auto* new_pairs = alloc<size_t, 2>(neighbors.pairs, neighbors.length, new_size);

    int32_t (*new_shifts)[3] = nullptr;
    if (options.return_shifts) {
        new_shifts = alloc<int32_t, 3>(neighbors.shifts, neighbors.length, new_size);
    }

    double* new_distances = nullptr;
    if (options.return_distances) {
        new_distances = alloc<double>(neighbors.distances, neighbors.length, new_size);
    }

    double (*new_vectors)[3] = nullptr;
    if (options.return_vectors) {
        new_vectors = alloc<double, 3>(neighbors.vectors, neighbors.length, new_size);
    }

    if (
        (new_pairs == nullptr) ||
        (options.return_shifts && new_shifts == nullptr) ||
        (options.return_distances && new_distances == nullptr) ||
        (options.return_vectors && new_vectors == nullptr)
    ) {
        std::free(new_pairs);
        std::free(new_shifts);
        std::free(new_distances);
        std::free(new_vectors);
        throw std::runtime_error("could not allocate memory for growing neighbor list");
    }

    this->neighbors.pairs = new_pairs;
    this->neighbors.shifts = new_shifts;
    this->neighbors.distances = new_distances;
    this->neighbors.vectors = new_vectors;

    this->capacity = new_size;
}

void GrowableNeighborList::reset() {
#ifndef NDEBUG
    auto size = this->neighbors.length;
    // set all allocated data to a bit pattern that maps to NaN for double
    std::memset(this->neighbors.pairs, 0b11111111, size * sizeof(size_t[2]));

    if (this->neighbors.shifts != nullptr) {
        std::memset(this->neighbors.shifts, 0b11111111, size * sizeof(int32_t[3]));
    }

    if (this->neighbors.distances != nullptr) {
        std::memset(this->neighbors.distances, 0b11111111, size * sizeof(double));
    }

    if (this->neighbors.vectors != nullptr) {
        std::memset(this->neighbors.vectors, 0b11111111, size * sizeof(double[3]));
    }
#endif

    // reset length (but keep the capacity where it's at)
    this->neighbors.length = 0;

    // allocate/deallocate pointers as required
    auto* shifts = this->neighbors.shifts;
    if (this->options.return_shifts && shifts == nullptr) {
        shifts = alloc<int32_t, 3>(shifts, 0, capacity);
    } else if (!this->options.return_shifts && shifts != nullptr) {
        std::free(shifts);
        shifts = nullptr;
    }

    auto* distances = this->neighbors.distances;
    if (this->options.return_distances && distances == nullptr) {
        distances = alloc<double>(distances, 0, capacity);
    } else if (!this->options.return_distances && distances != nullptr) {
        std::free(distances);
        distances = nullptr;
    }

    auto* vectors = this->neighbors.vectors;
    if (this->options.return_vectors && vectors == nullptr) {
        vectors = alloc<double, 3>(vectors, 0, capacity);
    } else if (!this->options.return_vectors && vectors != nullptr) {
        std::free(vectors);
        vectors = nullptr;
    }

    this->neighbors.shifts = shifts;
    this->neighbors.distances = distances;
    this->neighbors.vectors = vectors;
}

void GrowableNeighborList::sort() {
    if (this->length() == 0) {
        return;
    }

    // step 1: sort an array of indices, comparing the pairs at the indices
    auto indices = std::vector<int64_t>(this->length(), 0);
    std::iota(std::begin(indices), std::end(indices), 0);

    struct compare_pairs {
        compare_pairs(size_t (*pairs_)[2]):
            pairs(pairs_) {}

        bool operator()(int64_t a, int64_t b) const {
            return pairs[a][0] < pairs[b][0];
        }

        size_t (*pairs)[2];
    };

    std::sort(
        std::begin(indices),
        std::end(indices),
        compare_pairs(this->neighbors.pairs)
    );

    // step 2: move all data according to the sorted indices.
    auto* sorted_pairs = alloc<size_t, 2>(nullptr, 0, this->capacity);

    int32_t (*sorted_shifts)[3] = nullptr;
    if (options.return_shifts) {
        sorted_shifts = alloc<int32_t, 3>(nullptr, 0, this->capacity);
    }

    double* sorted_distances = nullptr;
    if (options.return_distances) {
        sorted_distances = alloc<double>(nullptr, 0, this->capacity);
    }

    double (*sorted_vectors)[3] = nullptr;
    if (options.return_vectors) {
        sorted_vectors = alloc<double, 3>(nullptr, 0, this->capacity);
    }

    if (
        (sorted_pairs == nullptr) ||
        (options.return_shifts && sorted_shifts == nullptr) ||
        (options.return_distances && sorted_distances == nullptr) ||
        (options.return_vectors && sorted_vectors == nullptr)
    ) {
        std::free(sorted_pairs);
        std::free(sorted_shifts);
        std::free(sorted_distances);
        std::free(sorted_vectors);
        throw std::runtime_error("could not allocate memory for sorting neighbor list");
    }

    for (size_t i = 0; i < this->neighbors.length; i++) {
        auto from = static_cast<size_t>(indices[i]);
        sorted_pairs[i][0] = this->neighbors.pairs[from][0];
        sorted_pairs[i][1] = this->neighbors.pairs[from][1];

        if (options.return_shifts) {
            sorted_shifts[i][0] = this->neighbors.shifts[from][0];
            sorted_shifts[i][1] = this->neighbors.shifts[from][1];
            sorted_shifts[i][2] = this->neighbors.shifts[from][2];
        }

        if (options.return_distances) {
            sorted_distances[i] = this->neighbors.distances[from];
        }

        if (options.return_vectors) {
            sorted_vectors[i][0] = this->neighbors.vectors[from][0];
            sorted_vectors[i][1] = this->neighbors.vectors[from][1];
            sorted_vectors[i][2] = this->neighbors.vectors[from][2];
        }
    }

    std::free(this->neighbors.pairs);
    this->neighbors.pairs = sorted_pairs;

    if (options.return_shifts) {
        std::free(this->neighbors.shifts);
        this->neighbors.shifts = sorted_shifts;
    }

    if (options.return_distances) {
        std::free(this->neighbors.distances);
        this->neighbors.distances = sorted_distances;
    }

    if (options.return_vectors) {
        std::free(this->neighbors.vectors);
        this->neighbors.vectors = sorted_vectors;
    }
}

void PLMD::metatomic::vesin::cpu::free_neighbors(VesinNeighborList& neighbors) {
    assert(neighbors.device.type == VesinCPU);

    std::free(neighbors.pairs);
    std::free(neighbors.shifts);
    std::free(neighbors.vectors);
    std::free(neighbors.distances);
}
#include <stdexcept>

#ifndef VESIN_CUDA_HPP
#define VESIN_CUDA_HPP


namespace PLMD {
namespace metatomic {
namespace vesin {
namespace cuda {

#ifndef VESIN_CUDA_AT_LEAST_PAIRS_PER_POINT
/// Default value for the number of pairs per points in the CUDA implementation.
/// Unless `VESIN_CUDA_MAX_PAIRS_PER_POINT` is set in the environement, the
/// maximal number of pairs is `n_points *
/// max(VESIN_CUDA_AT_LEAST_PAIRS_PER_POINT, cutoff^3)`. This can be overriden
/// at compile time.
#define VESIN_CUDA_AT_LEAST_PAIRS_PER_POINT 128
#endif

/// @brief Buffers for cell list-based neighbor search
struct CellListBuffers {
    size_t max_points = 0; // Capacity for point-related arrays
    size_t max_cells = 0;  // Capacity for cell-related arrays

    // Per-particle arrays
    int32_t* cell_indices = nullptr;    // [max_points] linear cell index per particle
    int32_t* particle_shifts = nullptr; // [max_points * 3] shift applied to wrap into cell

    // Per-cell arrays
    int32_t* cell_counts = nullptr;  // [max_cells] number of particles in each cell
    int32_t* cell_starts = nullptr;  // [max_cells] starting index in sorted arrays
    int32_t* cell_offsets = nullptr; // [max_cells] working copy for scatter

    // Sorted particle data (for coalesced memory access)
    double* sorted_positions = nullptr;     // [max_points * 3]
    int32_t* sorted_indices = nullptr;      // [max_points] original particle indices
    int32_t* sorted_shifts = nullptr;       // [max_points * 3] shifts for sorted particles
    int32_t* sorted_cell_indices = nullptr; // [max_points] cell indices in sorted order

    // Cell grid parameters (computed on device)
    double* inv_box = nullptr;        // [9] inverse box matrix
    int32_t* n_cells = nullptr;       // [3] number of cells in each direction
    int32_t* n_search = nullptr;      // [3] search range in each direction
    int32_t* n_cells_total = nullptr; // [1] total number of cells

    double* bounding_min = nullptr; // [3] per-dimension min for non-periodic axes
    double* bounding_max = nullptr; // [3] per-dimension max for non-periodic axes
};

struct CudaNeighborListExtras {
    size_t* length_ptr = nullptr;      // GPU-side counter
    size_t capacity = 0;               // Current capacity per device
    size_t max_pairs = 0;              // Maximum number of pairs that can be stored; depends on VESIN_CUDA_MAX_PAIRS_PER_POINT
    int32_t* cell_check_ptr = nullptr; // GPU-side status code for checking cell
    int32_t* overflow_flag = nullptr;  // GPU-side flag to detect overflow of pair buffers
    int32_t allocated_device_id = -1;  // which device are we currently allocated on

    // Pinned host memory for async D2H copy (Approach 2)
    size_t* pinned_length_ptr = nullptr;

    // Cell list buffers (allocated on demand for large systems)
    CellListBuffers cell_list;

    // Buffers for optimized brute force kernels
    double* box_diag = nullptr;      // [3] diagonal elements for orthogonal boxes
    double* inv_box_brute = nullptr; // [9] inverse box matrix for general boxes

    // Temporary buffers for on-device sorting
    size_t* sort_pairs_tmp = nullptr;     // [sort_capacity * 2]
    int32_t* sort_shifts_tmp = nullptr;   // [sort_capacity * 3]
    double* sort_distances_tmp = nullptr; // [sort_capacity]
    double* sort_vectors_tmp = nullptr;   // [sort_capacity * 3]
    size_t sort_capacity = 0;

    ~CudaNeighborListExtras();
};

/// @brief Frees GPU memory associated with a VesinNeighborList.
///
/// This function should be called to release all CUDA-allocated memory
/// tied to the given neighbor list. It does not delete the structure itself,
/// only the device-side memory buffers.
///
/// @param neighbors Reference to the VesinNeighborList to clean up.
void free_neighbors(VesinNeighborList& neighbors);

/// @brief Computes the neighbor list on the GPU.
///
/// This function only works under Minimum Image Convention for now.
///
/// This function generates a neighbor list for a set of points within a
/// periodic simulation box using GPU acceleration. The output is stored in a
/// `VesinNeighborList` structure, which must be initialized for GPU usage.
///
/// @param points Pointer to an array of 3D points (shape: [n_points][3]).
/// @param n_points Number of points (atoms, particles, etc.).
/// @param box 3×3 matrix defining the bounding box of the system.
/// @param periodic Array of three booleans indicating periodicity in each dimension.
/// @param options Struct holding parameters such as cutoff, symmetry, etc.
/// @param neighbors Output neighbor list (device memory will be allocated as
/// needed).
void neighbors(
    const double (*points)[3],
    size_t n_points,
    const double box[3][3],
    const bool periodic[3],
    VesinOptions options,
    VesinNeighborList& neighbors
);

/// Get the `CudaNeighborListExtras` stored inside `VesinNeighborList`'s opaque pointer
CudaNeighborListExtras* get_cuda_extras(VesinNeighborList* neighbors);

} // namespace cuda
} // namespace vesin
} // namespace metatomic
} // namespace PLMD

#endif // VESIN_CUDA_HPP

void PLMD::metatomic::vesin::cuda::free_neighbors(VesinNeighborList& neighbors) {
    throw std::runtime_error("CUDA neighbor list generation is not included in this build of vesin");
}

void PLMD::metatomic::vesin::cuda::neighbors(
    const double (*points)[3],
    size_t n_points,
    const double box[3][3],
    const bool periodic[3],
    VesinOptions options,
    VesinNeighborList& neighbors
) {
    throw std::runtime_error("CUDA neighbor list generation is not included in this build of vesin");
}
#include <cstring>
#include <iostream>
#include <string>


// used to store dynamically allocated error messages before giving a pointer
// to them back to the user
thread_local std::string LAST_ERROR;

extern "C" int vesin_neighbors(
    const double (*points)[3],
    size_t n_points,
    const double box[3][3],
    const bool periodic[3],
    VesinDevice device,
    VesinOptions options,
    VesinNeighborList* neighbors,
    const char** error_message
) {
    if (error_message == nullptr) {
        return EXIT_FAILURE;
    }

    if (points == nullptr) {
        *error_message = "`points` can not be a NULL pointer";
        return EXIT_FAILURE;
    }

    if (box == nullptr) {
        *error_message = "`cell` can not be a NULL pointer";
        return EXIT_FAILURE;
    }

    if (neighbors == nullptr) {
        *error_message = "`neighbors` can not be a NULL pointer";
        return EXIT_FAILURE;
    }

    if (!std::isfinite(options.cutoff) || options.cutoff <= 0) {
        *error_message = "cutoff must be a finite, positive number";
        return EXIT_FAILURE;
    }

    if (options.cutoff <= 1e-6) {
        *error_message = "cutoff is too small";
        return EXIT_FAILURE;
    }

    if (neighbors->device.type != VesinUnknownDevice && neighbors->device.type != device.type) {
        *error_message = "`neighbors` device and data `device` do not match, free the neighbors first";
        return EXIT_FAILURE;
    }

    if (device.type == VesinUnknownDevice) {
        *error_message = "got an unknown device type";
        return EXIT_FAILURE;
    }

    if (neighbors->device.type == VesinUnknownDevice) {
        // initialize the device
        neighbors->device = device;
    } else if (neighbors->device.type != device.type) {
        *error_message = "`neighbors.device` and `device` do not match, free the neighbors first";
        return EXIT_FAILURE;
    }

    try {
        if (device.type == VesinCPU) {
            auto matrix = PLMD::metatomic::vesin::Matrix{{{
                {{box[0][0], box[0][1], box[0][2]}},
                {{box[1][0], box[1][1], box[1][2]}},
                {{box[2][0], box[2][1], box[2][2]}},
            }}};

            auto box = PLMD::metatomic::vesin::BoundingBox(matrix, periodic);
            box.make_bounding_for(points, n_points);

            PLMD::metatomic::vesin::cpu::neighbors(
                reinterpret_cast<const PLMD::metatomic::vesin::Vector*>(points),
                n_points,
                std::move(box),
                options,
                *neighbors
            );
        } else if (device.type == VesinCUDA) {
            PLMD::metatomic::vesin::cuda::neighbors(
                points,
                n_points,
                box,
                periodic,
                options,
                *neighbors
            );
        } else {
            throw std::runtime_error("unknown device " + std::to_string(device.type));
        }
    } catch (const std::bad_alloc&) {
        LAST_ERROR = "failed to allocate memory";
        *error_message = LAST_ERROR.c_str();
        return EXIT_FAILURE;
    } catch (const std::exception& e) {
        LAST_ERROR = e.what();
        *error_message = LAST_ERROR.c_str();
        return EXIT_FAILURE;
    } catch (...) {
        *error_message = "fatal error: unknown type thrown as exception";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

extern "C" void vesin_free(VesinNeighborList* neighbors) {
    if (neighbors == nullptr) {
        return;
    }

    try {
        if (neighbors->device.type == VesinUnknownDevice) {
            // nothing to do
        } else if (neighbors->device.type == VesinCPU) {
            PLMD::metatomic::vesin::cpu::free_neighbors(*neighbors);
        } else if (neighbors->device.type == VesinCUDA) {
            PLMD::metatomic::vesin::cuda::free_neighbors(*neighbors);
        } else {
            throw std::runtime_error("unknown device " + std::to_string(neighbors->device.type) + " when freeing memory");
        }
    } catch (const std::exception& e) {
        std::cerr << "error in vesin_free: " << e.what() << std::endl;
        return;
    } catch (...) {
        std::cerr << "fatal error in vesin_free, unknown type thrown as exception" << std::endl;
        return;
    }

    std::memset(neighbors, 0, sizeof(VesinNeighborList));
}
