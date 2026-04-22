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
// automatically generated 
// vesin version: 0.5.3

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <new>
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
    BoundingBox(Matrix matrix, const bool periodic[3]):
        matrix_(matrix),
        periodic_({periodic[0], periodic[1], periodic[2]}) {

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
    }

    const Matrix& matrix() const {
        return this->matrix_;
    }

    bool periodic(size_t spatial) const {
        return this->periodic_[spatial];
    }

    /// Convert a vector from cartesian coordinates to fractional coordinates
    Vector cartesian_to_fractional(Vector cartesian) const {
        return cartesian * inverse_;
    }

    /// Convert a vector from fractional coordinates to cartesian coordinates
    Vector fractional_to_cartesian(Vector fractional) const {
        return fractional * matrix_;
    }

    /// Get the three distances between faces of the bounding box
    Vector distances_between_faces() const {
        auto a = Vector{matrix_[0]};
        auto b = Vector{matrix_[1]};
        auto c = Vector{matrix_[2]};

        // Plans normal vectors
        auto na = b.cross(c).normalize();
        auto nb = c.cross(a).normalize();
        auto nc = a.cross(b).normalize();

        return Vector{
            std::abs(na.dot(a)),
            std::abs(nb.dot(b)),
            std::abs(nc.dot(c)),
        };
    }

private:
    Matrix matrix_;
    Matrix inverse_;
    std::array<bool, 3> periodic_;
};

/// A cell shift represents the displacement along cell axis between the actual
/// position of an atom and a periodic image of this atom.
///
/// The cell shift can be used to reconstruct the vector between two points,
/// wrapped inside the unit cell.
struct CellShift: public std::array<int32_t, 3> {
    /// Compute the shift vector in cartesian coordinates, using the given cell
    /// matrix (stored in row major order).
    Vector cartesian(Matrix cell) const {
        auto vector = Vector{
            static_cast<double>((*this)[0]),
            static_cast<double>((*this)[1]),
            static_cast<double>((*this)[2]),
        };
        return vector * cell;
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
    BoundingBox cell,
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
    BoundingBox cell,
    VesinOptions options,
    VesinNeighborList& raw_neighbors
) {
    if (options.algorithm == VesinAutoAlgorithm || options.algorithm == VesinCellList) {
        // all good, this is the only thing we implement
    } else {
        throw std::runtime_error("only VesinAutoAlgorithm and VesinCellList are supported on CPU");
    }

    auto cell_list = CellList(cell, options.cutoff);

    for (size_t i = 0; i < n_points; i++) {
        cell_list.add_point(i, points[i]);
    }

    auto cell_matrix = cell.matrix();
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

        auto vector = points[second] - points[first] + shift.cartesian(cell_matrix);
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
    box_(box) {
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

        // don't look for neighboring cells if we have only one cell and no
        // periodic boundary condition
        if (n_cells[spatial] == 1 && !box.periodic(spatial)) {
            n_search_[spatial] = 0;
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
        if (box_.periodic(spatial)) {
            auto result = divmod(cell_index[spatial], cells_shape_[spatial]);
            shift[spatial] = std::get<0>(result);
            cell_index[spatial] = std::get<1>(result);
        } else {
            shift[spatial] = 0;
            cell_index[spatial] = std::clamp(cell_index[spatial], 0, static_cast<int32_t>(cells_shape_[spatial] - 1));
        }
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
static scalar_t (*alloc(scalar_t (*ptr)[N], size_t size, size_t new_size))[N] {
    auto* new_ptr = reinterpret_cast<scalar_t(*)[N]>(std::realloc(ptr, new_size * sizeof(scalar_t[N])));

    if (new_ptr == nullptr) {
        throw std::bad_alloc();
    }

    // initialize with a bit pattern that maps to NaN for double
    std::memset(new_ptr + size, 0b11111111, (new_size - size) * sizeof(scalar_t[N]));

    return new_ptr;
}

template <typename scalar_t>
static scalar_t* alloc(scalar_t* ptr, size_t size, size_t new_size) {
    auto* new_ptr = reinterpret_cast<scalar_t*>(std::realloc(ptr, new_size * sizeof(scalar_t)));

    if (new_ptr == nullptr) {
        throw std::bad_alloc();
    }

    // initialize with a bit pattern that maps to NaN for double
    std::memset(new_ptr + size, 0b11111111, (new_size - size) * sizeof(scalar_t));

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

    this->neighbors.pairs = new_pairs;
    this->neighbors.shifts = new_shifts;
    this->neighbors.distances = new_distances;
    this->neighbors.vectors = new_vectors;

    this->capacity = new_size;
}

void GrowableNeighborList::reset() {
    // set all allocated data to zero
    auto size = this->neighbors.length;
    std::memset(this->neighbors.pairs, 0, size * sizeof(size_t[2]));

    if (this->neighbors.shifts != nullptr) {
        std::memset(this->neighbors.shifts, 0, size * sizeof(int32_t[3]));
    }

    if (this->neighbors.distances != nullptr) {
        std::memset(this->neighbors.distances, 0, size * sizeof(double));
    }

    if (this->neighbors.vectors != nullptr) {
        std::memset(this->neighbors.vectors, 0, size * sizeof(double[3]));
    }

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
            if (pairs[a][0] == pairs[b][0]) {
                return pairs[a][1] < pairs[b][1];
            } else {
                return pairs[a][0] < pairs[b][0];
            }
        }

        size_t (*pairs)[2];
    };

    std::sort(std::begin(indices), std::end(indices), compare_pairs(this->neighbors.pairs));

    // step 2: permute all data according to the sorted indices.
    int64_t cur = 0;
    int64_t is_sorted_up_to = 0;
    // data in `from` should go to `cur`
    auto from = indices[cur];

    size_t tmp_pair[2] = {0};
    double tmp_distance = 0;
    double tmp_vector[3] = {0};
    int32_t tmp_shift[3] = {0};

    while (cur < this->length()) {
        // move data from `cur` to temporary
        std::swap(tmp_pair, this->neighbors.pairs[cur]);
        if (options.return_distances) {
            std::swap(tmp_distance, this->neighbors.distances[cur]);
        }
        if (options.return_vectors) {
            std::swap(tmp_vector, this->neighbors.vectors[cur]);
        }
        if (options.return_shifts) {
            std::swap(tmp_shift, this->neighbors.shifts[cur]);
        }

        from = indices[cur];
        do {
            if (from == cur) {
                // permutation loop of a single entry, i.e. this value stayed
                // where is already was
                break;
            }
            // move data from `from` to `cur`
            std::swap(this->neighbors.pairs[cur], this->neighbors.pairs[from]);
            if (options.return_distances) {
                std::swap(this->neighbors.distances[cur], this->neighbors.distances[from]);
            }
            if (options.return_vectors) {
                std::swap(this->neighbors.vectors[cur], this->neighbors.vectors[from]);
            }
            if (options.return_shifts) {
                std::swap(this->neighbors.shifts[cur], this->neighbors.shifts[from]);
            }

            // mark this spot as already visited
            indices[cur] = -1;

            // update the indices
            cur = from;
            from = indices[cur];
        } while (indices[from] != -1);

        // we found a full loop of permutation, we can put tmp into `cur`
        std::swap(this->neighbors.pairs[cur], tmp_pair);
        if (options.return_distances) {
            std::swap(this->neighbors.distances[cur], tmp_distance);
        }
        if (options.return_vectors) {
            std::swap(this->neighbors.vectors[cur], tmp_vector);
        }
        if (options.return_shifts) {
            std::swap(this->neighbors.shifts[cur], tmp_shift);
        }

        indices[cur] = -1;

        // look for the next loop of permutation
        cur = is_sorted_up_to;
        while (indices[cur] == -1) {
            cur += 1;
            is_sorted_up_to += 1;
            if (cur == this->length()) {
                break;
            }
        }
    }
}

void PLMD::metatomic::vesin::cpu::free_neighbors(VesinNeighborList& neighbors) {
    assert(neighbors.device.type == VesinCPU);

    std::free(neighbors.pairs);
    std::free(neighbors.shifts);
    std::free(neighbors.vectors);
    std::free(neighbors.distances);
}
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <optional>
#include <stdexcept>

#define NOMINMAX
// gpu-lite - Combined Header
// A lightweight C++ library for dynamic CUDA runtime compilation and kernel caching
#ifndef GPULITE_HPP
#define GPULITE_HPP

#if defined(_MSC_VER)
  // MSVC historically reports __cplusplus wrong unless /Zc:__cplusplus is enabled,
  // so prefer _MSVC_LANG there.
  #if !defined(_MSVC_LANG) || _MSVC_LANG < 201703L
    #error "This project requires C++17 or newer (/std:c++17)."
  #endif
#else
  #if __cplusplus < 201703L
    #error "This project requires C++17 or newer (-std=c++17)."
  #endif
#endif

// =============================================================================
// CUDA Types Wrapper - Minimal CUDA type definitions for build-time independence
// =============================================================================

#include <cstddef>

#ifdef __cplusplus
extern "C" {
#endif

// CUDA Driver API types (from cuda.h)
typedef int CUresult;
typedef int CUdevice;
typedef struct CUctx_st* CUcontext;
typedef struct CUmod_st* CUmodule;
typedef struct CUfunc_st* CUfunction;
typedef struct CUstream_st* CUstream;
typedef unsigned long long CUdeviceptr;

// CUDA Runtime API types (from cuda_runtime.h)
typedef int cudaError_t;
typedef void* cudaStream_t;

// NVRTC types (from nvrtc.h)
typedef int nvrtcResult;
typedef struct _nvrtcProgram* nvrtcProgram;

// dim3 structure for kernel launch parameters
struct dim3 {
    unsigned int x, y, z;

    dim3(unsigned int x = 1, unsigned int y = 1, unsigned int z = 1) : x(x), y(y), z(z) {}
};

// CUDA memory copy kinds
typedef enum cudaMemcpyKind {
    cudaMemcpyHostToHost = 0,
    cudaMemcpyHostToDevice = 1,
    cudaMemcpyDeviceToHost = 2,
    cudaMemcpyDeviceToDevice = 3,
    cudaMemcpyDefault = 4
} cudaMemcpyKind;

//global definition for CUDA memory types
typedef enum cudaMemoryType {
      cudaMemoryTypeUnregistered = 0,
      cudaMemoryTypeHost = 1,
      cudaMemoryTypeDevice = 2,
      cudaMemoryTypeManaged = 3
} cudaMemoryType;


// CUDA pointer attributes structure
typedef struct cudaPointerAttributes {
    enum cudaMemoryType type;
    int device;
    void* devicePointer;
    void* hostPointer;
} cudaPointerAttributes;

// CUDA Driver API constants
enum {
    CUDA_SUCCESS = 0,
    CUDA_ERROR_NOT_INITIALIZED = 3
};

// CUDA Runtime API constants
enum {
    cudaSuccess = 0
};

// CUDA Host Alloc flags
enum {
    cudaHostAllocDefault = 0x00,
    cudaHostAllocPortable = 0x01,
    cudaHostAllocMapped = 0x02,
    cudaHostAllocWriteCombined = 0x04
};

// NVRTC constants
enum {
    NVRTC_SUCCESS = 0
};

// CUDA device attributes
typedef enum CUdevice_attribute {
    CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK = 8,
    CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK_OPTIN = 97,
    CU_DEVICE_ATTRIBUTE_RESERVED_SHARED_MEMORY_PER_BLOCK = 83,
    CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR = 75,
    CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR = 76
} CUdevice_attribute;

// CUDA function attributes
typedef enum CUfunction_attribute {
    CU_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES = 8,
    CU_FUNC_ATTRIBUTE_NUM_REGS = 4
} CUfunction_attribute;

// CUDA pointer attributes for cuPointerGetAttribute
typedef enum CUpointer_attribute {
    CU_POINTER_ATTRIBUTE_CONTEXT = 1,
    CU_POINTER_ATTRIBUTE_MEMORY_TYPE = 2
} CUpointer_attribute;

// CUDA memory types
enum {
    CU_MEMORYTYPE_HOST = 0x01,
    CU_MEMORYTYPE_DEVICE = 0x02,
    CU_MEMORYTYPE_ARRAY = 0x03,
    CU_MEMORYTYPE_UNIFIED = 0x04
};

typedef enum CUjit_option_enum {
    CU_JIT_GENERATE_DEBUG_INFO = 11,
} CUjit_option;

#ifdef __cplusplus
}
#endif

// =============================================================================
// Dynamic CUDA - Dynamic loading of CUDA runtime libraries
// =============================================================================

#if defined(__linux__) || defined(__APPLE__)
#include <dlfcn.h>
#include <unistd.h>  // for getcwd
#elif defined(_WIN32)
#include <windows.h>
#include <libloaderapi.h>

#include <direct.h>  // for _getcwd
#define getcwd _getcwd

#include <filesystem>
#else
#error "Platform not supported"
#endif

#include <stdexcept>
#include <mutex>
#include <string>
#include <unordered_map>
#include <sstream>
#include <list>
#include <vector>

#define NVRTC_SAFE_CALL(x)                                                                         \
    do {                                                                                           \
        nvrtcResult result = x;                                                                    \
        if (result != NVRTC_SUCCESS) {                                                             \
            std::ostringstream errorMsg;                                                           \
            errorMsg << "\nerror: " #x " failed with error "                                       \
                     << NVRTC_INSTANCE.nvrtcGetErrorString(result) << '\n'                         \
                     << "File: " << __FILE__ << '\n'                                               \
                     << "Line: " << static_cast<int>(__LINE__) << '\n';                            \
            throw std::runtime_error(errorMsg.str());                                              \
        }                                                                                          \
    } while (0)

#define CUDADRIVER_SAFE_CALL(x)                                                                    \
    do {                                                                                           \
        CUresult result = x;                                                                       \
        if (result != CUDA_SUCCESS) {                                                              \
            const char* msg;                                                                       \
            CUDA_DRIVER_INSTANCE.cuGetErrorName(result, &msg);                                     \
            std::ostringstream errorMsg;                                                           \
            errorMsg << "\nerror: " #x " failed with error " << (msg ? msg : "Unknown error")      \
                     << '\n'                                                                       \
                     << "File: " << __FILE__ << '\n'                                               \
                     << "Line: " << static_cast<int>(__LINE__) << '\n';                            \
            throw std::runtime_error(errorMsg.str());                                              \
        }                                                                                          \
    } while (0)

#define CUDART_SAFE_CALL(call)                                                                     \
    do {                                                                                           \
        cudaError_t cudaStatus = (call);                                                           \
        if (cudaStatus != cudaSuccess) {                                                           \
            std::ostringstream errorMsg;                                                           \
            const char* error = CUDART_INSTANCE.cudaGetErrorString(cudaStatus);                    \
            errorMsg << "\nfailed with error " << (error ? error : "Unknown error") << '\n'        \
                     << "File: " << __FILE__ << '\n'                                               \
                     << "Line: " << static_cast<int>(__LINE__) << '\n';                            \
            throw std::runtime_error(errorMsg.str());                                              \
        }                                                                                          \
    } while (0)

// Define a template to dynamically load symbols
template <typename FuncType> FuncType load(void* handle, const char* functionName) {
#if defined(__linux__) || defined(__APPLE__)
    auto func = reinterpret_cast<FuncType>(dlsym(handle, functionName));
#elif defined(_WIN32)
    auto func = reinterpret_cast<FuncType>(GetProcAddress(static_cast<HMODULE>(handle), functionName));
#endif
    if (!func) {
        throw std::runtime_error(std::string("Failed to load function: ") + functionName);
    }
    return func;
}

#ifdef _WIN32

namespace details {

inline std::wstring GetEnvVar(const wchar_t* name) {
    DWORD n = GetEnvironmentVariableW(name, nullptr, 0);
    if (n == 0) return L"";
    std::wstring val(n, L'\0');
    GetEnvironmentVariableW(name, val.data(), n);
    if (!val.empty() && val.back() == L'\0') val.pop_back();
    return val;
}

// Parse versions from filenames like cudart64_90.dll, cudart64_12.dll, cudart64_120.dll
inline int ParseCudartVersionScore(std::wstring prefix, const std::wstring& filename) {
    // Return a score; higher = preferred. Unknown parse => 0.

    prefix = prefix + L"_";
    // We try to parse digits after `prefix` and before ".dll".
    if (filename.rfind(prefix, 0) != 0) {
        return 0;
    }

    size_t start = prefix.size();
    size_t end = filename.find(L".dll");
    if (end == std::wstring::npos || end <= start) {
        return 0;
    }

    std::wstring num = filename.substr(start, end - start);
    if (num.empty()) {
        return 0;
    }
    for (wchar_t c : num) {
        if (c < L'0' || c > L'9') {
            return 0;
        }
    }

    // e.g. "90" -> 90, "12" -> 12, "120" -> 120
    return std::stoi(num);
}


inline std::vector<std::filesystem::path> CandidateCudaDirs() {
    std::vector<std::filesystem::path> dirs;

    // 1) CUDA_PATH\bin and CUDA_PATH\bin\x64
    std::wstring cudaPath = GetEnvVar(L"CUDA_PATH");
    if (!cudaPath.empty()) {
        dirs.push_back(std::filesystem::path(cudaPath) / L"bin");
        dirs.push_back(std::filesystem::path(cudaPath) / L"bin" / L"x64");
    }

    // 2) Search in PATH
    std::wstring pathEnv = GetEnvVar(L"PATH");
    if (!pathEnv.empty()) {
        size_t start = 0;
        size_t end = pathEnv.find(L';');
        while (end != std::wstring::npos) {
            std::wstring token = pathEnv.substr(start, end - start);
            if (!token.empty()) {
                dirs.push_back(std::filesystem::path(token));
            }
            start = end + 1;
            end = pathEnv.find(L';', start);
        }
        std::wstring token = pathEnv.substr(start);
        if (!token.empty()) {
            dirs.push_back(std::filesystem::path(token));
        }
    }

    // 3) Default toolkit install root (scan versions)
    //    C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v*\bin
    std::filesystem::path root = L"C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA";

    std::error_code ec;
    std::filesystem::directory_options options = std::filesystem::directory_options::skip_permission_denied;
    for (auto& entry: std::filesystem::directory_iterator(root, options, ec)) {
        if (ec) {
            break;
        }

        // folders like v12.4, v13.0, etc.
        if (!entry.is_directory(ec) || ec) {
            continue;
        }

        dirs.push_back(entry.path() / L"bin");
        dirs.push_back(entry.path() / L"bin" / L"x64");
    }

    // De-dup + keep only existing dirs
    std::sort(dirs.begin(), dirs.end());
    dirs.erase(std::unique(dirs.begin(), dirs.end()), dirs.end());
    dirs.erase(
        std::remove_if(dirs.begin(), dirs.end(), [](const std::filesystem::path& p) {
            assert(!p.empty());

            std::error_code ec;
            if (!std::filesystem::exists(p, ec) || ec) {
                return true;
            }

            if (!std::filesystem::is_directory(p, ec) || ec) {
                return true;
            }

            return false;

        }),
        dirs.end()
    );
    return dirs;
}

inline std::optional<std::filesystem::path> FindBestCudaDll(const std::wstring& prefix) {
    auto dirs = CandidateCudaDirs();

    struct Match {
        std::filesystem::path path;
        int score;
    };
    std::vector<Match> matches;

    for (const auto& d: dirs) {
        std::error_code ec;
        std::filesystem::directory_options options = std::filesystem::directory_options::skip_permission_denied;
        for (auto& e: std::filesystem::directory_iterator(d, options, ec)) {
            if (ec) {
                break;
            }

            if (!e.is_regular_file(ec) || ec) {
                continue;
            }

            auto name = e.path().filename().wstring();
            // Must look like <prefix>*.dll
            if (name.size() < prefix.size() + 4) {
                continue;
            }

            if (name.rfind(prefix, 0) != 0) {
                continue;
            }

            if (e.path().extension().wstring() != L".dll") {
                continue;
            }

            int score = ParseCudartVersionScore(prefix, name);
            // Prefer versioned DLLs; still accept plain "<prefix>.dll" with score 1
            if (name == prefix + L".dll") {
                score = std::max(score, 1);
            }

            matches.push_back({ e.path(), score });
        }
    }

    if (matches.empty()) {
        return std::nullopt;
    }

    // Prefer highest score; if tie, prefer shortest path (arbitrary stable tie-break)
    std::sort(matches.begin(), matches.end(), [](const Match& a, const Match& b) {
        if (a.score != b.score) return a.score > b.score;
        return a.path.wstring().size() < b.path.wstring().size();
    });

    return matches.front().path;
}

} // namespace details

#endif

/*
This class allows us to dynamically load the CUDA runtime and reference the functions contained
within the libcudart.so library (see CUDA Runtime API:
https://docs.nvidia.com/cuda/cuda-runtime-api/index.html).
*/
class CUDART {
  public:
    static CUDART& instance() {
        static CUDART instance;
        return instance;
    }

    bool loaded() { return cudartHandle != nullptr; }

    using cudaGetDeviceCount_t = cudaError_t (*)(int*);
    using cudaGetDevice_t = cudaError_t (*)(int*);
    using cudaSetDevice_t = cudaError_t (*)(int);
    using cudaMalloc_t = cudaError_t (*)(void**, size_t);
    using cudaMemcpy_t = cudaError_t (*)(void*, const void*, size_t, cudaMemcpyKind);
    using cudaMemset_t = cudaError_t (*)(void*, int, size_t);
    using cudaGetErrorName_t = const char* (*)(cudaError_t);
    using cudaGetErrorString_t = const char* (*)(cudaError_t);
    using cudaDeviceSynchronize_t = cudaError_t (*)(void);
    using cudaPointerGetAttributes_t = cudaError_t (*)(cudaPointerAttributes*, const void*);
    using cudaFree_t = cudaError_t (*)(void*);
    using cudaRuntimeGetVersion_t = cudaError_t (*)(int*);
    using cudaStreamCreate_t = cudaError_t (*)(cudaStream_t*);
    using cudaStreamDestroy_t = cudaError_t (*)(cudaStream_t);
    using cudaStreamSynchronize_t = cudaError_t (*)(cudaStream_t);
    using cudaHostAlloc_t = cudaError_t (*)(void**, size_t, unsigned int);
    using cudaFreeHost_t = cudaError_t (*)(void*);
    using cudaHostGetDevicePointer_t = cudaError_t (*)(void**, void*, unsigned int);
    using cudaMemcpyAsync_t = cudaError_t (*)(void*, const void*, size_t, cudaMemcpyKind, cudaStream_t);

    cudaGetDeviceCount_t cudaGetDeviceCount;
    cudaGetDevice_t cudaGetDevice;
    cudaSetDevice_t cudaSetDevice;
    cudaMalloc_t cudaMalloc;
    cudaMemset_t cudaMemset;
    cudaMemcpy_t cudaMemcpy;
    cudaGetErrorName_t cudaGetErrorName;
    cudaGetErrorString_t cudaGetErrorString;
    cudaDeviceSynchronize_t cudaDeviceSynchronize;
    cudaPointerGetAttributes_t cudaPointerGetAttributes;
    cudaFree_t cudaFree;
    cudaRuntimeGetVersion_t cudaRuntimeGetVersion;
    cudaStreamCreate_t cudaStreamCreate;
    cudaStreamDestroy_t cudaStreamDestroy;
    cudaStreamSynchronize_t cudaStreamSynchronize;
    cudaHostAlloc_t cudaHostAlloc;
    cudaFreeHost_t cudaFreeHost;
    cudaHostGetDevicePointer_t cudaHostGetDevicePointer;
    cudaMemcpyAsync_t cudaMemcpyAsync;

    CUDART() {
#if defined(__linux__) || defined(__APPLE__)
        static const char* CANDIDATES[] = {
            "libcudart.so",
            "libcudart.so.11",
            "libcudart.so.12",
            "libcudart.so.13",
            "libcudart.so.14",
            "libcudart.so.15",
        };
        for (auto* candidate: CANDIDATES) {
            cudartHandle = dlopen(candidate, RTLD_NOW);
            if (cudartHandle) {
                break;
            }
        }
#elif defined(_WIN32)
        auto dllPathOpt = details::FindBestCudaDll(L"cudart64");
        if (dllPathOpt) {
            auto dllPath = *dllPathOpt;
            auto dir = dllPath.parent_path();
            // add the directory containing the DLL to the search path
            SetDllDirectoryW(dir.c_str());

            cudartHandle = LoadLibraryExW(
                dllPath.c_str(),
                nullptr,
                LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR |
                LOAD_LIBRARY_SEARCH_DEFAULT_DIRS |
                LOAD_LIBRARY_SEARCH_USER_DIRS
            );
        }
#else
#error "Platform not supported"
#endif
        if (cudartHandle) {
            // load cudart function pointers using template
            cudaGetDeviceCount = load<cudaGetDeviceCount_t>(cudartHandle, "cudaGetDeviceCount");
            cudaGetDevice = load<cudaGetDevice_t>(cudartHandle, "cudaGetDevice");
            cudaSetDevice = load<cudaSetDevice_t>(cudartHandle, "cudaSetDevice");
            cudaMalloc = load<cudaMalloc_t>(cudartHandle, "cudaMalloc");
            cudaMemset = load<cudaMemset_t>(cudartHandle, "cudaMemset");
            cudaMemcpy = load<cudaMemcpy_t>(cudartHandle, "cudaMemcpy");
            cudaGetErrorName = load<cudaGetErrorName_t>(cudartHandle, "cudaGetErrorName");
            cudaGetErrorString = load<cudaGetErrorString_t>(cudartHandle, "cudaGetErrorString");
            cudaDeviceSynchronize =
                load<cudaDeviceSynchronize_t>(cudartHandle, "cudaDeviceSynchronize");
            cudaPointerGetAttributes =
                load<cudaPointerGetAttributes_t>(cudartHandle, "cudaPointerGetAttributes");
            cudaFree = load<cudaFree_t>(cudartHandle, "cudaFree");
            cudaRuntimeGetVersion =
                load<cudaRuntimeGetVersion_t>(cudartHandle, "cudaRuntimeGetVersion");
            cudaStreamCreate = load<cudaStreamCreate_t>(cudartHandle, "cudaStreamCreate");
            cudaStreamDestroy = load<cudaStreamDestroy_t>(cudartHandle, "cudaStreamDestroy");
            cudaStreamSynchronize = load<cudaStreamSynchronize_t>(cudartHandle, "cudaStreamSynchronize");
            cudaHostAlloc = load<cudaHostAlloc_t>(cudartHandle, "cudaHostAlloc");
            cudaFreeHost = load<cudaFreeHost_t>(cudartHandle, "cudaFreeHost");
            cudaHostGetDevicePointer = load<cudaHostGetDevicePointer_t>(cudartHandle, "cudaHostGetDevicePointer");
            cudaMemcpyAsync = load<cudaMemcpyAsync_t>(cudartHandle, "cudaMemcpyAsync");
        }
    }

    ~CUDART() {
#if defined(__linux__) || defined(__APPLE__)
        if (cudartHandle) {
            dlclose(cudartHandle);
        }
#elif defined(_WIN32)
        if (cudartHandle) {
            FreeLibrary(static_cast<HMODULE>(cudartHandle));
        }
#else
#error "Platform not supported"
#endif
    }

    // Prevent copying
    CUDART(const CUDART&) = delete;
    CUDART& operator=(const CUDART&) = delete;

    void* cudartHandle = nullptr;
};

/*
This class allows us to dynamically load the CUDA Driver and reference the functions contained
within the libcuda.so library (CUDA Driver API:
https://docs.nvidia.com/cuda/cuda-driver-api/index.html).
*/
class CUDADriver {

  public:
    static CUDADriver& instance() {
        static CUDADriver instance;
        return instance;
    }

    bool loaded() { return cudaHandle != nullptr; }

    using cuInit_t = CUresult (*)(unsigned int);
    using cuDeviceGetCount_t = CUresult (*)(int*);
    using cuDevicePrimaryCtxRetain_t = CUresult (*)(CUcontext*, CUdevice);
    using cuDevicePrimaryCtxRelease_t = CUresult (*)(CUdevice);
    using cuCtxCreate_t = CUresult (*)(CUcontext*, unsigned int, CUdevice);
    using cuCtxDestroy_t = CUresult (*)(CUcontext);
    using cuCtxGetCurrent_t = CUresult (*)(CUcontext*);
    using cuCtxSetCurrent_t = CUresult (*)(CUcontext);
    using cuModuleLoadDataEx_t = CUresult (*)(CUmodule*, const void*, unsigned int, CUjit_option*, void**);
    using cuModuleGetFunction_t = CUresult (*)(CUfunction*, CUmodule, const char*);
    using cuFuncSetAttribute_t = CUresult (*)(CUfunction, CUfunction_attribute, int);
    using cuFuncGetAttribute_t = CUresult (*)(int*, CUfunction_attribute, CUfunction);
    using cuCtxGetDevice_t = CUresult (*)(CUdevice*);
    using cuDeviceGetAttribute_t = CUresult (*)(int*, CUdevice_attribute, CUdevice);
    using cuDeviceGetName_t = CUresult (*)(char*, int, CUdevice);
    using cuDeviceTotalMem_t = CUresult (*)(size_t*, CUdevice);
    using cuLaunchKernel_t = CUresult (*)(
        CUfunction,
        unsigned int,
        unsigned int,
        unsigned int,
        unsigned int,
        unsigned int,
        unsigned int,
        size_t,
        CUstream,
        void**,
        void*
    );
    using cuStreamCreate_t = CUresult (*)(CUstream*, unsigned int);
    using cuStreamDestroy_t = CUresult (*)(CUstream);
    using cuCtxSynchronize_t = CUresult (*)(void);
    using cuGetErrorName_t = CUresult (*)(CUresult, const char**);
    using cuCtxPushCurrent_t = CUresult (*)(CUcontext);
    using cuPointerGetAttribute_t = CUresult (*)(void*, CUpointer_attribute, CUdeviceptr);

    cuInit_t cuInit;
    cuDeviceGetCount_t cuDeviceGetCount;
    cuCtxCreate_t cuCtxCreate;
    cuCtxDestroy_t cuCtxDestroy;
    cuDevicePrimaryCtxRetain_t cuDevicePrimaryCtxRetain;
    cuDevicePrimaryCtxRelease_t cuDevicePrimaryCtxRelease;
    cuCtxGetCurrent_t cuCtxGetCurrent;
    cuCtxSetCurrent_t cuCtxSetCurrent;
    cuModuleLoadDataEx_t cuModuleLoadDataEx;
    cuModuleGetFunction_t cuModuleGetFunction;
    cuFuncSetAttribute_t cuFuncSetAttribute;
    cuFuncGetAttribute_t cuFuncGetAttribute;
    cuCtxGetDevice_t cuCtxGetDevice;
    cuDeviceGetAttribute_t cuDeviceGetAttribute;
    cuDeviceGetName_t cuDeviceGetName;
    cuDeviceTotalMem_t cuDeviceTotalMem;
    cuLaunchKernel_t cuLaunchKernel;
    cuStreamCreate_t cuStreamCreate;
    cuStreamDestroy_t cuStreamDestroy;
    cuGetErrorName_t cuGetErrorName;
    cuCtxSynchronize_t cuCtxSynchronize;
    cuCtxPushCurrent_t cuCtxPushCurrent;
    cuPointerGetAttribute_t cuPointerGetAttribute;

    CUDADriver() {
#if defined(__linux__) || defined(__APPLE__)
        cudaHandle = dlopen("libcuda.so", RTLD_NOW);
#elif defined(_WIN32)
        cudaHandle = LoadLibraryA("nvcuda.dll");
#else
#error "Platform not supported"
#endif
        if (cudaHandle) {
            // Load CUDA driver function pointers using template
            cuInit = load<cuInit_t>(cudaHandle, "cuInit");
            cuDeviceGetCount = load<cuDeviceGetCount_t>(cudaHandle, "cuDeviceGetCount");
            cuCtxCreate = load<cuCtxCreate_t>(cudaHandle, "cuCtxCreate");
            cuCtxDestroy = load<cuCtxDestroy_t>(cudaHandle, "cuCtxDestroy");
            cuDevicePrimaryCtxRetain =
                load<cuDevicePrimaryCtxRetain_t>(cudaHandle, "cuDevicePrimaryCtxRetain");
            cuDevicePrimaryCtxRelease =
                load<cuDevicePrimaryCtxRelease_t>(cudaHandle, "cuDevicePrimaryCtxRelease");
            cuCtxGetCurrent = load<cuCtxGetCurrent_t>(cudaHandle, "cuCtxGetCurrent");
            cuCtxSetCurrent = load<cuCtxSetCurrent_t>(cudaHandle, "cuCtxSetCurrent");
            cuModuleLoadDataEx = load<cuModuleLoadDataEx_t>(cudaHandle, "cuModuleLoadDataEx");
            cuModuleGetFunction = load<cuModuleGetFunction_t>(cudaHandle, "cuModuleGetFunction");
            cuFuncSetAttribute = load<cuFuncSetAttribute_t>(cudaHandle, "cuFuncSetAttribute");
            cuFuncGetAttribute = load<cuFuncGetAttribute_t>(cudaHandle, "cuFuncGetAttribute");
            cuCtxGetDevice = load<cuCtxGetDevice_t>(cudaHandle, "cuCtxGetDevice");
            cuDeviceGetAttribute = load<cuDeviceGetAttribute_t>(cudaHandle, "cuDeviceGetAttribute");
            cuDeviceGetName = load<cuDeviceGetName_t>(cudaHandle, "cuDeviceGetName");
            cuDeviceTotalMem = load<cuDeviceTotalMem_t>(cudaHandle, "cuDeviceTotalMem");
            cuLaunchKernel = load<cuLaunchKernel_t>(cudaHandle, "cuLaunchKernel");
            cuStreamCreate = load<cuStreamCreate_t>(cudaHandle, "cuStreamCreate");
            cuStreamDestroy = load<cuStreamDestroy_t>(cudaHandle, "cuStreamDestroy");
            cuCtxSynchronize = load<cuCtxSynchronize_t>(cudaHandle, "cuCtxSynchronize");
            cuGetErrorName = load<cuGetErrorName_t>(cudaHandle, "cuGetErrorName");
            cuCtxPushCurrent = load<cuCtxPushCurrent_t>(cudaHandle, "cuCtxPushCurrent");
            cuPointerGetAttribute =
                load<cuPointerGetAttribute_t>(cudaHandle, "cuPointerGetAttribute");
        }
    }

    ~CUDADriver() {
#if defined(__linux__) || defined(__APPLE__)
        if (cudaHandle) {
            dlclose(cudaHandle);
        }
#elif defined(_WIN32)
        if (cudaHandle) {
            FreeLibrary(static_cast<HMODULE>(cudaHandle));
        }
#else
#error "Platform not supported"
#endif
    }

    // Prevent copying
    CUDADriver(const CUDADriver&) = delete;
    CUDADriver& operator=(const CUDADriver&) = delete;

    void* cudaHandle = nullptr;
};

/*
This class allows us to dynamically load NVRTC and reference the functions contained within the
libnvrtc.so library (see NVRTC API: https://docs.nvidia.com/cuda/nvrtc/index.html).
*/
class NVRTC {

  public:
    static NVRTC& instance() {
        static NVRTC instance;
        return instance;
    }

    bool loaded() { return nvrtcHandle != nullptr; }

    using nvrtcCreateProgram_t =
        nvrtcResult (*)(nvrtcProgram*, const char*, const char*, int, const char*[], const char*[]);
    using nvrtcCompileProgram_t = nvrtcResult (*)(nvrtcProgram, int, const char*[]);
    using nvrtcGetPTX_t = nvrtcResult (*)(nvrtcProgram, char*);
    using nvrtcGetPTXSize_t = nvrtcResult (*)(nvrtcProgram, size_t*);
    using nvrtcGetCUBIN_t = nvrtcResult (*)(nvrtcProgram, char*);
    using nvrtcGetCUBINSize_t = nvrtcResult (*)(nvrtcProgram, size_t*);
    using nvrtcGetProgramLog_t = nvrtcResult (*)(nvrtcProgram, char*);
    using nvrtcGetProgramLogSize_t = nvrtcResult (*)(nvrtcProgram, size_t*);
    using nvrtcAddNameExpression_t = nvrtcResult (*)(nvrtcProgram, const char* const);
    using nvrtcGetLoweredName_t = nvrtcResult (*)(nvrtcProgram, const char*, const char**);
    using nvrtcDestroyProgram_t = nvrtcResult (*)(nvrtcProgram*);
    using nvrtcGetErrorString_t = const char* (*)(nvrtcResult);

    nvrtcCreateProgram_t nvrtcCreateProgram;
    nvrtcCompileProgram_t nvrtcCompileProgram;
    nvrtcGetPTX_t nvrtcGetPTX;
    nvrtcGetPTXSize_t nvrtcGetPTXSize;
    nvrtcGetCUBIN_t nvrtcGetCUBIN;
    nvrtcGetCUBINSize_t nvrtcGetCUBINSize;
    nvrtcGetProgramLog_t nvrtcGetProgramLog;
    nvrtcGetProgramLogSize_t nvrtcGetProgramLogSize;
    nvrtcGetLoweredName_t nvrtcGetLoweredName;
    nvrtcAddNameExpression_t nvrtcAddNameExpression;
    nvrtcDestroyProgram_t nvrtcDestroyProgram;
    nvrtcGetErrorString_t nvrtcGetErrorString;

    NVRTC() {
#if defined(__linux__) || defined(__APPLE__)
        static const char* CANDIDATES[] = {
            "libnvrtc.so",
            "libnvrtc.so.11",
            "libnvrtc.so.12",
            "libnvrtc.so.13",
            "libnvrtc.so.14",
            "libnvrtc.so.15",
        };
        for (auto* candidate: CANDIDATES) {
            nvrtcHandle = dlopen(candidate, RTLD_NOW);
            if (nvrtcHandle != nullptr) {
                break;
            }
        }

#elif defined(_WIN32)
        auto dllPathOpt = details::FindBestCudaDll(L"nvrtc64");
        if (dllPathOpt) {
            auto dllPath = *dllPathOpt;
            // add the directory containing the DLL to the search path
            auto dir = dllPath.parent_path();
            SetDllDirectoryW(dir.c_str());

            nvrtcHandle = LoadLibraryExW(
                dllPath.c_str(),
                nullptr,
                LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR |
                LOAD_LIBRARY_SEARCH_DEFAULT_DIRS |
                LOAD_LIBRARY_SEARCH_USER_DIRS
            );
        }
#else
#error "Platform not supported"
#endif

        if (nvrtcHandle) {
            // Load NVRTC function pointers using template
            nvrtcCreateProgram = load<nvrtcCreateProgram_t>(nvrtcHandle, "nvrtcCreateProgram");
            nvrtcCompileProgram = load<nvrtcCompileProgram_t>(nvrtcHandle, "nvrtcCompileProgram");
            nvrtcGetPTX = load<nvrtcGetPTX_t>(nvrtcHandle, "nvrtcGetPTX");
            nvrtcGetPTXSize = load<nvrtcGetPTXSize_t>(nvrtcHandle, "nvrtcGetPTXSize");
            nvrtcGetCUBIN = load<nvrtcGetCUBIN_t>(nvrtcHandle, "nvrtcGetCUBIN");
            nvrtcGetCUBINSize = load<nvrtcGetCUBINSize_t>(nvrtcHandle, "nvrtcGetCUBINSize");
            nvrtcGetProgramLog = load<nvrtcGetProgramLog_t>(nvrtcHandle, "nvrtcGetProgramLog");
            nvrtcGetProgramLogSize =
                load<nvrtcGetProgramLogSize_t>(nvrtcHandle, "nvrtcGetProgramLogSize");
            nvrtcGetLoweredName = load<nvrtcGetLoweredName_t>(nvrtcHandle, "nvrtcGetLoweredName");
            nvrtcAddNameExpression =
                load<nvrtcAddNameExpression_t>(nvrtcHandle, "nvrtcAddNameExpression");
            nvrtcDestroyProgram = load<nvrtcDestroyProgram_t>(nvrtcHandle, "nvrtcDestroyProgram");
            nvrtcGetErrorString = load<nvrtcGetErrorString_t>(nvrtcHandle, "nvrtcGetErrorString");
        }
    }

    ~NVRTC() {
#if defined(__linux__) || defined(__APPLE__)
        if (nvrtcHandle) {
            dlclose(nvrtcHandle);
        }
#elif defined(_WIN32)
        if (nvrtcHandle) {
            FreeLibrary(static_cast<HMODULE>(nvrtcHandle));
        }
#else
#error "Platform not supported"
#endif
    }

    // Prevent copying
    NVRTC(const NVRTC&) = delete;
    NVRTC& operator=(const NVRTC&) = delete;

    void* nvrtcHandle = nullptr;
};

#define CUDART_INSTANCE CUDART::instance()
#define CUDA_DRIVER_INSTANCE CUDADriver::instance()
#define NVRTC_INSTANCE NVRTC::instance()

// Convenience macros for cleaner API - hides the global instance from users
// Usage: GPULITE_CUDART_CALL(cudaMalloc(&ptr, size))
//        GPULITE_DRIVER_CALL(cuCtxGetCurrent(&ctx))
#define GPULITE_CUDART_CALL(func) CUDART_SAFE_CALL(CUDART_INSTANCE.func)
#define GPULITE_DRIVER_CALL(func) CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.func)

// =============================================================================
// CUDA Kernel Cache Manager - Runtime compilation and caching system
// =============================================================================

#include <fstream>
#include <vector>
#include <unordered_map>
#include <memory>
#include <typeinfo>

#if defined(__GNUC__) || defined(__clang__)
#include <cxxabi.h>
#endif

// Helper function to demangle the type name if necessary
inline std::string demangleTypeName(const std::string& name) {
#if defined(__GNUC__) || defined(__clang__)
    int status = 0;
    std::unique_ptr<char, void (*)(void*)> demangled_name(
        abi::__cxa_demangle(name.c_str(), 0, 0, &status), std::free
    );
    return (status == 0) ? demangled_name.get() : name;
#else
    // not ideal, but better than nothing on other compilers
    return name;
#endif
}

// Base case: No template arguments, return function name without any type information
inline std::string getKernelName(const std::string& fn_name) { return fn_name; }

// Function to get type name of a single type
template <typename T> std::string typeName() { return demangleTypeName(typeid(T).name()); }

// Variadic template function to build type list
template <typename T, typename... Ts> void buildTemplateTypes(std::string& base) {
    base += typeName<T>(); // Add the first type
    // If there are more types, add a comma and recursively call for the remaining types
    if constexpr (sizeof...(Ts) > 0) {
        base += ", ";
        buildTemplateTypes<Ts...>(base); // Recursively call for the rest of the types
    }
}

// Helper function to start building the types
template <typename T, typename... Ts> std::string buildTemplateTypes() {
    std::string result;
    buildTemplateTypes<T, Ts...>(result); // Use recursive variadic template
    return result;
}

/*
Function to get the kernel name with the list of templated types if any:
*/
template <typename T, typename... Ts> std::string getKernelName(const std::string& fn_name) {
    std::string type_list = buildTemplateTypes<T, Ts...>(); // Build type list
    return fn_name + "<" + type_list + ">"; // Return function name with type list in angle brackets
}

// Function to load CUDA source code from a file
inline std::string load_cuda_source(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    std::ostringstream ss;
    ss << file.rdbuf();
    return ss.str();
}

/*
Container class for the cached kernels. Provides functionality for launching compiled kernels as
well as automatically resizing dynamic shared memory allocations, when needed. Kernels are compiled
on first launch.
*/
class CachedKernel {

  public:
    CachedKernel(
        std::string kernel_name,
        std::string kernel_code,
        std::string source_name,
        std::vector<std::string> options
    ) {
        this->kernel_name = kernel_name;
        this->kernel_code = kernel_code;
        this->source_name = source_name;
        this->options = options;
    }

    CachedKernel() = default;

    // Copy constructor
    CachedKernel(const CachedKernel&) = default;

    // Copy assignment operator
    CachedKernel& operator=(const CachedKernel&) = default;

    inline void setFuncAttribute(CUfunction_attribute attribute, int value) const {
        CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuFuncSetAttribute(function, attribute, value));
    }

    int getFuncAttribute(CUfunction_attribute attribute) const {
        int value;
        CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuFuncGetAttribute(&value, attribute, function));
        return value;
    }

    /*
    launches the kernel, and optionally synchronizes until control can be passed back to host.
    */
    void launch(
        dim3 grid,
        dim3 block,
        size_t shared_mem_size,
        void* cuda_stream,
        std::vector<void*> args,
        bool synchronize = true
    ) {

        if (!compiled) {
            this->compileKernel(args);
        }

        CUcontext currentContext = nullptr;
        // Get current context
        CUresult result = CUDA_DRIVER_INSTANCE.cuCtxGetCurrent(&currentContext);

        if (result != CUDA_SUCCESS || !currentContext) {
            throw std::runtime_error("CachedKernel::launch error getting current context.");
        }

        if (currentContext != context) {
            CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuCtxSetCurrent(context));
        }

        this->checkAndAdjustSharedMem(shared_mem_size);

        CUstream cstream = reinterpret_cast<CUstream>(cuda_stream);

        CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuLaunchKernel(
            function,
            grid.x,
            grid.y,
            grid.z,
            block.x,
            block.y,
            block.z,
            shared_mem_size,
            cstream,
            args.data(),
            0
        ));

        if (synchronize) {
            CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuCtxSynchronize());
        }

        if (currentContext != context) {
            CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuCtxSetCurrent(currentContext));
        }
    }

  private:
    /*
    The default shared memory space on most recent NVIDIA cards is defaulted
    49152 bytes. This method
    attempts to adjust the shared memory to fit the requested configuration if
    the kernel launch parameters exceeds the default 49152 bytes.
    */
    void checkAndAdjustSharedMem(int query_shared_mem_size) {
        if (current_smem_size == 0) {
            CUdevice cuDevice;
            CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuCtxGetDevice(&cuDevice));

            CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuDeviceGetAttribute(
                &max_smem_size_optin, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK_OPTIN, cuDevice
            ));

            int reserved_smem_per_block = 0;

            CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuDeviceGetAttribute(
                &reserved_smem_per_block, CU_DEVICE_ATTRIBUTE_RESERVED_SHARED_MEMORY_PER_BLOCK, cuDevice
            ));

            int curr_max_smem_per_block = 0;

            CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuDeviceGetAttribute(
                &curr_max_smem_per_block, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, cuDevice
            ));

            current_smem_size = (curr_max_smem_per_block - reserved_smem_per_block);
        }

        if (query_shared_mem_size > current_smem_size) {

            if (query_shared_mem_size > max_smem_size_optin) {
                throw std::runtime_error(
                    "CachedKernel::launch requested more smem than available on card."
                );
            } else {
                CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuFuncSetAttribute(
                    function, CU_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES, query_shared_mem_size
                ));
                current_smem_size = query_shared_mem_size;
            }
        }
    }

    /*
        Compiles the kernel "kernel_name" located in source file "kernel_code", which additional
        parameters "options" passed to NVRTC_INSTANCE. Will auto-detect the compute capability of
       the available card. args for the launch need to be queried as we need to grab the CUcontext
       in which these ptrs exist.
        */
    void compileKernel(std::vector<void*>& kernel_args) {

        this->initCudaDriver();

        CUcontext currentContext = nullptr;

        for (size_t ptr_id = 0; ptr_id < kernel_args.size(); ptr_id++) {
            unsigned int memtype = 0;
            CUdeviceptr device_ptr = *reinterpret_cast<CUdeviceptr*>(kernel_args[ptr_id]);

            CUresult res = CUDA_DRIVER_INSTANCE.cuPointerGetAttribute(
                &memtype, CU_POINTER_ATTRIBUTE_MEMORY_TYPE, device_ptr
            );

            if (res == CUDA_SUCCESS && memtype == CU_MEMORYTYPE_DEVICE) {
                CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuPointerGetAttribute(
                    &currentContext, CU_POINTER_ATTRIBUTE_CONTEXT, device_ptr
                ));

                if (currentContext) {
                    break;
                }
            }
        }

        CUcontext query = nullptr;
        CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuCtxGetCurrent(&query));

        if (query != currentContext) {
            CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuCtxSetCurrent(currentContext));
        }

        CUdevice cuDevice;
        CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuCtxGetDevice(&cuDevice));

        // Check if debug option is enabled
        const bool enableDebug = std::any_of(
            this->options.cbegin(), this->options.cend(),
            [](const std::string& opt) {
                return opt == "-G" || opt == "--device-debug";
            }
        );

        // When debugging, write source to a real file so cuda-gdb can find it
        std::string effective_source_name = this->source_name;
        if (enableDebug) {
            // Create a debug source file in the current working directory
            // Use absolute path so cuda-gdb can reliably find it
            char cwd[4096];
            if (getcwd(cwd, sizeof(cwd)) != nullptr) {
                effective_source_name = std::string(cwd) + "/" + this->source_name;
            }

            std::ofstream debug_source_file(effective_source_name);
            if (debug_source_file.is_open()) {
                debug_source_file << this->kernel_code;
                debug_source_file.close();
            } else {
                throw std::runtime_error(
                    "Failed to write debug source file: " + effective_source_name
                );
            }
        }

        nvrtcProgram prog;

        NVRTC_SAFE_CALL(NVRTC_INSTANCE.nvrtcCreateProgram(
            &prog, this->kernel_code.c_str(), effective_source_name.c_str(), 0, nullptr, nullptr
        ));

        NVRTC_SAFE_CALL(NVRTC_INSTANCE.nvrtcAddNameExpression(prog, this->kernel_name.c_str()));

        std::vector<const char*> c_options;
        c_options.reserve(this->options.size());
        for (const auto& option : this->options) {
            c_options.push_back(option.c_str());
        }

        int major = 0;
        int minor = 0;
        CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuDeviceGetAttribute(
            &major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, cuDevice
        ));
        CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuDeviceGetAttribute(
            &minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, cuDevice
        ));
        int arch = major * 10 + minor;
        std::string smbuf = "--gpu-architecture=sm_" + std::to_string(arch);

        c_options.push_back(smbuf.c_str());

        nvrtcResult compileResult =
            NVRTC_INSTANCE.nvrtcCompileProgram(prog, c_options.size(), c_options.data());
        if (compileResult != NVRTC_SUCCESS) {
            size_t logSize;
            NVRTC_SAFE_CALL(NVRTC_INSTANCE.nvrtcGetProgramLogSize(prog, &logSize));
            std::string log(logSize, '\0');
            NVRTC_SAFE_CALL(NVRTC_INSTANCE.nvrtcGetProgramLog(prog, &log[0]));
            throw std::runtime_error(
                "KernelFactory::compileAndCacheKernel: Failed to compile CUDA program:\n" + log
            );
        }

        // fetch CUBIN
        size_t cubinSize = 0;
        NVRTC_SAFE_CALL(NVRTC_INSTANCE.nvrtcGetCUBINSize(prog, &cubinSize));
        std::vector<char> cubin(cubinSize);
        NVRTC_SAFE_CALL(NVRTC_INSTANCE.nvrtcGetCUBIN(prog, cubin.data()));

        // load the module from cubin
        CUmodule module = nullptr;
        CUresult cuResult;

        if (enableDebug) {
            // Load with JIT debug info
            CUjit_option opts[1];
            opts[0] = CU_JIT_GENERATE_DEBUG_INFO;
            void** vals = new void*[1];
            vals[0] = (void*)(size_t)1;
            cuResult = CUDA_DRIVER_INSTANCE.cuModuleLoadDataEx(
                &module, cubin.data(), 1, opts, vals
            );
            delete[] vals;
        } else {
            // Load without JIT options
            cuResult = CUDA_DRIVER_INSTANCE.cuModuleLoadDataEx(
                &module, cubin.data(), 0, 0, 0
            );
        }

        if (cuResult != CUDA_SUCCESS) {
            throw std::runtime_error(
                "KernelFactory::compileAndCacheKernel: Failed to load PTX code into CUDA "
                "module "
                "(error code: " +
                std::to_string(cuResult) + ")"
            );
        }

        const char* lowered_name;
        NVRTC_SAFE_CALL(
            NVRTC_INSTANCE.nvrtcGetLoweredName(prog, this->kernel_name.c_str(), &lowered_name)
        );
        CUfunction kernel;
        CUDADRIVER_SAFE_CALL(CUDA_DRIVER_INSTANCE.cuModuleGetFunction(&kernel, module, lowered_name));

        this->module = module;
        this->function = kernel;
        this->context = currentContext;
        this->compiled = true;

        NVRTC_SAFE_CALL(NVRTC_INSTANCE.nvrtcDestroyProgram(&prog));
    }

    void initCudaDriver() {

        int deviceCount = 0;
        // Check if CUDA has already been initialized
        CUresult res = CUDA_DRIVER_INSTANCE.cuDeviceGetCount(&deviceCount);
        if (res == CUDA_ERROR_NOT_INITIALIZED) {
            // CUDA hasn't been initialized, so we initialize it now
            res = CUDA_DRIVER_INSTANCE.cuInit(0);
            if (res != CUDA_SUCCESS) {
                throw std::runtime_error(
                    "KernelFactory::initCudaDriver: Failed to initialize CUDA CUDA_DRIVER_INSTANCE."
                );
                return;
            }
        }
    }

    int current_smem_size = 0;
    int max_smem_size_optin = 0;
    CUmodule module = nullptr;
    CUfunction function = nullptr;
    CUcontext context = nullptr;
    bool compiled = false;

    std::string kernel_name;
    std::string kernel_code;
    std::string source_name;
    std::vector<std::string> options;
};

/*
Factory class to create and store compiled cuda kernels for caching as a simple name-based hashmap.
Allows both compiling from a source file, or for compiling from a variable containing CUDA code.
*/
class KernelFactory {
  public:
    KernelFactory(const KernelFactory&) = delete;
    KernelFactory& operator=(const KernelFactory&) = delete;

    KernelFactory(KernelFactory&&) = default;
    KernelFactory& operator=(KernelFactory&&) = default;

    // Get the singleton instance of the KernelFactory for a given CUDA device.
    // This ensures that each CUDA device has its own kernel cache.
    static KernelFactory& instance(CUdevice device) {
        static std::list<KernelFactory> INSTANCES;
        for (size_t i = INSTANCES.size(); i < device + 1; i++) {
            INSTANCES.emplace_back(KernelFactory());
        }

        // get the element at index "device" in the list and return it
        auto it = INSTANCES.begin();
        std::advance(it, device);
        return *it;
    }

    void cacheKernel(
        const std::string& kernel_name,
        const std::string& source_path,
        const std::string& source_name,
        const std::vector<std::string>& options
    ) {
        std::lock_guard<std::mutex> kernel_cache_lock(kernel_cache_mutex_);
        kernel_cache_[kernel_name] =
            std::make_unique<CachedKernel>(kernel_name, source_path, source_name, options);
    }

    bool hasKernel(const std::string& kernel_name) {
        std::lock_guard<std::mutex> kernel_cache_lock(kernel_cache_mutex_);
        return kernel_cache_.find(kernel_name) != kernel_cache_.end();
    }

    CachedKernel* getKernel(const std::string& kernel_name) {
        std::lock_guard<std::mutex> kernel_cache_lock(kernel_cache_mutex_);
        auto it = kernel_cache_.find(kernel_name);
        if (it != kernel_cache_.end()) {
            return it->second.get();
        } else {
            throw std::runtime_error("Kernel not found in cache.");
        }
    }

    /*
    Tries to retrieve the kernel "kernel_name". If not found, instantiate it and save to cache.
    */
    CachedKernel* createFromSource(
        const std::string& kernel_name,
        const std::string& source_path,
        const std::string& source_name,
        const std::vector<std::string>& options
    ) {
        if (!this->hasKernel(kernel_name)) {
            std::string kernel_code = load_cuda_source(source_path);
            this->cacheKernel(kernel_name, kernel_code, source_name, options);
        }
        return this->getKernel(kernel_name);
    }

    /*
    Tries to retrieve the kernel "kernel_name". If not found, instantiate it and save to cache.
    */
    CachedKernel* create(
        const std::string& kernel_name,
        const std::string& source_variable,
        const std::string& source_name,
        const std::vector<std::string>& options
    ) {
        if (!this->hasKernel(kernel_name)) {
            this->cacheKernel(kernel_name, source_variable, source_name, options);
        }

        return this->getKernel(kernel_name);
    }

  private:
    KernelFactory() {}
    std::unordered_map<std::string, std::unique_ptr<CachedKernel>> kernel_cache_;

    static std::mutex kernel_cache_mutex_;
};

inline std::mutex KernelFactory::kernel_cache_mutex_;

#endif // GPULITE_HPP

#ifndef VESIN_CUDA_HPP
#define VESIN_CUDA_HPP


namespace PLMD {
namespace metatomic {
namespace vesin {
namespace cuda {

#ifndef VESIN_DEFAULT_CUDA_MAX_PAIRS_PER_POINT
/// @brief Default maximum number of pairs per point on the GPU (can be
/// overridden).
#define VESIN_DEFAULT_CUDA_MAX_PAIRS_PER_POINT 256
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

using namespace PLMD::metatomic::vesin::cuda;

// NVTX for profiling (optional, enabled if available)
#ifdef VESIN_ENABLE_NVTX
#include <nvtx3/nvToolsExt.h>
#define NVTX_PUSH(name) nvtxRangePushA(name)
#define NVTX_POP() nvtxRangePop()
#else
#define NVTX_PUSH(name) \
    do {                \
    } while (0)
#define NVTX_POP() \
    do {           \
    } while (0)
#endif

static const char* CUDA_BRUTEFORCE_CODE =
R"======(#define NWARPS 4
)======"
R"======(#define WARP_SIZE 32
)======"
R"======(
)======"
R"======(__device__ inline size_t atomicAdd_size_t(size_t* address, size_t val) {
)======"
R"======(    return static_cast<size_t>(atomicAdd(
)======"
R"======(        reinterpret_cast<unsigned long long*>(address),
)======"
R"======(        static_cast<unsigned long long>(val)
)======"
R"======(    ));
)======"
R"======(}
)======"
R"======(
)======"
R"======(// Vector math helpers for double3
)======"
R"======(__device__ inline double3 operator-(const double3& a, const double3& b) {
)======"
R"======(    return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
)======"
R"======(}
)======"
R"======(
)======"
R"======(__device__ inline double dot(const double3& a, const double3& b) {
)======"
R"======(    return a.x * b.x + a.y * b.y + a.z * b.z;
)======"
R"======(}
)======"
R"======(
)======"
R"======(__device__ inline double3 cross(const double3& a, const double3& b) {
)======"
R"======(    return make_double3(
)======"
R"======(        a.y * b.z - a.z * b.y,
)======"
R"======(        a.z * b.x - a.x * b.z,
)======"
R"======(        a.x * b.y - a.y * b.x
)======"
R"======(    );
)======"
R"======(}
)======"
R"======(
)======"
R"======(__device__ inline double norm(const double3& a) {
)======"
R"======(    return sqrt(dot(a, a));
)======"
R"======(}
)======"
R"======(
)======"
R"======(__device__ inline double3 normalize(const double3& a) {
)======"
R"======(    double n = norm(a);
)======"
R"======(    return make_double3(a.x / n, a.y / n, a.z / n);
)======"
R"======(}
)======"
R"======(
)======"
R"======(__device__ void invert_matrix(const double3 box[3], double3 inverse[3]) {
)======"
R"======(    double a = box[0].x, b = box[0].y, c = box[0].z;
)======"
R"======(    double d = box[1].x, e = box[1].y, f = box[1].z;
)======"
R"======(    double g = box[2].x, h = box[2].y, i = box[2].z;
)======"
R"======(
)======"
R"======(    double det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
)======"
R"======(    double invdet = 1.0 / det;
)======"
R"======(
)======"
R"======(    inverse[0] = make_double3(
)======"
R"======(        (e * i - f * h) * invdet,
)======"
R"======(        (c * h - b * i) * invdet,
)======"
R"======(        (b * f - c * e) * invdet
)======"
R"======(    );
)======"
R"======(    inverse[1] = make_double3(
)======"
R"======(        (f * g - d * i) * invdet,
)======"
R"======(        (a * i - c * g) * invdet,
)======"
R"======(        (c * d - a * f) * invdet
)======"
R"======(    );
)======"
R"======(    inverse[2] = make_double3(
)======"
R"======(        (d * h - e * g) * invdet,
)======"
R"======(        (b * g - a * h) * invdet,
)======"
R"======(        (a * e - b * d) * invdet
)======"
R"======(    );
)======"
R"======(}
)======"
R"======(
)======"
R"======(// Helper to compute Cartesian vector from fractional coordinates
)======"
R"======(// Using row convention: cart = frac @ box (frac as row vector times box matrix)
)======"
R"======(// cart[j] = sum_i(frac[i] * box[i].j)
)======"
R"======(__device__ inline double3 frac_to_cart(const double3& frac, const double3 box[3]) {
)======"
R"======(    return make_double3(
)======"
R"======(        frac.x * box[0].x + frac.y * box[1].x + frac.z * box[2].x,
)======"
R"======(        frac.x * box[0].y + frac.y * box[1].y + frac.z * box[2].y,
)======"
R"======(        frac.x * box[0].z + frac.y * box[1].z + frac.z * box[2].z
)======"
R"======(    );
)======"
R"======(}
)======"
R"======(
)======"
R"======(__device__ void apply_periodic_boundary(
)======"
R"======(    double3& vector,
)======"
R"======(    int3& shift,
)======"
R"======(    const double3 box[3],
)======"
R"======(    const double3 inv_box[3],
)======"
R"======(    const bool* periodic,
)======"
R"======(    bool is_orthogonal
)======"
R"======() {
)======"
R"======(    // Compute fractional coordinates using row convention: frac = vector @ inv_box
)======"
R"======(    // frac[i] = sum_j(vector[j] * inv_box[j].i)
)======"
R"======(    double3 fractional = make_double3(
)======"
R"======(        vector.x * inv_box[0].x + vector.y * inv_box[1].x + vector.z * inv_box[2].x,
)======"
R"======(        vector.x * inv_box[0].y + vector.y * inv_box[1].y + vector.z * inv_box[2].y,
)======"
R"======(        vector.x * inv_box[0].z + vector.y * inv_box[1].z + vector.z * inv_box[2].z
)======"
R"======(    );
)======"
R"======(
)======"
R"======(    // Compute the initial wrapping to bring fractional coords into [-0.5, 0.5]
)======"
R"======(    // The multiplication by `periodic` sets the wrap to zero for non-periodic directions
)======"
R"======(    int3 wrap = make_int3(
)======"
R"======(        static_cast<int>(periodic[0]) * static_cast<int>(round(fractional.x)),
)======"
R"======(        static_cast<int>(periodic[1]) * static_cast<int>(round(fractional.y)),
)======"
R"======(        static_cast<int>(periodic[2]) * static_cast<int>(round(fractional.z))
)======"
R"======(    );
)======"
R"======(
)======"
R"======(    if (!is_orthogonal) {
)======"
R"======(        // For non-orthogonal cells, simple rounding may not find the true minimum image.
)======"
R"======(        // Search all 27 neighboring images to find the one with minimum distance.
)======"
R"======(        double min_dist2 = 1e30;
)======"
R"======(        int3 best_wrap = wrap;
)======"
R"======(
)======"
R"======(        for (int dx = -1; dx <= 1; dx++) {
)======"
R"======(            for (int dy = -1; dy <= 1; dy++) {
)======"
R"======(                for (int dz = -1; dz <= 1; dz++) {
)======"
R"======(                    int3 test_wrap = make_int3(
)======"
R"======(                        (wrap.x + dx) * static_cast<int>(periodic[0]),
)======"
R"======(                        (wrap.y + dy) * static_cast<int>(periodic[1]),
)======"
R"======(                        (wrap.z + dz) * static_cast<int>(periodic[2])
)======"
R"======(                    );
)======"
R"======(
)======"
R"======(                    double3 test_frac = make_double3(
)======"
R"======(                        fractional.x - test_wrap.x,
)======"
R"======(                        fractional.y - test_wrap.y,
)======"
R"======(                        fractional.z - test_wrap.z
)======"
R"======(                    );
)======"
R"======(
)======"
R"======(                    double3 test_vec = frac_to_cart(test_frac, box);
)======"
R"======(                    double dist2 = dot(test_vec, test_vec);
)======"
R"======(
)======"
R"======(                    if (dist2 < min_dist2) {
)======"
R"======(                        min_dist2 = dist2;
)======"
R"======(                        best_wrap = test_wrap;
)======"
R"======(                    }
)======"
R"======(                }
)======"
R"======(            }
)======"
R"======(        }
)======"
R"======(        wrap = best_wrap;
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    // The stored shift follows the convention: vector = rj - ri + shift @ box
)======"
R"======(    // Since we compute wrapped = vector - wrap @ box, the shift is -wrap
)======"
R"======(    shift = make_int3(-wrap.x, -wrap.y, -wrap.z);
)======"
R"======(
)======"
R"======(    fractional = make_double3(
)======"
R"======(        fractional.x - wrap.x,
)======"
R"======(        fractional.y - wrap.y,
)======"
R"======(        fractional.z - wrap.z
)======"
R"======(    );
)======"
R"======(
)======"
R"======(    vector = frac_to_cart(fractional, box);
)======"
R"======(}
)======"
R"======(
)======"
R"======(__global__ void compute_mic_neighbours_full_impl(
)======"
R"======(    const double* positions,
)======"
R"======(    const double* box,
)======"
R"======(    const bool* periodic,
)======"
R"======(    size_t n_points,
)======"
R"======(    double cutoff,
)======"
R"======(    size_t* length,
)======"
R"======(    size_t* pair_indices,
)======"
R"======(    int* shifts,
)======"
R"======(    double* distances,
)======"
R"======(    double* vectors,
)======"
R"======(    bool return_shifts,
)======"
R"======(    bool return_distances,
)======"
R"======(    bool return_vectors
)======"
R"======() {
)======"
R"======(    __shared__ double3 shared_box[3];
)======"
R"======(    __shared__ double3 shared_inv_box[3];
)======"
R"======(    __shared__ bool shared_is_orthogonal;
)======"
R"======(
)======"
R"======(    const int warp_id = threadIdx.x / WARP_SIZE;
)======"
R"======(    const int thread_id = threadIdx.x % WARP_SIZE;
)======"
R"======(
)======"
R"======(    const size_t point_i = blockIdx.x * NWARPS + warp_id;
)======"
R"======(    const double cutoff2 = cutoff * cutoff;
)======"
R"======(
)======"
R"======(    // Load current box to shared memory
)======"
R"======(    if (threadIdx.x < 3) {
)======"
R"======(        shared_box[threadIdx.x] = make_double3(
)======"
R"======(            box[threadIdx.x * 3],
)======"
R"======(            box[threadIdx.x * 3 + 1],
)======"
R"======(            box[threadIdx.x * 3 + 2]
)======"
R"======(        );
)======"
R"======(    }
)======"
R"======(    __syncthreads();
)======"
R"======(
)======"
R"======(    // Overwrite non-periodic directions with unit vectors orthogonal to the periodic subspace
)======"
R"======(    if (threadIdx.x == 0) {
)======"
R"======(        // Collect periodic / non-periodic indices
)======"
R"======(        int n_periodic = 0;
)======"
R"======(        int periodic_idx_1 = -1;
)======"
R"======(        int periodic_idx_2 = -1;
)======"
R"======(        for (int i = 0; i < 3; ++i) {
)======"
R"======(            if (periodic[i]) {
)======"
R"======(                n_periodic += 1;
)======"
R"======(                if (periodic_idx_1 == -1) {
)======"
R"======(                    periodic_idx_1 = i;
)======"
R"======(                } else if (periodic_idx_2 == -1) {
)======"
R"======(                    periodic_idx_2 = i;
)======"
R"======(                }
)======"
R"======(            }
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        if (n_periodic == 0) {
)======"
R"======(            // Fully non-periodic: any orthonormal basis is fine
)======"
R"======(            shared_box[0] = make_double3(1.0, 0.0, 0.0);
)======"
R"======(            shared_box[1] = make_double3(0.0, 1.0, 0.0);
)======"
R"======(            shared_box[2] = make_double3(0.0, 0.0, 1.0);
)======"
R"======(        } else if (n_periodic == 1) {
)======"
R"======(            // 1D periodic: build an orthonormal pair spanning the plane orthogonal to the periodic vector
)======"
R"======(            double3 a = shared_box[periodic_idx_1];
)======"
R"======(            double3 b = make_double3(0, 1, 0);
)======"
R"======(            if (fabs(dot(normalize(a), b)) > 0.9) {
)======"
R"======(                b = make_double3(0, 0, 1);
)======"
R"======(            }
)======"
R"======(            double3 c = normalize(cross(a, b));
)======"
R"======(            b = normalize(cross(c, a));
)======"
R"======(
)======"
R"======(            shared_box[(periodic_idx_1 + 1) % 3] = b;
)======"
R"======(            shared_box[(periodic_idx_1 + 2) % 3] = c;
)======"
R"======(        } else if (n_periodic == 2) {
)======"
R"======(            // 2D periodic: set the sole non-periodic direction to the plane normal
)======"
R"======(            double3 a = shared_box[periodic_idx_1];
)======"
R"======(            double3 b = shared_box[periodic_idx_2];
)======"
R"======(            double3 c = normalize(cross(a, b));
)======"
R"======(
)======"
R"======(            int non_periodic_idx = 3 - periodic_idx_1 - periodic_idx_2;
)======"
R"======(            shared_box[non_periodic_idx] = c;
)======"
R"======(        }
)======"
R"======(        // n_periodic == 3: fully periodic, keep shared_box as-is
)======"
R"======(
)======"
R"======(        invert_matrix(shared_box, shared_inv_box);
)======"
R"======(
)======"
R"======(        // Check orthogonality: all off-diagonal dot products should be ~0
)======"
R"======(        double tol = 1e-10;
)======"
R"======(        double ab = fabs(dot(shared_box[0], shared_box[1]));
)======"
R"======(        double ac = fabs(dot(shared_box[0], shared_box[2]));
)======"
R"======(        double bc = fabs(dot(shared_box[1], shared_box[2]));
)======"
R"======(        shared_is_orthogonal = (ab < tol) && (ac < tol) && (bc < tol);
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    // Ensure inv_box and is_orthogonal are ready
)======"
R"======(    __syncthreads();
)======"
R"======(
)======"
R"======(    if (point_i >= n_points) {
)======"
R"======(        return;
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    bool is_orthogonal = shared_is_orthogonal;
)======"
R"======(    double3 ri = make_double3(
)======"
R"======(        positions[point_i * 3],
)======"
R"======(        positions[point_i * 3 + 1],
)======"
R"======(        positions[point_i * 3 + 2]
)======"
R"======(    );
)======"
R"======(
)======"
R"======(    for (size_t j = thread_id; j < n_points; j += WARP_SIZE) {
)======"
R"======(        double3 rj = make_double3(
)======"
R"======(            positions[j * 3],
)======"
R"======(            positions[j * 3 + 1],
)======"
R"======(            positions[j * 3 + 2]
)======"
R"======(        );
)======"
R"======(
)======"
R"======(        double3 vector = rj - ri;
)======"
R"======(        int3 shift = make_int3(0, 0, 0);
)======"
R"======(        apply_periodic_boundary(vector, shift, shared_box, shared_inv_box, periodic, is_orthogonal);
)======"
R"======(
)======"
R"======(        double distance2 = dot(vector, vector);
)======"
R"======(        bool is_valid = (distance2 < cutoff2 && distance2 > 0.0);
)======"
R"======(
)======"
R"======(        if (is_valid) {
)======"
R"======(            size_t current_pair = atomicAdd_size_t(&length[0], 1);
)======"
R"======(            pair_indices[current_pair * 2] = point_i;
)======"
R"======(            pair_indices[current_pair * 2 + 1] = j;
)======"
R"======(
)======"
R"======(            if (return_shifts) {
)======"
R"======(                shifts[current_pair * 3] = shift.x;
)======"
R"======(                shifts[current_pair * 3 + 1] = shift.y;
)======"
R"======(                shifts[current_pair * 3 + 2] = shift.z;
)======"
R"======(            }
)======"
R"======(            if (return_vectors) {
)======"
R"======(                vectors[current_pair * 3] = vector.x;
)======"
R"======(                vectors[current_pair * 3 + 1] = vector.y;
)======"
R"======(                vectors[current_pair * 3 + 2] = vector.z;
)======"
R"======(            }
)======"
R"======(            if (return_distances) {
)======"
R"======(                distances[current_pair] = sqrt(distance2);
)======"
R"======(            }
)======"
R"======(        }
)======"
R"======(    }
)======"
R"======(}
)======"
R"======(
)======"
R"======(__global__ void compute_mic_neighbours_half_impl(
)======"
R"======(    const double* positions,
)======"
R"======(    const double* box,
)======"
R"======(    const bool* periodic,
)======"
R"======(    size_t n_points,
)======"
R"======(    double cutoff,
)======"
R"======(    size_t* length,
)======"
R"======(    size_t* pair_indices,
)======"
R"======(    int* shifts,
)======"
R"======(    double* distances,
)======"
R"======(    double* vectors,
)======"
R"======(    bool return_shifts,
)======"
R"======(    bool return_distances,
)======"
R"======(    bool return_vectors
)======"
R"======() {
)======"
R"======(    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
)======"
R"======(    const size_t num_all_pairs = n_points * (n_points - 1) / 2;
)======"
R"======(    const double cutoff2 = cutoff * cutoff;
)======"
R"======(
)======"
R"======(    __shared__ double3 shared_box[3];
)======"
R"======(    __shared__ double3 shared_inv_box[3];
)======"
R"======(    __shared__ bool shared_is_orthogonal;
)======"
R"======(
)======"
R"======(    // Load current box to shared memory
)======"
R"======(    if (threadIdx.x < 3) {
)======"
R"======(        shared_box[threadIdx.x] = make_double3(
)======"
R"======(            box[threadIdx.x * 3],
)======"
R"======(            box[threadIdx.x * 3 + 1],
)======"
R"======(            box[threadIdx.x * 3 + 2]
)======"
R"======(        );
)======"
R"======(    }
)======"
R"======(    __syncthreads();
)======"
R"======(
)======"
R"======(    // Overwrite non-periodic directions with unit vectors orthogonal to the periodic subspace
)======"
R"======(    if (threadIdx.x == 0) {
)======"
R"======(        int n_periodic = 0;
)======"
R"======(        int periodic_idx_1 = -1;
)======"
R"======(        int periodic_idx_2 = -1;
)======"
R"======(        for (int i = 0; i < 3; ++i) {
)======"
R"======(            if (periodic[i]) {
)======"
R"======(                n_periodic += 1;
)======"
R"======(                if (periodic_idx_1 == -1) {
)======"
R"======(                    periodic_idx_1 = i;
)======"
R"======(                } else if (periodic_idx_2 == -1) {
)======"
R"======(                    periodic_idx_2 = i;
)======"
R"======(                }
)======"
R"======(            }
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        if (n_periodic == 0) {
)======"
R"======(            shared_box[0] = make_double3(1.0, 0.0, 0.0);
)======"
R"======(            shared_box[1] = make_double3(0.0, 1.0, 0.0);
)======"
R"======(            shared_box[2] = make_double3(0.0, 0.0, 1.0);
)======"
R"======(        } else if (n_periodic == 1) {
)======"
R"======(            double3 a = shared_box[periodic_idx_1];
)======"
R"======(            double3 b = make_double3(0, 1, 0);
)======"
R"======(            if (fabs(dot(normalize(a), b)) > 0.9) {
)======"
R"======(                b = make_double3(0, 0, 1);
)======"
R"======(            }
)======"
R"======(            double3 c = normalize(cross(a, b));
)======"
R"======(            b = normalize(cross(c, a));
)======"
R"======(
)======"
R"======(            shared_box[(periodic_idx_1 + 1) % 3] = b;
)======"
R"======(            shared_box[(periodic_idx_1 + 2) % 3] = c;
)======"
R"======(        } else if (n_periodic == 2) {
)======"
R"======(            double3 a = shared_box[periodic_idx_1];
)======"
R"======(            double3 b = shared_box[periodic_idx_2];
)======"
R"======(            double3 c = normalize(cross(a, b));
)======"
R"======(
)======"
R"======(            int non_periodic_idx = 3 - periodic_idx_1 - periodic_idx_2;
)======"
R"======(            shared_box[non_periodic_idx] = c;
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        invert_matrix(shared_box, shared_inv_box);
)======"
R"======(
)======"
R"======(        double tol = 1e-10;
)======"
R"======(        double ab = fabs(dot(shared_box[0], shared_box[1]));
)======"
R"======(        double ac = fabs(dot(shared_box[0], shared_box[2]));
)======"
R"======(        double bc = fabs(dot(shared_box[1], shared_box[2]));
)======"
R"======(        shared_is_orthogonal = (ab < tol) && (ac < tol) && (bc < tol);
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    __syncthreads();
)======"
R"======(
)======"
R"======(    if (index >= num_all_pairs) {
)======"
R"======(        return;
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    bool is_orthogonal = shared_is_orthogonal;
)======"
R"======(
)======"
R"======(    size_t point_j = floor((sqrt(8.0 * index + 1.0) + 1.0) / 2.0);
)======"
R"======(    if (point_j * (point_j - 1) > 2 * index) {
)======"
R"======(        point_j--;
)======"
R"======(    }
)======"
R"======(    const size_t point_i = index - point_j * (point_j - 1) / 2;
)======"
R"======(
)======"
R"======(    double3 ri = make_double3(
)======"
R"======(        positions[point_i * 3],
)======"
R"======(        positions[point_i * 3 + 1],
)======"
R"======(        positions[point_i * 3 + 2]
)======"
R"======(    );
)======"
R"======(    double3 rj = make_double3(
)======"
R"======(        positions[point_j * 3],
)======"
R"======(        positions[point_j * 3 + 1],
)======"
R"======(        positions[point_j * 3 + 2]
)======"
R"======(    );
)======"
R"======(
)======"
R"======(    double3 vector = rj - ri;
)======"
R"======(    int3 shift = make_int3(0, 0, 0);
)======"
R"======(    apply_periodic_boundary(vector, shift, shared_box, shared_inv_box, periodic, is_orthogonal);
)======"
R"======(
)======"
R"======(    double distance2 = dot(vector, vector);
)======"
R"======(    bool is_valid = (distance2 < cutoff2 && distance2 > 0.0);
)======"
R"======(
)======"
R"======(    if (is_valid) {
)======"
R"======(        size_t pair_index = atomicAdd_size_t(&length[0], 1);
)======"
R"======(        pair_indices[pair_index * 2] = point_i;
)======"
R"======(        pair_indices[pair_index * 2 + 1] = point_j;
)======"
R"======(
)======"
R"======(        if (return_shifts) {
)======"
R"======(            shifts[pair_index * 3] = shift.x;
)======"
R"======(            shifts[pair_index * 3 + 1] = shift.y;
)======"
R"======(            shifts[pair_index * 3 + 2] = shift.z;
)======"
R"======(        }
)======"
R"======(        if (return_vectors) {
)======"
R"======(            vectors[pair_index * 3] = vector.x;
)======"
R"======(            vectors[pair_index * 3 + 1] = vector.y;
)======"
R"======(            vectors[pair_index * 3 + 2] = vector.z;
)======"
R"======(        }
)======"
R"======(        if (return_distances) {
)======"
R"======(            distances[pair_index] = sqrt(distance2);
)======"
R"======(        }
)======"
R"======(    }
)======"
R"======(}
)======"
R"======(
)======"
R"======(// ============================================================================
)======"
R"======(// Optimized brute force kernels with precomputed box parameters
)======"
R"======(// These avoid per-block initialization by having inv_box, is_orthogonal passed in
)======"
R"======(// ============================================================================
)======"
R"======(
)======"
R"======(// Simple PBC for orthogonal boxes (most common case)
)======"
R"======(// For non-periodic directions, no wrapping is applied
)======"
R"======(__device__ inline void apply_pbc_orthogonal(
)======"
R"======(    double3& d,
)======"
R"======(    int3& shift,
)======"
R"======(    const double3& box_diag,
)======"
R"======(    const bool* periodic
)======"
R"======() {
)======"
R"======(    shift = make_int3(0, 0, 0);
)======"
R"======(    if (periodic[0] && box_diag.x > 0) {
)======"
R"======(        int s = static_cast<int>(round(d.x / box_diag.x));
)======"
R"======(        d.x -= s * box_diag.x;
)======"
R"======(        shift.x = -s;
)======"
R"======(    }
)======"
R"======(    if (periodic[1] && box_diag.y > 0) {
)======"
R"======(        int s = static_cast<int>(round(d.y / box_diag.y));
)======"
R"======(        d.y -= s * box_diag.y;
)======"
R"======(        shift.y = -s;
)======"
R"======(    }
)======"
R"======(    if (periodic[2] && box_diag.z > 0) {
)======"
R"======(        int s = static_cast<int>(round(d.z / box_diag.z));
)======"
R"======(        d.z -= s * box_diag.z;
)======"
R"======(        shift.z = -s;
)======"
R"======(    }
)======"
R"======(}
)======"
R"======(
)======"
R"======(// General PBC using precomputed inverse box
)======"
R"======(__device__ inline void apply_pbc_general(
)======"
R"======(    double3& vector,
)======"
R"======(    int3& shift,
)======"
R"======(    const double3 box[3],
)======"
R"======(    const double3 inv_box[3],
)======"
R"======(    const bool* periodic
)======"
R"======() {
)======"
R"======(    double3 frac = make_double3(
)======"
R"======(        vector.x * inv_box[0].x + vector.y * inv_box[1].x + vector.z * inv_box[2].x,
)======"
R"======(        vector.x * inv_box[0].y + vector.y * inv_box[1].y + vector.z * inv_box[2].y,
)======"
R"======(        vector.x * inv_box[0].z + vector.y * inv_box[1].z + vector.z * inv_box[2].z
)======"
R"======(    );
)======"
R"======(
)======"
R"======(    int3 wrap = make_int3(
)======"
R"======(        periodic[0] ? static_cast<int>(round(frac.x)) : 0,
)======"
R"======(        periodic[1] ? static_cast<int>(round(frac.y)) : 0,
)======"
R"======(        periodic[2] ? static_cast<int>(round(frac.z)) : 0
)======"
R"======(    );
)======"
R"======(
)======"
R"======(    frac.x -= wrap.x;
)======"
R"======(    frac.y -= wrap.y;
)======"
R"======(    frac.z -= wrap.z;
)======"
R"======(
)======"
R"======(    vector = frac_to_cart(frac, box);
)======"
R"======(
)======"
R"======(    shift = make_int3(-wrap.x, -wrap.y, -wrap.z);
)======"
R"======(}
)======"
R"======(
)======"
R"======(__global__ void brute_force_half_orthogonal(
)======"
R"======(    const double* __restrict__ positions,
)======"
R"======(    const double* __restrict__ box_diag,
)======"
R"======(    const bool* __restrict__ periodic,
)======"
R"======(    size_t n_points,
)======"
R"======(    double cutoff2,
)======"
R"======(    size_t* length,
)======"
R"======(    size_t* pair_indices,
)======"
R"======(    int* shifts,
)======"
R"======(    double* distances,
)======"
R"======(    double* vectors,
)======"
R"======(    bool return_shifts,
)======"
R"======(    bool return_distances,
)======"
R"======(    bool return_vectors,
)======"
R"======(    size_t max_pairs,
)======"
R"======(    int* overflow_flag
)======"
R"======() {
)======"
R"======(    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
)======"
R"======(    const size_t num_all_pairs = n_points * (n_points - 1) / 2;
)======"
R"======(
)======"
R"======(    if (index >= num_all_pairs) {
)======"
R"======(        return;
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    size_t j = static_cast<size_t>(floor((sqrt(8.0 * index + 1.0) + 1.0) / 2.0));
)======"
R"======(    if (j * (j - 1) > 2 * index) {
)======"
R"======(        j--;
)======"
R"======(    }
)======"
R"======(    const size_t i = index - j * (j - 1) / 2;
)======"
R"======(
)======"
R"======(    const double3* pos3 = reinterpret_cast<const double3*>(positions);
)======"
R"======(    double3 pi = pos3[i];
)======"
R"======(    double3 pj = pos3[j];
)======"
R"======(    double3 d = pj - pi;
)======"
R"======(    double3 L = make_double3(box_diag[0], box_diag[1], box_diag[2]);
)======"
R"======(
)======"
R"======(    int3 s;
)======"
R"======(    apply_pbc_orthogonal(d, s, L, periodic);
)======"
R"======(
)======"
R"======(    double dist2 = dot(d, d);
)======"
R"======(
)======"
R"======(    if (dist2 < cutoff2 && dist2 > 0.0) {
)======"
R"======(        size_t idx = atomicAdd_size_t(length, 1UL);
)======"
R"======(
)======"
R"======(        // Check if we are about to exceed max_pairs
)======"
R"======(        if (idx + 1 > max_pairs) {
)======"
R"======(            atomicExch(overflow_flag, 1);
)======"
R"======(            return;
)======"
R"======(        }
)======"
R"======(        pair_indices[idx * 2] = i;
)======"
R"======(        pair_indices[idx * 2 + 1] = j;
)======"
R"======(        if (return_shifts) {
)======"
R"======(            shifts[idx * 3] = s.x;
)======"
R"======(            shifts[idx * 3 + 1] = s.y;
)======"
R"======(            shifts[idx * 3 + 2] = s.z;
)======"
R"======(        }
)======"
R"======(        if (return_vectors) {
)======"
R"======(            vectors[idx * 3] = d.x;
)======"
R"======(            vectors[idx * 3 + 1] = d.y;
)======"
R"======(            vectors[idx * 3 + 2] = d.z;
)======"
R"======(        }
)======"
R"======(        if (return_distances) {
)======"
R"======(            distances[idx] = sqrt(dist2);
)======"
R"======(        }
)======"
R"======(    }
)======"
R"======(}
)======"
R"======(
)======"
R"======(// Triangular indexing: one thread per unordered pair, outputs both (i,j) and (j,i)
)======"
R"======(__global__ void brute_force_full_orthogonal(
)======"
R"======(    const double* __restrict__ positions,
)======"
R"======(    const double* __restrict__ box_diag,
)======"
R"======(    const bool* __restrict__ periodic,
)======"
R"======(    size_t n_points,
)======"
R"======(    double cutoff2,
)======"
R"======(    size_t* length,
)======"
R"======(    size_t* pair_indices,
)======"
R"======(    int* shifts,
)======"
R"======(    double* distances,
)======"
R"======(    double* vectors,
)======"
R"======(    bool return_shifts,
)======"
R"======(    bool return_distances,
)======"
R"======(    bool return_vectors,
)======"
R"======(    size_t max_pairs,
)======"
R"======(    int* overflow_flag
)======"
R"======() {
)======"
R"======(    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
)======"
R"======(    const size_t num_half_pairs = n_points * (n_points - 1) / 2;
)======"
R"======(
)======"
R"======(    if (index >= num_half_pairs) {
)======"
R"======(        return;
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    size_t j = static_cast<size_t>(floor((sqrt(8.0 * index + 1.0) + 1.0) / 2.0));
)======"
R"======(    if (j * (j - 1) > 2 * index) {
)======"
R"======(        j--;
)======"
R"======(    }
)======"
R"======(    const size_t i = index - j * (j - 1) / 2;
)======"
R"======(
)======"
R"======(    const double3* pos3 = reinterpret_cast<const double3*>(positions);
)======"
R"======(    double3 pi = pos3[i];
)======"
R"======(    double3 pj = pos3[j];
)======"
R"======(    double3 d = pj - pi;
)======"
R"======(    double3 L = make_double3(box_diag[0], box_diag[1], box_diag[2]);
)======"
R"======(
)======"
R"======(    int3 s;
)======"
R"======(    apply_pbc_orthogonal(d, s, L, periodic);
)======"
R"======(
)======"
R"======(    double dist2 = dot(d, d);
)======"
R"======(
)======"
R"======(    if (dist2 < cutoff2) {
)======"
R"======(        size_t idx = atomicAdd_size_t(length, 2UL);
)======"
R"======(
)======"
R"======(        // Check if we are about to exceed max_pairs
)======"
R"======(        if (idx + 2 > max_pairs) {
)======"
R"======(            atomicExch(overflow_flag, 1);
)======"
R"======(            return;
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        pair_indices[idx * 2] = i;
)======"
R"======(        pair_indices[idx * 2 + 1] = j;
)======"
R"======(        pair_indices[(idx + 1) * 2] = j;
)======"
R"======(        pair_indices[(idx + 1) * 2 + 1] = i;
)======"
R"======(
)======"
R"======(        if (return_shifts) {
)======"
R"======(            shifts[idx * 3] = s.x;
)======"
R"======(            shifts[idx * 3 + 1] = s.y;
)======"
R"======(            shifts[idx * 3 + 2] = s.z;
)======"
R"======(            shifts[(idx + 1) * 3] = -s.x;
)======"
R"======(            shifts[(idx + 1) * 3 + 1] = -s.y;
)======"
R"======(            shifts[(idx + 1) * 3 + 2] = -s.z;
)======"
R"======(        }
)======"
R"======(        if (return_vectors) {
)======"
R"======(            vectors[idx * 3] = d.x;
)======"
R"======(            vectors[idx * 3 + 1] = d.y;
)======"
R"======(            vectors[idx * 3 + 2] = d.z;
)======"
R"======(            vectors[(idx + 1) * 3] = -d.x;
)======"
R"======(            vectors[(idx + 1) * 3 + 1] = -d.y;
)======"
R"======(            vectors[(idx + 1) * 3 + 2] = -d.z;
)======"
R"======(        }
)======"
R"======(        if (return_distances) {
)======"
R"======(            double dist = sqrt(dist2);
)======"
R"======(            distances[idx] = dist;
)======"
R"======(            distances[idx + 1] = dist;
)======"
R"======(        }
)======"
R"======(    }
)======"
R"======(}
)======"
R"======(
)======"
R"======(__global__ void brute_force_half_general(
)======"
R"======(    const double* __restrict__ positions,
)======"
R"======(    const double* __restrict__ box,
)======"
R"======(    const double* __restrict__ inv_box,
)======"
R"======(    const bool* __restrict__ periodic,
)======"
R"======(    size_t n_points,
)======"
R"======(    double cutoff2,
)======"
R"======(    size_t* length,
)======"
R"======(    size_t* pair_indices,
)======"
R"======(    int* shifts,
)======"
R"======(    double* distances,
)======"
R"======(    double* vectors,
)======"
R"======(    bool return_shifts,
)======"
R"======(    bool return_distances,
)======"
R"======(    bool return_vectors,
)======"
R"======(    size_t max_pairs,
)======"
R"======(    int* overflow_flag
)======"
R"======() {
)======"
R"======(    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
)======"
R"======(    const size_t num_all_pairs = n_points * (n_points - 1) / 2;
)======"
R"======(
)======"
R"======(    // Load box into double3 arrays
)======"
R"======(    __shared__ double3 shared_box[3];
)======"
R"======(    __shared__ double3 shared_inv_box[3];
)======"
R"======(
)======"
R"======(    if (threadIdx.x < 3) {
)======"
R"======(        shared_box[threadIdx.x] = make_double3(
)======"
R"======(            box[threadIdx.x * 3],
)======"
R"======(            box[threadIdx.x * 3 + 1],
)======"
R"======(            box[threadIdx.x * 3 + 2]
)======"
R"======(        );
)======"
R"======(        shared_inv_box[threadIdx.x] = make_double3(
)======"
R"======(            inv_box[threadIdx.x * 3],
)======"
R"======(            inv_box[threadIdx.x * 3 + 1],
)======"
R"======(            inv_box[threadIdx.x * 3 + 2]
)======"
R"======(        );
)======"
R"======(    }
)======"
R"======(    __syncthreads();
)======"
R"======(
)======"
R"======(    if (index >= num_all_pairs) {
)======"
R"======(        return;
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    size_t j = static_cast<size_t>(floor((sqrt(8.0 * index + 1.0) + 1.0) / 2.0));
)======"
R"======(    if (j * (j - 1) > 2 * index) {
)======"
R"======(        j--;
)======"
R"======(    }
)======"
R"======(    const size_t i = index - j * (j - 1) / 2;
)======"
R"======(
)======"
R"======(    const double3* pos3 = reinterpret_cast<const double3*>(positions);
)======"
R"======(    double3 pi = pos3[i];
)======"
R"======(    double3 pj = pos3[j];
)======"
R"======(    double3 vector = pj - pi;
)======"
R"======(    int3 shift;
)======"
R"======(
)======"
R"======(    apply_pbc_general(vector, shift, shared_box, shared_inv_box, periodic);
)======"
R"======(
)======"
R"======(    double dist2 = dot(vector, vector);
)======"
R"======(
)======"
R"======(    if (dist2 < cutoff2 && dist2 > 0.0) {
)======"
R"======(        size_t idx = atomicAdd_size_t(length, 1UL);
)======"
R"======(
)======"
R"======(        // Check if we are about to exceed max_pairs
)======"
R"======(        if (idx + 1 > max_pairs) {
)======"
R"======(            atomicExch(overflow_flag, 1);
)======"
R"======(            return;
)======"
R"======(        }
)======"
R"======(        pair_indices[idx * 2] = i;
)======"
R"======(        pair_indices[idx * 2 + 1] = j;
)======"
R"======(        if (return_shifts) {
)======"
R"======(            shifts[idx * 3] = shift.x;
)======"
R"======(            shifts[idx * 3 + 1] = shift.y;
)======"
R"======(            shifts[idx * 3 + 2] = shift.z;
)======"
R"======(        }
)======"
R"======(        if (return_vectors) {
)======"
R"======(            vectors[idx * 3] = vector.x;
)======"
R"======(            vectors[idx * 3 + 1] = vector.y;
)======"
R"======(            vectors[idx * 3 + 2] = vector.z;
)======"
R"======(        }
)======"
R"======(        if (return_distances) {
)======"
R"======(            distances[idx] = sqrt(dist2);
)======"
R"======(        }
)======"
R"======(    }
)======"
R"======(}
)======"
R"======(
)======"
R"======(// Optimized full-list kernel for GENERAL boxes
)======"
R"======(// NNPOps-style triangular indexing: one thread per unordered pair, outputs both (i,j) and (j,i)
)======"
R"======(// Uses double3 for vectorized position loads
)======"
R"======(__global__ void brute_force_full_general(
)======"
R"======(    const double* __restrict__ positions,
)======"
R"======(    const double* __restrict__ box,
)======"
R"======(    const double* __restrict__ inv_box,
)======"
R"======(    const bool* __restrict__ periodic,
)======"
R"======(    size_t n_points,
)======"
R"======(    double cutoff2,
)======"
R"======(    size_t* length,
)======"
R"======(    size_t* pair_indices,
)======"
R"======(    int* shifts,
)======"
R"======(    double* distances,
)======"
R"======(    double* vectors,
)======"
R"======(    bool return_shifts,
)======"
R"======(    bool return_distances,
)======"
R"======(    bool return_vectors,
)======"
R"======(    size_t max_pairs,
)======"
R"======(    int* overflow_flag
)======"
R"======() {
)======"
R"======(    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
)======"
R"======(    const size_t num_half_pairs = n_points * (n_points - 1) / 2;
)======"
R"======(
)======"
R"======(    // Load box into double3 arrays
)======"
R"======(    __shared__ double3 shared_box[3];
)======"
R"======(    __shared__ double3 shared_inv_box[3];
)======"
R"======(
)======"
R"======(    if (threadIdx.x < 3) {
)======"
R"======(        shared_box[threadIdx.x] = make_double3(
)======"
R"======(            box[threadIdx.x * 3],
)======"
R"======(            box[threadIdx.x * 3 + 1],
)======"
R"======(            box[threadIdx.x * 3 + 2]
)======"
R"======(        );
)======"
R"======(        shared_inv_box[threadIdx.x] = make_double3(
)======"
R"======(            inv_box[threadIdx.x * 3],
)======"
R"======(            inv_box[threadIdx.x * 3 + 1],
)======"
R"======(            inv_box[threadIdx.x * 3 + 2]
)======"
R"======(        );
)======"
R"======(    }
)======"
R"======(    __syncthreads();
)======"
R"======(
)======"
R"======(    if (index >= num_half_pairs) {
)======"
R"======(        return;
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    // NNPOps-style triangular indexing for half-list
)======"
R"======(    size_t j = static_cast<size_t>(floor((sqrt(8.0 * index + 1.0) + 1.0) / 2.0));
)======"
R"======(    if (j * (j - 1) > 2 * index) {
)======"
R"======(        j--;
)======"
R"======(    }
)======"
R"======(    const size_t i = index - j * (j - 1) / 2;
)======"
R"======(
)======"
R"======(    const double3* pos3 = reinterpret_cast<const double3*>(positions);
)======"
R"======(    double3 pi = pos3[i];
)======"
R"======(    double3 pj = pos3[j];
)======"
R"======(    double3 vector = pj - pi;
)======"
R"======(    int3 shift;
)======"
R"======(
)======"
R"======(    apply_pbc_general(vector, shift, shared_box, shared_inv_box, periodic);
)======"
R"======(
)======"
R"======(    double dist2 = dot(vector, vector);
)======"
R"======(
)======"
R"======(    if (dist2 < cutoff2) {
)======"
R"======(        size_t idx = atomicAdd_size_t(length, 2UL);
)======"
R"======(
)======"
R"======(        // Check if we are about to exceed max_pairs
)======"
R"======(        if (idx + 2 > max_pairs) {
)======"
R"======(            atomicExch(overflow_flag, 1);
)======"
R"======(            return;
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        pair_indices[idx * 2] = i;
)======"
R"======(        pair_indices[idx * 2 + 1] = j;
)======"
R"======(        pair_indices[(idx + 1) * 2] = j;
)======"
R"======(        pair_indices[(idx + 1) * 2 + 1] = i;
)======"
R"======(
)======"
R"======(        if (return_shifts) {
)======"
R"======(            shifts[idx * 3] = shift.x;
)======"
R"======(            shifts[idx * 3 + 1] = shift.y;
)======"
R"======(            shifts[idx * 3 + 2] = shift.z;
)======"
R"======(            shifts[(idx + 1) * 3] = -shift.x;
)======"
R"======(            shifts[(idx + 1) * 3 + 1] = -shift.y;
)======"
R"======(            shifts[(idx + 1) * 3 + 2] = -shift.z;
)======"
R"======(        }
)======"
R"======(        if (return_vectors) {
)======"
R"======(            vectors[idx * 3] = vector.x;
)======"
R"======(            vectors[idx * 3 + 1] = vector.y;
)======"
R"======(            vectors[idx * 3 + 2] = vector.z;
)======"
R"======(            vectors[(idx + 1) * 3] = -vector.x;
)======"
R"======(            vectors[(idx + 1) * 3 + 1] = -vector.y;
)======"
R"======(            vectors[(idx + 1) * 3 + 2] = -vector.z;
)======"
R"======(        }
)======"
R"======(        if (return_distances) {
)======"
R"======(            double dist = sqrt(dist2);
)======"
R"======(            distances[idx] = dist;
)======"
R"======(            distances[idx + 1] = dist;
)======"
R"======(        }
)======"
R"======(    }
)======"
R"======(}
)======"
R"======(
)======"
R"======(// Status flags for mic_box_check
)======"
R"======(// bit 0: error (cutoff too large)
)======"
R"======(// bit 1: is_orthogonal
)======"
R"======(#define BOX_STATUS_ERROR 1
)======"
R"======(#define BOX_STATUS_ORTHOGONAL 2
)======"
R"======(
)======"
R"======(__global__ void mic_box_check(
)======"
R"======(    const double* box,
)======"
R"======(    const bool* periodic,
)======"
R"======(    const double cutoff,
)======"
R"======(    int* status,
)======"
R"======(    double* box_diag,   // Output: [Lx, Ly, Lz] for orthogonal boxes (can be nullptr)
)======"
R"======(    double* inv_box_out // Output: 9-element inverse box matrix (can be nullptr)
)======"
R"======() {
)======"
R"======(    __shared__ double3 shared_box[3];
)======"
R"======(
)======"
R"======(    if (threadIdx.x < 3) {
)======"
R"======(        shared_box[threadIdx.x] = make_double3(
)======"
R"======(            box[threadIdx.x * 3],
)======"
R"======(            box[threadIdx.x * 3 + 1],
)======"
R"======(            box[threadIdx.x * 3 + 2]
)======"
R"======(        );
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    __syncthreads();
)======"
R"======(
)======"
R"======(    if (threadIdx.x == 0) {
)======"
R"======(        double3 a = shared_box[0];
)======"
R"======(        double3 b = shared_box[1];
)======"
R"======(        double3 c = shared_box[2];
)======"
R"======(
)======"
R"======(        double a_norm = norm(a);
)======"
R"======(        double b_norm = norm(b);
)======"
R"======(        double c_norm = norm(c);
)======"
R"======(
)======"
R"======(        // Count periodic directions
)======"
R"======(        int n_periodic = 0;
)======"
R"======(        if (periodic[0]) {
)======"
R"======(            n_periodic++;
)======"
R"======(        }
)======"
R"======(        if (periodic[1]) {
)======"
R"======(            n_periodic++;
)======"
R"======(        }
)======"
R"======(        if (periodic[2]) {
)======"
R"======(            n_periodic++;
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        double ab_dot = dot(a, b);
)======"
R"======(        double ac_dot = dot(a, c);
)======"
R"======(        double bc_dot = dot(b, c);
)======"
R"======(
)======"
R"======(        double tol = 1e-6;
)======"
R"======(        // Treat fully non-periodic systems as orthogonal (no PBC needed)
)======"
R"======(        // Also treat systems with zero-norm vectors as orthogonal (degenerate case)
)======"
R"======(        bool is_orthogonal = (n_periodic == 0) ||
)======"
R"======(                             (a_norm < tol || b_norm < tol || c_norm < tol) ||
)======"
R"======(                             ((fabs(ab_dot) < tol * a_norm * b_norm) &&
)======"
R"======(                              (fabs(ac_dot) < tol * a_norm * c_norm) &&
)======"
R"======(                              (fabs(bc_dot) < tol * b_norm * c_norm));
)======"
R"======(
)======"
R"======(        if (box_diag != nullptr) {
)======"
R"======(            box_diag[0] = a_norm;
)======"
R"======(            box_diag[1] = b_norm;
)======"
R"======(            box_diag[2] = c_norm;
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        if (inv_box_out != nullptr && !is_orthogonal) {
)======"
R"======(            double3 inv_box[3];
)======"
R"======(            invert_matrix(shared_box, inv_box);
)======"
R"======(            inv_box_out[0] = inv_box[0].x;
)======"
R"======(            inv_box_out[1] = inv_box[0].y;
)======"
R"======(            inv_box_out[2] = inv_box[0].z;
)======"
R"======(            inv_box_out[3] = inv_box[1].x;
)======"
R"======(            inv_box_out[4] = inv_box[1].y;
)======"
R"======(            inv_box_out[5] = inv_box[1].z;
)======"
R"======(            inv_box_out[6] = inv_box[2].x;
)======"
R"======(            inv_box_out[7] = inv_box[2].y;
)======"
R"======(            inv_box_out[8] = inv_box[2].z;
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        double min_dim = 1e30;
)======"
R"======(        if (is_orthogonal) {
)======"
R"======(            if (periodic[0]) {
)======"
R"======(                min_dim = a_norm;
)======"
R"======(            }
)======"
R"======(
)======"
R"======(            if (periodic[1]) {
)======"
R"======(                min_dim = fmin(min_dim, b_norm);
)======"
R"======(            }
)======"
R"======(
)======"
R"======(            if (periodic[2]) {
)======"
R"======(                min_dim = fmin(min_dim, c_norm);
)======"
R"======(            }
)======"
R"======(        } else {
)======"
R"======(            // General case
)======"
R"======(            double3 bc_cross = cross(b, c);
)======"
R"======(            double3 ac_cross = cross(a, c);
)======"
R"======(            double3 ab_cross = cross(a, b);
)======"
R"======(
)======"
R"======(            double bc_norm = norm(bc_cross);
)======"
R"======(            double ac_norm = norm(ac_cross);
)======"
R"======(            double ab_norm = norm(ab_cross);
)======"
R"======(
)======"
R"======(            double V = fabs(dot(a, bc_cross));
)======"
R"======(
)======"
R"======(            double d_a = V / bc_norm;
)======"
R"======(            double d_b = V / ac_norm;
)======"
R"======(            double d_c = V / ab_norm;
)======"
R"======(
)======"
R"======(            if (periodic[0]) {
)======"
R"======(                min_dim = d_a;
)======"
R"======(            }
)======"
R"======(
)======"
R"======(            if (periodic[1]) {
)======"
R"======(                min_dim = fmin(min_dim, d_b);
)======"
R"======(            }
)======"
R"======(
)======"
R"======(            if (periodic[2]) {
)======"
R"======(                min_dim = fmin(min_dim, d_c);
)======"
R"======(            }
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        int result = 0;
)======"
R"======(        if (cutoff * 2.0 > min_dim) {
)======"
R"======(            result |= BOX_STATUS_ERROR;
)======"
R"======(        }
)======"
R"======(        if (is_orthogonal) {
)======"
R"======(            result |= BOX_STATUS_ORTHOGONAL;
)======"
R"======(        }
)======"
R"======(        status[0] = result;
)======"
R"======(    }
)======"
R"======(}
)======"
    ;

static const char* CUDA_CELL_LIST_CODE =
R"======(// Cell list neighbor finding: bin particles into cells, then search neighboring cells.
)======"
R"======(// Particles are sorted by cell for memory coalescing. Multiple threads per particle
)======"
R"======(// cooperate on neighbor search. Output buffering reduces atomic contention.
)======"
R"======(
)======"
R"======(__device__ inline size_t atomicAdd_size_t(size_t* address, size_t val) {
)======"
R"======(    unsigned long long* address_as_ull = reinterpret_cast<unsigned long long*>(address);
)======"
R"======(    return static_cast<size_t>(atomicAdd(address_as_ull, static_cast<unsigned long long>(val)));
)======"
R"======(}
)======"
R"======(
)======"
R"======(// Vector math helpers for double3
)======"
R"======(__device__ inline double3 operator-(const double3& a, const double3& b) {
)======"
R"======(    return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
)======"
R"======(}
)======"
R"======(
)======"
R"======(__device__ inline double dot3(const double3& a, const double3& b) {
)======"
R"======(    return a.x * b.x + a.y * b.y + a.z * b.z;
)======"
R"======(}
)======"
R"======(
)======"
R"======(// Compute inv_box, n_cells, n_search from box matrix and cutoff (single thread)
)======"
R"======(__global__ void compute_cell_grid_params(
)======"
R"======(    const double* __restrict__ box,
)======"
R"======(    const bool* __restrict__ periodic,
)======"
R"======(    double cutoff,
)======"
R"======(    int max_cells,
)======"
R"======(    size_t n_points,
)======"
R"======(    int min_particles_per_cell,
)======"
R"======(    double* __restrict__ inv_box,
)======"
R"======(    int* __restrict__ n_cells,
)======"
R"======(    int* __restrict__ n_search,
)======"
R"======(    int* __restrict__ n_cells_total
)======"
R"======() {
)======"
R"======(    if (threadIdx.x != 0 || blockIdx.x != 0) {
)======"
R"======(        return;
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    // Box matrix elements
)======"
R"======(    double a = box[0], b = box[1], c = box[2];
)======"
R"======(    double d = box[3], e = box[4], f = box[5];
)======"
R"======(    double g = box[6], h = box[7], i = box[8];
)======"
R"======(
)======"
R"======(    // Determinant
)======"
R"======(    double det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
)======"
R"======(    double invdet = 1.0 / det;
)======"
R"======(
)======"
R"======(    // Inverse box matrix
)======"
R"======(    inv_box[0] = (e * i - f * h) * invdet;
)======"
R"======(    inv_box[1] = (c * h - b * i) * invdet;
)======"
R"======(    inv_box[2] = (b * f - c * e) * invdet;
)======"
R"======(    inv_box[3] = (f * g - d * i) * invdet;
)======"
R"======(    inv_box[4] = (a * i - c * g) * invdet;
)======"
R"======(    inv_box[5] = (c * d - a * f) * invdet;
)======"
R"======(    inv_box[6] = (d * h - e * g) * invdet;
)======"
R"======(    inv_box[7] = (b * g - a * h) * invdet;
)======"
R"======(    inv_box[8] = (a * e - b * d) * invdet;
)======"
R"======(
)======"
R"======(    // Box vectors
)======"
R"======(    double va[3] = {box[0], box[1], box[2]};
)======"
R"======(    double vb[3] = {box[3], box[4], box[5]};
)======"
R"======(    double vc[3] = {box[6], box[7], box[8]};
)======"
R"======(
)======"
R"======(    // Cross products for face normals
)======"
R"======(    double bc[3] = {vb[1] * vc[2] - vb[2] * vc[1], vb[2] * vc[0] - vb[0] * vc[2], vb[0] * vc[1] - vb[1] * vc[0]};
)======"
R"======(    double ca[3] = {vc[1] * va[2] - vc[2] * va[1], vc[2] * va[0] - vc[0] * va[2], vc[0] * va[1] - vc[1] * va[0]};
)======"
R"======(    double ab[3] = {va[1] * vb[2] - va[2] * vb[1], va[2] * vb[0] - va[0] * vb[2], va[0] * vb[1] - va[1] * vb[0]};
)======"
R"======(
)======"
R"======(    double bc_norm = sqrt(bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]);
)======"
R"======(    double ca_norm = sqrt(ca[0] * ca[0] + ca[1] * ca[1] + ca[2] * ca[2]);
)======"
R"======(    double ab_norm = sqrt(ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2]);
)======"
R"======(
)======"
R"======(    // Distances between opposite faces
)======"
R"======(    double dist_a = fabs(va[0] * bc[0] + va[1] * bc[1] + va[2] * bc[2]) / bc_norm;
)======"
R"======(    double dist_b = fabs(vb[0] * ca[0] + vb[1] * ca[1] + vb[2] * ca[2]) / ca_norm;
)======"
R"======(    double dist_c = fabs(vc[0] * ab[0] + vc[1] * ab[1] + vc[2] * ab[2]) / ab_norm;
)======"
R"======(    double distances[3] = {dist_a, dist_b, dist_c};
)======"
R"======(
)======"
R"======(    // Compute number of cells based on cutoff (one cell per cutoff distance)
)======"
R"======(    n_cells[0] = max(1, (int)floor(distances[0] / cutoff));
)======"
R"======(    n_cells[1] = max(1, (int)floor(distances[1] / cutoff));
)======"
R"======(    n_cells[2] = max(1, (int)floor(distances[2] / cutoff));
)======"
R"======(
)======"
R"======(    int total = n_cells[0] * n_cells[1] * n_cells[2];
)======"
R"======(
)======"
R"======(    // Compute effective max cells based on minimum particles per cell target
)======"
R"======(    // This ensures we have enough work per cell for good GPU utilization
)======"
R"======(    int max_cells_from_particles = max(1, (int)(n_points / min_particles_per_cell));
)======"
R"======(    int effective_max_cells = min(max_cells, max_cells_from_particles);
)======"
R"======(
)======"
R"======(    // Limit total cells to effective maximum
)======"
R"======(    if (total > effective_max_cells) {
)======"
R"======(        double ratio = cbrt((double)effective_max_cells / total);
)======"
R"======(        n_cells[0] = max(1, (int)(n_cells[0] * ratio));
)======"
R"======(        n_cells[1] = max(1, (int)(n_cells[1] * ratio));
)======"
R"======(        n_cells[2] = max(1, (int)(n_cells[2] * ratio));
)======"
R"======(        total = n_cells[0] * n_cells[1] * n_cells[2];
)======"
R"======(    }
)======"
R"======(    n_cells_total[0] = total;
)======"
R"======(
)======"
R"======(    // Compute search range - how many cells to search in each direction
)======"
R"======(    // When cells are larger than cutoff, we need to search more cells
)======"
R"======(    for (int dim = 0; dim < 3; dim++) {
)======"
R"======(        double cell_size = distances[dim] / n_cells[dim];
)======"
R"======(        n_search[dim] = max(1, (int)ceil(cutoff / cell_size));
)======"
R"======(        if (!periodic[dim] && n_cells[dim] == 1) {
)======"
R"======(            n_search[dim] = 0;
)======"
R"======(        }
)======"
R"======(    }
)======"
R"======(}
)======"
R"======(
)======"
R"======(// Map particles to cells via fractional coords, record periodic wrap shifts
)======"
R"======(__global__ void assign_cell_indices(
)======"
R"======(    const double* __restrict__ positions,
)======"
R"======(    const double* __restrict__ inv_box,
)======"
R"======(    const bool* __restrict__ periodic,
)======"
R"======(    const int* __restrict__ n_cells,
)======"
R"======(    size_t n_points,
)======"
R"======(    int* __restrict__ cell_indices,
)======"
R"======(    int* __restrict__ particle_shifts
)======"
R"======() {
)======"
R"======(    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
)======"
R"======(    if (i >= n_points) {
)======"
R"======(        return;
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    // Vectorized position load
)======"
R"======(    const double3* pos3 = reinterpret_cast<const double3*>(positions);
)======"
R"======(    const double3 pos = pos3[i];
)======"
R"======(
)======"
R"======(    // frac = pos @ inv_box (inv_box stored row-major: rows are inv_box[0..2], inv_box[3..5], inv_box[6..8])
)======"
R"======(    double frac[3];
)======"
R"======(    frac[0] = pos.x * inv_box[0] + pos.y * inv_box[3] + pos.z * inv_box[6];
)======"
R"======(    frac[1] = pos.x * inv_box[1] + pos.y * inv_box[4] + pos.z * inv_box[7];
)======"
R"======(    frac[2] = pos.x * inv_box[2] + pos.y * inv_box[5] + pos.z * inv_box[8];
)======"
R"======(
)======"
R"======(    int cell_idx[3];
)======"
R"======(    int shift[3];
)======"
R"======(
)======"
R"======(    for (int d = 0; d < 3; d++) {
)======"
R"======(        if (periodic[d]) {
)======"
R"======(            shift[d] = static_cast<int>(floor(frac[d]));
)======"
R"======(            frac[d] -= shift[d];
)======"
R"======(        } else {
)======"
R"======(            shift[d] = 0;
)======"
R"======(            frac[d] = fmax(0.0, fmin(frac[d], 1.0 - 1e-10));
)======"
R"======(        }
)======"
R"======(        cell_idx[d] = static_cast<int>(frac[d] * n_cells[d]);
)======"
R"======(        cell_idx[d] = max(0, min(n_cells[d] - 1, cell_idx[d]));
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    cell_indices[i] = cell_idx[0] + cell_idx[1] * n_cells[0] + cell_idx[2] * n_cells[0] * n_cells[1];
)======"
R"======(    particle_shifts[i * 3 + 0] = shift[0];
)======"
R"======(    particle_shifts[i * 3 + 1] = shift[1];
)======"
R"======(    particle_shifts[i * 3 + 2] = shift[2];
)======"
R"======(}
)======"
R"======(
)======"
R"======(// Count particles per cell (histogram)
)======"
R"======(__global__ void count_particles_per_cell(
)======"
R"======(    const int* __restrict__ cell_indices,
)======"
R"======(    size_t n_points,
)======"
R"======(    int* __restrict__ cell_counts
)======"
R"======() {
)======"
R"======(    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
)======"
R"======(    if (i >= n_points) {
)======"
R"======(        return;
)======"
R"======(    }
)======"
R"======(    atomicAdd(&cell_counts[cell_indices[i]], 1);
)======"
R"======(}
)======"
R"======(
)======"
R"======(// Exclusive prefix sum of cell counts -> cell_starts (single block, uses shared mem)
)======"
R"======(__global__ void prefix_sum_cells(
)======"
R"======(    int* __restrict__ cell_counts,
)======"
R"======(    int* __restrict__ cell_starts,
)======"
R"======(    const int* __restrict__ n_cells_total_ptr
)======"
R"======() {
)======"
R"======(    extern __shared__ int shared[];
)======"
R"======(    int tid = threadIdx.x;
)======"
R"======(    int nthreads = blockDim.x;
)======"
R"======(    int n_cells_total = n_cells_total_ptr[0];
)======"
R"======(
)======"
R"======(    // Each thread computes sum for its chunk
)======"
R"======(    int chunk_size = (n_cells_total + nthreads - 1) / nthreads;
)======"
R"======(    int start = tid * chunk_size;
)======"
R"======(    int end = min(start + chunk_size, n_cells_total);
)======"
R"======(
)======"
R"======(    // Local scan within chunk
)======"
R"======(    int local_sum = 0;
)======"
R"======(    for (int i = start; i < end; i++) {
)======"
R"======(        int val = cell_counts[i];
)======"
R"======(        cell_starts[i] = local_sum;
)======"
R"======(        local_sum += val;
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    // Store chunk totals in shared memory
)======"
R"======(    shared[tid] = local_sum;
)======"
R"======(    __syncthreads();
)======"
R"======(
)======"
R"======(    // Thread 0 computes prefix sum of chunk totals
)======"
R"======(    if (tid == 0) {
)======"
R"======(        int sum = 0;
)======"
R"======(        for (int i = 0; i < nthreads; i++) {
)======"
R"======(            int val = shared[i];
)======"
R"======(            shared[i] = sum;
)======"
R"======(            sum += val;
)======"
R"======(        }
)======"
R"======(    }
)======"
R"======(    __syncthreads();
)======"
R"======(
)======"
R"======(    // Add chunk offset to local results
)======"
R"======(    int offset = shared[tid];
)======"
R"======(    for (int i = start; i < end; i++) {
)======"
R"======(        cell_starts[i] += offset;
)======"
R"======(    }
)======"
R"======(}
)======"
R"======(
)======"
R"======(// Reorder particles by cell for coalesced access in neighbor search
)======"
R"======(__global__ void scatter_particles(
)======"
R"======(    const double* __restrict__ positions,
)======"
R"======(    const int* __restrict__ cell_indices,
)======"
R"======(    const int* __restrict__ particle_shifts,
)======"
R"======(    int* __restrict__ cell_offsets,
)======"
R"======(    size_t n_points,
)======"
R"======(    double* __restrict__ sorted_positions,
)======"
R"======(    int* __restrict__ sorted_indices,
)======"
R"======(    int* __restrict__ sorted_shifts,
)======"
R"======(    int* __restrict__ sorted_cell_indices
)======"
R"======() {
)======"
R"======(    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
)======"
R"======(    if (i >= n_points) {
)======"
R"======(        return;
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    int cell = cell_indices[i];
)======"
R"======(    int slot = atomicAdd(&cell_offsets[cell], 1);
)======"
R"======(
)======"
R"======(    // Vectorized position copy
)======"
R"======(    const double3* pos_in = reinterpret_cast<const double3*>(positions);
)======"
R"======(    double3* pos_out = reinterpret_cast<double3*>(sorted_positions);
)======"
R"======(    pos_out[slot] = pos_in[i];
)======"
R"======(
)======"
R"======(    sorted_indices[slot] = static_cast<int>(i);
)======"
R"======(    sorted_shifts[slot * 3 + 0] = particle_shifts[i * 3 + 0];
)======"
R"======(    sorted_shifts[slot * 3 + 1] = particle_shifts[i * 3 + 1];
)======"
R"======(    sorted_shifts[slot * 3 + 2] = particle_shifts[i * 3 + 2];
)======"
R"======(    sorted_cell_indices[slot] = cell;
)======"
R"======(}
)======"
R"======(
)======"
R"======(// THREADS_PER_PARTICLE threads cooperate on each particle's neighbor search
)======"
R"======(#define THREADS_PER_PARTICLE 8
)======"
R"======(// Buffer this many pairs before writing to global memory (reduces atomics)
)======"
R"======(#define MAX_BUFFERED_PAIRS 8
)======"
R"======(
)======"
R"======(// Main neighbor search kernel: each particle searches neighboring cells,
)======"
R"======(// threads within a group split the work across neighbor cells.
)======"
R"======(// Uses output buffering to batch writes and reduce atomic contention.
)======"
R"======(__global__ void find_neighbors_optimized(
)======"
R"======(    const double* __restrict__ sorted_positions,
)======"
R"======(    const int* __restrict__ sorted_indices,
)======"
R"======(    const int* __restrict__ sorted_shifts,
)======"
R"======(    const int* __restrict__ sorted_cell_indices,
)======"
R"======(    const int* __restrict__ cell_starts,
)======"
R"======(    const int* __restrict__ cell_counts,
)======"
R"======(    const double* __restrict__ box,
)======"
R"======(    const bool* __restrict__ periodic,
)======"
R"======(    const int* __restrict__ n_cells,
)======"
R"======(    const int* __restrict__ n_search,
)======"
R"======(    size_t n_points,
)======"
R"======(    double cutoff,
)======"
R"======(    bool full_list,
)======"
R"======(    size_t* length,
)======"
R"======(    size_t* pair_indices,
)======"
R"======(    int* shifts_out,
)======"
R"======(    double* distances,
)======"
R"======(    double* vectors,
)======"
R"======(    bool return_shifts,
)======"
R"======(    bool return_distances,
)======"
R"======(    bool return_vectors,
)======"
R"======(    size_t max_pairs,
)======"
R"======(    int* overflow_flag
)======"
R"======() {
)======"
R"======(    // Thread organization: 32 threads/warp, THREADS_PER_PARTICLE threads per particle
)======"
R"======(    // Example with THREADS_PER_PARTICLE=8: 4 particles per warp, each gets 8 threads
)======"
R"======(    const int warp_id = threadIdx.x / 32;
)======"
R"======(    const int lane_id = threadIdx.x % 32;
)======"
R"======(    const int particle_in_warp = lane_id / THREADS_PER_PARTICLE; // which particle in warp
)======"
R"======(    const int thread_in_group = lane_id % THREADS_PER_PARTICLE;  // thread's role within group
)======"
R"======(    const int particles_per_warp = 32 / THREADS_PER_PARTICLE;
)======"
R"======(
)======"
R"======(    const int warps_per_block = blockDim.x / 32;
)======"
R"======(    const size_t base_particle = (size_t)(blockIdx.x * warps_per_block + warp_id) * particles_per_warp;
)======"
R"======(    const size_t i = base_particle + particle_in_warp; // global particle index
)======"
R"======(
)======"
R"======(    if (i >= n_points) {
)======"
R"======(        return;
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    const double cutoff2 = cutoff * cutoff;
)======"
R"======(
)======"
R"======(    const int nc_x = __ldg(&n_cells[0]);
)======"
R"======(    const int nc_y = __ldg(&n_cells[1]);
)======"
R"======(    const int nc_z = __ldg(&n_cells[2]);
)======"
R"======(    const int nc_xy = nc_x * nc_y;
)======"
R"======(    const int ns_x = __ldg(&n_search[0]);
)======"
R"======(    const int ns_y = __ldg(&n_search[1]);
)======"
R"======(    const int ns_z = __ldg(&n_search[2]);
)======"
R"======(    const bool pbc_x = periodic[0];
)======"
R"======(    const bool pbc_y = periodic[1];
)======"
R"======(    const bool pbc_z = periodic[2];
)======"
R"======(
)======"
R"======(    // Load box matrix rows as double3
)======"
R"======(    const double3* box3 = reinterpret_cast<const double3*>(box);
)======"
R"======(    const double3 box_row0 = box3[0]; // box[0], box[1], box[2]
)======"
R"======(    const double3 box_row1 = box3[1]; // box[3], box[4], box[5]
)======"
R"======(    const double3 box_row2 = box3[2]; // box[6], box[7], box[8]
)======"
R"======(
)======"
R"======(    // Vectorized position load
)======"
R"======(    const double3* pos3 = reinterpret_cast<const double3*>(sorted_positions);
)======"
R"======(    const double3 ri = pos3[i];
)======"
R"======(    const int orig_i = __ldg(&sorted_indices[i]);
)======"
R"======(    const int shift_i_x = __ldg(&sorted_shifts[i * 3 + 0]);
)======"
R"======(    const int shift_i_y = __ldg(&sorted_shifts[i * 3 + 1]);
)======"
R"======(    const int shift_i_z = __ldg(&sorted_shifts[i * 3 + 2]);
)======"
R"======(
)======"
R"======(    const int cell_i = __ldg(&sorted_cell_indices[i]);
)======"
R"======(    const int cell_iz = cell_i / nc_xy;
)======"
R"======(    const int cell_iy = (cell_i % nc_xy) / nc_x;
)======"
R"======(    const int cell_ix = cell_i % nc_x;
)======"
R"======(
)======"
R"======(    // Per-thread output buffer to reduce atomic contention
)======"
R"======(    int buffered_count = 0;
)======"
R"======(    int buffered_j[MAX_BUFFERED_PAIRS];
)======"
R"======(    int buffered_shift[MAX_BUFFERED_PAIRS * 3];
)======"
R"======(    double buffered_dist[MAX_BUFFERED_PAIRS];
)======"
R"======(    double buffered_vec[MAX_BUFFERED_PAIRS * 3];
)======"
R"======(
)======"
R"======(    // Threads in group split neighbor cells: thread 0 does cells 0,8,16,...; thread 1 does 1,9,17,...
)======"
R"======(    int total_neighbor_cells = (2 * ns_x + 1) * (2 * ns_y + 1) * (2 * ns_z + 1);
)======"
R"======(
)======"
R"======(    for (int cell_idx = thread_in_group; cell_idx < total_neighbor_cells; cell_idx += THREADS_PER_PARTICLE) {
)======"
R"======(        // Convert linear cell_idx to 3D offset (dx,dy,dz) from particle's cell
)======"
R"======(        int temp = cell_idx;
)======"
R"======(        int dx = (temp % (2 * ns_x + 1)) - ns_x;
)======"
R"======(        temp /= (2 * ns_x + 1);
)======"
R"======(        int dy = (temp % (2 * ns_y + 1)) - ns_y;
)======"
R"======(        int dz = (temp / (2 * ns_y + 1)) - ns_z;
)======"
R"======(
)======"
R"======(        int cell_jx = cell_ix + dx;
)======"
R"======(        int cell_jy = cell_iy + dy;
)======"
R"======(        int cell_jz = cell_iz + dz;
)======"
R"======(        int cell_shift_x = 0, cell_shift_y = 0, cell_shift_z = 0;
)======"
R"======(
)======"
R"======(        // Wrap cell indices for PBC, track shift; skip out-of-bounds for non-PBC
)======"
R"======(        if (pbc_x) {
)======"
R"======(            while (cell_jx < 0) {
)======"
R"======(                cell_jx += nc_x;
)======"
R"======(                cell_shift_x -= 1;
)======"
R"======(            }
)======"
R"======(            while (cell_jx >= nc_x) {
)======"
R"======(                cell_jx -= nc_x;
)======"
R"======(                cell_shift_x += 1;
)======"
R"======(            }
)======"
R"======(        } else {
)======"
R"======(            if (cell_jx < 0 || cell_jx >= nc_x) {
)======"
R"======(                continue;
)======"
R"======(            }
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        if (pbc_y) {
)======"
R"======(            while (cell_jy < 0) {
)======"
R"======(                cell_jy += nc_y;
)======"
R"======(                cell_shift_y -= 1;
)======"
R"======(            }
)======"
R"======(            while (cell_jy >= nc_y) {
)======"
R"======(                cell_jy -= nc_y;
)======"
R"======(                cell_shift_y += 1;
)======"
R"======(            }
)======"
R"======(        } else {
)======"
R"======(            if (cell_jy < 0 || cell_jy >= nc_y) {
)======"
R"======(                continue;
)======"
R"======(            }
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        if (pbc_z) {
)======"
R"======(            while (cell_jz < 0) {
)======"
R"======(                cell_jz += nc_z;
)======"
R"======(                cell_shift_z -= 1;
)======"
R"======(            }
)======"
R"======(            while (cell_jz >= nc_z) {
)======"
R"======(                cell_jz -= nc_z;
)======"
R"======(                cell_shift_z += 1;
)======"
R"======(            }
)======"
R"======(        } else {
)======"
R"======(            if (cell_jz < 0 || cell_jz >= nc_z) {
)======"
R"======(                continue;
)======"
R"======(            }
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        const int cell_j = cell_jx + cell_jy * nc_x + cell_jz * nc_xy;
)======"
R"======(        const int start_j = __ldg(&cell_starts[cell_j]);
)======"
R"======(        const int count_j = __ldg(&cell_counts[cell_j]);
)======"
R"======(
)======"
R"======(        for (int k = start_j; k < start_j + count_j; k++) {
)======"
R"======(            const int orig_j = __ldg(&sorted_indices[k]);
)======"
R"======(
)======"
R"======(            const int shift_j_x = __ldg(&sorted_shifts[k * 3 + 0]);
)======"
R"======(            const int shift_j_y = __ldg(&sorted_shifts[k * 3 + 1]);
)======"
R"======(            const int shift_j_z = __ldg(&sorted_shifts[k * 3 + 2]);
)======"
R"======(
)======"
R"======(            const int total_shift_x = shift_i_x - shift_j_x + cell_shift_x;
)======"
R"======(            const int total_shift_y = shift_i_y - shift_j_y + cell_shift_y;
)======"
R"======(            const int total_shift_z = shift_i_z - shift_j_z + cell_shift_z;
)======"
R"======(
)======"
R"======(            const bool shift_is_zero = (total_shift_x == 0 && total_shift_y == 0 && total_shift_z == 0);
)======"
R"======(
)======"
R"======(            if (orig_i == orig_j && shift_is_zero) {
)======"
R"======(                continue;
)======"
R"======(            }
)======"
R"======(
)======"
R"======(            // Half-list: only keep pair once (i<j, or i==j with positive shift direction)
)======"
R"======(            if (!full_list) {
)======"
R"======(                if (orig_i > orig_j) {
)======"
R"======(                    continue;
)======"
R"======(                }
)======"
R"======(                if (orig_i == orig_j) {
)======"
R"======(                    // For self-images, use lexicographic ordering on shift vector
)======"
R"======(                    int shift_sum = total_shift_x + total_shift_y + total_shift_z;
)======"
R"======(                    if (shift_sum < 0) {
)======"
R"======(                        continue;
)======"
R"======(                    }
)======"
R"======(                    if (shift_sum == 0) {
)======"
R"======(                        if (total_shift_z < 0 || (total_shift_z == 0 && total_shift_y < 0)) {
)======"
R"======(                            continue;
)======"
R"======(                        }
)======"
R"======(                    }
)======"
R"======(                }
)======"
R"======(            }
)======"
R"======(
)======"
R"======(            // Vectorized position load for particle j
)======"
R"======(            const double3 rj = pos3[k];
)======"
R"======(
)======"
R"======(            // Convert integer shift to Cartesian displacement: shift @ box
)======"
R"======(            const double3 shift_cart = make_double3(
)======"
R"======(                total_shift_x * box_row0.x + total_shift_y * box_row1.x + total_shift_z * box_row2.x,
)======"
R"======(                total_shift_x * box_row0.y + total_shift_y * box_row1.y + total_shift_z * box_row2.y,
)======"
R"======(                total_shift_x * box_row0.z + total_shift_y * box_row1.z + total_shift_z * box_row2.z
)======"
R"======(            );
)======"
R"======(
)======"
R"======(            // Vector from i to j (accounting for PBC)
)======"
R"======(            const double3 vec = make_double3(
)======"
R"======(                rj.x - ri.x + shift_cart.x,
)======"
R"======(                rj.y - ri.y + shift_cart.y,
)======"
R"======(                rj.z - ri.z + shift_cart.z
)======"
R"======(            );
)======"
R"======(
)======"
R"======(            const double dist2 = dot3(vec, vec);
)======"
R"======(
)======"
R"======(            if (dist2 < cutoff2 && dist2 > 0.0) {
)======"
R"======(                buffered_j[buffered_count] = orig_j;
)======"
R"======(                buffered_shift[buffered_count * 3 + 0] = total_shift_x;
)======"
R"======(                buffered_shift[buffered_count * 3 + 1] = total_shift_y;
)======"
R"======(                buffered_shift[buffered_count * 3 + 2] = total_shift_z;
)======"
R"======(                buffered_dist[buffered_count] = sqrt(dist2);
)======"
R"======(                buffered_vec[buffered_count * 3 + 0] = vec.x;
)======"
R"======(                buffered_vec[buffered_count * 3 + 1] = vec.y;
)======"
R"======(                buffered_vec[buffered_count * 3 + 2] = vec.z;
)======"
R"======(                buffered_count++;
)======"
R"======(
)======"
R"======(                // Flush buffer when full
)======"
R"======(                if (buffered_count >= MAX_BUFFERED_PAIRS) {
)======"
R"======(                    size_t base_idx = atomicAdd_size_t(length, buffered_count);
)======"
R"======(
)======"
R"======(                    // Check if we are about to exceed max_pairs
)======"
R"======(                    if (base_idx + buffered_count > max_pairs) {
)======"
R"======(                        atomicExch(overflow_flag, 1);
)======"
R"======(                        return;
)======"
R"======(                    }
)======"
R"======(
)======"
R"======(                    for (int b = 0; b < buffered_count; b++) {
)======"
R"======(                        pair_indices[(base_idx + b) * 2] = orig_i;
)======"
R"======(                        pair_indices[(base_idx + b) * 2 + 1] = buffered_j[b];
)======"
R"======(                        if (return_shifts) {
)======"
R"======(                            shifts_out[(base_idx + b) * 3 + 0] = buffered_shift[b * 3 + 0];
)======"
R"======(                            shifts_out[(base_idx + b) * 3 + 1] = buffered_shift[b * 3 + 1];
)======"
R"======(                            shifts_out[(base_idx + b) * 3 + 2] = buffered_shift[b * 3 + 2];
)======"
R"======(                        }
)======"
R"======(                        if (return_distances) {
)======"
R"======(                            distances[base_idx + b] = buffered_dist[b];
)======"
R"======(                        }
)======"
R"======(                        if (return_vectors) {
)======"
R"======(                            vectors[(base_idx + b) * 3 + 0] = buffered_vec[b * 3 + 0];
)======"
R"======(                            vectors[(base_idx + b) * 3 + 1] = buffered_vec[b * 3 + 1];
)======"
R"======(                            vectors[(base_idx + b) * 3 + 2] = buffered_vec[b * 3 + 2];
)======"
R"======(                        }
)======"
R"======(                    }
)======"
R"======(                    buffered_count = 0;
)======"
R"======(                }
)======"
R"======(            }
)======"
R"======(        }
)======"
R"======(    }
)======"
R"======(
)======"
R"======(    // Flush remaining buffered pairs
)======"
R"======(    if (buffered_count > 0) {
)======"
R"======(        size_t base_idx = atomicAdd_size_t(length, buffered_count);
)======"
R"======(
)======"
R"======(        // Check if we are about to exceed max_pairs
)======"
R"======(        if (base_idx + buffered_count > max_pairs) {
)======"
R"======(            atomicExch(overflow_flag, 1);
)======"
R"======(            return;
)======"
R"======(        }
)======"
R"======(
)======"
R"======(        for (int b = 0; b < buffered_count; b++) {
)======"
R"======(            pair_indices[(base_idx + b) * 2] = orig_i;
)======"
R"======(            pair_indices[(base_idx + b) * 2 + 1] = buffered_j[b];
)======"
R"======(            if (return_shifts) {
)======"
R"======(                shifts_out[(base_idx + b) * 3 + 0] = buffered_shift[b * 3 + 0];
)======"
R"======(                shifts_out[(base_idx + b) * 3 + 1] = buffered_shift[b * 3 + 1];
)======"
R"======(                shifts_out[(base_idx + b) * 3 + 2] = buffered_shift[b * 3 + 2];
)======"
R"======(            }
)======"
R"======(
)======"
R"======(            if (return_distances) {
)======"
R"======(                distances[base_idx + b] = buffered_dist[b];
)======"
R"======(            }
)======"
R"======(
)======"
R"======(            if (return_vectors) {
)======"
R"======(                vectors[(base_idx + b) * 3 + 0] = buffered_vec[b * 3 + 0];
)======"
R"======(                vectors[(base_idx + b) * 3 + 1] = buffered_vec[b * 3 + 1];
)======"
R"======(                vectors[(base_idx + b) * 3 + 2] = buffered_vec[b * 3 + 2];
)======"
R"======(            }
)======"
R"======(        }
)======"
R"======(    }
)======"
R"======(}
)======"
    ;

// Maximum number of cells (limited by single-block prefix sum)
static constexpr size_t MAX_CELLS = 8192;
// Minimum particles per cell target for good GPU utilization
// Higher values = fewer cells = more work per block but larger search range
static constexpr size_t MIN_PARTICLES_PER_CELL = 128;

// Helper functions for CPU-side vector math
static inline double cpu_dot3(const double* a, const double* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline double cpu_norm3(const double* v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static inline void cpu_cross3(const double* a, const double* b, double* result) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

static inline void cpu_invert_matrix(const double* m, double* inv) {
    double det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);

    double inv_det = 1.0 / det;

    inv[0] = (m[4] * m[8] - m[5] * m[7]) * inv_det;
    inv[1] = (m[2] * m[7] - m[1] * m[8]) * inv_det;
    inv[2] = (m[1] * m[5] - m[2] * m[4]) * inv_det;
    inv[3] = (m[5] * m[6] - m[3] * m[8]) * inv_det;
    inv[4] = (m[0] * m[8] - m[2] * m[6]) * inv_det;
    inv[5] = (m[2] * m[3] - m[0] * m[5]) * inv_det;
    inv[6] = (m[3] * m[7] - m[4] * m[6]) * inv_det;
    inv[7] = (m[1] * m[6] - m[0] * m[7]) * inv_det;
    inv[8] = (m[0] * m[4] - m[1] * m[3]) * inv_det;
}

/// CPU-side box check that avoids GPU kernel launch overhead
/// Returns: {is_valid, is_orthogonal}
/// Also fills box_diag_out[3] and inv_box_out[9] if provided
static std::pair<bool, bool> cpu_box_check(
    const double h_box[9],
    const bool h_periodic[3],
    double cutoff,
    double* box_diag_out, // [3] output, can be nullptr
    double* inv_box_out   // [9] output, can be nullptr
) {
    const double* a = &h_box[0];
    const double* b = &h_box[3];
    const double* c = &h_box[6];

    double a_norm = cpu_norm3(a);
    double b_norm = cpu_norm3(b);
    double c_norm = cpu_norm3(c);

    // Count periodic directions
    size_t n_periodic = 0;
    if (h_periodic[0]) {
        n_periodic++;
    }
    if (h_periodic[1]) {
        n_periodic++;
    }
    if (h_periodic[2]) {
        n_periodic++;
    }

    double ab_dot = cpu_dot3(a, b);
    double ac_dot = cpu_dot3(a, c);
    double bc_dot = cpu_dot3(b, c);

    double tol = 1e-6;
    // Treat fully non-periodic systems as orthogonal
    // Also treat systems with zero-norm vectors as orthogonal (degenerate case)
    bool is_orthogonal = (n_periodic == 0) ||
                         (a_norm < tol || b_norm < tol || c_norm < tol) ||
                         ((std::fabs(ab_dot) < tol * a_norm * b_norm) &&
                          (std::fabs(ac_dot) < tol * a_norm * c_norm) &&
                          (std::fabs(bc_dot) < tol * b_norm * c_norm));

    // Output box diagonal (lengths)
    if (box_diag_out != nullptr) {
        box_diag_out[0] = a_norm;
        box_diag_out[1] = b_norm;
        box_diag_out[2] = c_norm;
    }

    // Compute and output inverse box (needed for general PBC)
    if ((inv_box_out != nullptr) && !is_orthogonal) {
        cpu_invert_matrix(h_box, inv_box_out);
    }

    // Compute minimum dimension for cutoff check
    double min_dim = 1e30;
    if (is_orthogonal) {
        if (h_periodic[0]) {
            min_dim = a_norm;
        }
        if (h_periodic[1]) {
            min_dim = std::fmin(min_dim, b_norm);
        }
        if (h_periodic[2]) {
            min_dim = std::fmin(min_dim, c_norm);
        }
    } else {
        // General case: compute perpendicular distances
        double bc_cross[3];
        double ac_cross[3];
        double ab_cross[3];
        cpu_cross3(b, c, bc_cross);
        cpu_cross3(a, c, ac_cross);
        cpu_cross3(a, b, ab_cross);

        double bc_norm = cpu_norm3(bc_cross);
        double ac_norm = cpu_norm3(ac_cross);
        double ab_norm = cpu_norm3(ab_cross);

        double V = std::fabs(cpu_dot3(a, bc_cross));

        double d_a = V / bc_norm;
        double d_b = V / ac_norm;
        double d_c = V / ab_norm;

        if (h_periodic[0]) {
            min_dim = d_a;
        }
        if (h_periodic[1]) {
            min_dim = std::fmin(min_dim, d_b);
        }
        if (h_periodic[2]) {
            min_dim = std::fmin(min_dim, d_c);
        }
    }

    bool is_valid = (cutoff * 2.0 <= min_dim);
    return {is_valid, is_orthogonal};
}

static std::optional<cudaPointerAttributes> getPtrAttributes(const void* ptr) {
    if (ptr == nullptr) {
        return std::nullopt;
    }

    try {
        cudaPointerAttributes attr;
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaPointerGetAttributes(&attr, ptr));
        return attr;
    } catch (const std::runtime_error& e) {
        return std::nullopt;
    }
}

static bool is_device_ptr(const std::optional<cudaPointerAttributes>& maybe_attr, const char* name) {
    if (maybe_attr) {
        const cudaPointerAttributes& attr = *maybe_attr;
        return (attr.type == cudaMemoryTypeDevice);
    } else {
        throw std::runtime_error(
            "failed to resolve attributes for pointer: " + std::string(name)
        );
    }
}

static int32_t get_device_id(const void* ptr) {
    if (ptr == nullptr) {
        return -1;
    }

    auto maybe_attr = getPtrAttributes(ptr);
    if (maybe_attr) {
        const cudaPointerAttributes& attr = *maybe_attr;
        if (attr.type != cudaMemoryTypeDevice) {
            return -1;
        }
        return attr.device;
    }
    return -1;
}

static void free_cell_list_buffers(CellListBuffers& cl) {
    CUDART_INSTANCE.cudaFree(cl.cell_indices);
    CUDART_INSTANCE.cudaFree(cl.particle_shifts);
    CUDART_INSTANCE.cudaFree(cl.cell_counts);
    CUDART_INSTANCE.cudaFree(cl.cell_starts);
    CUDART_INSTANCE.cudaFree(cl.cell_offsets);
    CUDART_INSTANCE.cudaFree(cl.sorted_positions);
    CUDART_INSTANCE.cudaFree(cl.sorted_indices);
    CUDART_INSTANCE.cudaFree(cl.sorted_shifts);
    CUDART_INSTANCE.cudaFree(cl.sorted_cell_indices);
    CUDART_INSTANCE.cudaFree(cl.inv_box);
    CUDART_INSTANCE.cudaFree(cl.n_cells);
    CUDART_INSTANCE.cudaFree(cl.n_search);
    CUDART_INSTANCE.cudaFree(cl.n_cells_total);

    cl = CellListBuffers();
}

CudaNeighborListExtras::~CudaNeighborListExtras() {
    if (this->length_ptr != nullptr) {
        CUDART_INSTANCE.cudaFree(this->length_ptr);
    }
    if (this->cell_check_ptr != nullptr) {
        CUDART_INSTANCE.cudaFree(this->cell_check_ptr);
    }
    if (this->overflow_flag != nullptr) {
        CUDART_INSTANCE.cudaFree(this->overflow_flag);
    }
    if (this->box_diag != nullptr) {
        CUDART_INSTANCE.cudaFree(this->box_diag);
    }
    if (this->inv_box_brute != nullptr) {
        CUDART_INSTANCE.cudaFree(this->inv_box_brute);
    }
    free_cell_list_buffers(this->cell_list);
}

PLMD::metatomic::vesin::cuda::CudaNeighborListExtras*
PLMD::metatomic::vesin::cuda::get_cuda_extras(VesinNeighborList* neighbors) {
    if (neighbors->opaque == nullptr) {
        neighbors->opaque = new PLMD::metatomic::vesin::cuda::CudaNeighborListExtras();
        auto* test = static_cast<PLMD::metatomic::vesin::cuda::CudaNeighborListExtras*>(neighbors->opaque);
    }
    return static_cast<PLMD::metatomic::vesin::cuda::CudaNeighborListExtras*>(neighbors->opaque);
}

static void reset(VesinNeighborList& neighbors) {
    auto* extras = PLMD::metatomic::vesin::cuda::get_cuda_extras(&neighbors);

    if ((neighbors.pairs != nullptr) && is_device_ptr(getPtrAttributes(neighbors.pairs), "pairs")) {
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaFree(neighbors.pairs));
    }
    if ((neighbors.shifts != nullptr) && is_device_ptr(getPtrAttributes(neighbors.shifts), "shifts")) {
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaFree(neighbors.shifts));
    }
    if ((neighbors.distances != nullptr) && is_device_ptr(getPtrAttributes(neighbors.distances), "distances")) {
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaFree(neighbors.distances));
    }
    if ((neighbors.vectors != nullptr) && is_device_ptr(getPtrAttributes(neighbors.vectors), "vectors")) {
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaFree(neighbors.vectors));
    }

    neighbors.pairs = nullptr;
    neighbors.shifts = nullptr;
    neighbors.distances = nullptr;
    neighbors.vectors = nullptr;
    extras->length_ptr = nullptr;

    // Free pinned memory if allocated
    if (extras->pinned_length_ptr != nullptr) {
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaFreeHost(extras->pinned_length_ptr));
        extras->pinned_length_ptr = nullptr;
    }

    // Free brute force buffers
    CUDART_SAFE_CALL(CUDART_INSTANCE.cudaFree(extras->box_diag));
    extras->box_diag = nullptr;

    CUDART_SAFE_CALL(CUDART_INSTANCE.cudaFree(extras->inv_box_brute));
    extras->inv_box_brute = nullptr;

    CUDART_SAFE_CALL(CUDART_INSTANCE.cudaFree(extras->overflow_flag));
    extras->overflow_flag = nullptr;

    free_cell_list_buffers(extras->cell_list);

    *extras = CudaNeighborListExtras();
}

void PLMD::metatomic::vesin::cuda::free_neighbors(VesinNeighborList& neighbors) {
    assert(neighbors.device.type == VesinCUDA);

    int32_t curr_device = -1;
    int32_t device_id = -1;

    if (neighbors.pairs != nullptr) {
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaGetDevice(&curr_device));
        device_id = get_device_id(neighbors.pairs);

        if ((device_id != -1) && curr_device != device_id) {
            CUDART_SAFE_CALL(CUDART_INSTANCE.cudaSetDevice(device_id));
        }
    }

    reset(neighbors);

    if ((device_id != -1) && curr_device != device_id) {
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaSetDevice(curr_device));
    }

    delete static_cast<PLMD::metatomic::vesin::cuda::CudaNeighborListExtras*>(neighbors.opaque);
    neighbors.opaque = nullptr;
}

void checkCuda() {
    std::string cuda_libname;
    std::string cudart_libname;
    std::string nvrtc_libname;
    std::string suggestion;
#if defined(__linux__)
    cuda_libname = "libcuda.so";
    cudart_libname = "libcudart.so(.*)";
    nvrtc_libname = "libnvrtc.so(.*)";
    suggestion = ("Try appending the directory containing this library to "
                  "your $LD_LIBRARY_PATH environment variable.");

#elif defined(_WIN32)
    cuda_libname = "nvcuda.dll";
    cudart_libname = "cudart64_*.dll";
    nvrtc_libname = "nvrtc64_*.dll";
    suggestion = ("Try adding the directory containing this library to your "
                  "system PATH, or making sure that CUDA_PATH is properly set "
                  "to your CUDA installation directory.");
#else
    cuda_libname = "cuda";
    cudart_libname = "cudart";
    nvrtc_libname = "nvrtc";
    suggestion = "Unsupported platform: unable to load CUDA libraries.";
#endif
    if (!CUDA_DRIVER_INSTANCE.loaded()) {
        throw std::runtime_error(
            "Failed to load " + cuda_libname + ". " + suggestion
        );
    }

    if (!CUDART_INSTANCE.loaded()) {
        throw std::runtime_error(
            "Failed to load " + cudart_libname + ". " + suggestion
        );
    }

    if (!NVRTC_INSTANCE.loaded()) {
        throw std::runtime_error(
            "Failed to load " + nvrtc_libname + ". " + suggestion
        );
    }
}

// Ensure cell list buffers are allocated with sufficient capacity
static void ensure_cell_list_buffers(
    CellListBuffers& cl,
    size_t n_points,
    size_t n_cells_total
) {
    bool need_realloc_points = (cl.max_points < n_points);
    bool need_realloc_cells = (cl.max_cells < n_cells_total);

    if (need_realloc_points) {
        // Free old point-related buffers
        CUDART_INSTANCE.cudaFree(cl.cell_indices);
        CUDART_INSTANCE.cudaFree(cl.particle_shifts);
        CUDART_INSTANCE.cudaFree(cl.sorted_positions);
        CUDART_INSTANCE.cudaFree(cl.sorted_indices);
        CUDART_INSTANCE.cudaFree(cl.sorted_shifts);
        CUDART_INSTANCE.cudaFree(cl.sorted_cell_indices);

        auto new_max = static_cast<size_t>(1.2 * static_cast<double>(n_points));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.cell_indices, sizeof(int32_t) * new_max));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.particle_shifts, sizeof(int32_t) * new_max * 3));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.sorted_positions, sizeof(double) * new_max * 3));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.sorted_indices, sizeof(int32_t) * new_max));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.sorted_shifts, sizeof(int32_t) * new_max * 3));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.sorted_cell_indices, sizeof(int32_t) * new_max));
        cl.max_points = new_max;
    }

    if (need_realloc_cells) {
        // Free old cell-related buffers
        CUDART_INSTANCE.cudaFree(cl.cell_counts);
        CUDART_INSTANCE.cudaFree(cl.cell_starts);
        CUDART_INSTANCE.cudaFree(cl.cell_offsets);

        auto new_max = static_cast<size_t>(1.2 * static_cast<double>(n_cells_total));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.cell_counts, sizeof(int32_t) * new_max));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.cell_starts, sizeof(int32_t) * new_max));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.cell_offsets, sizeof(int32_t) * new_max));
        cl.max_cells = new_max;
    }

    // Allocate cell grid parameter buffers (fixed size, only once)
    if (cl.inv_box == nullptr) {
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.inv_box, sizeof(double) * 9));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.n_cells, sizeof(int32_t) * 3));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.n_search, sizeof(int32_t) * 3));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&cl.n_cells_total, sizeof(int32_t)));
    }
}

void PLMD::metatomic::vesin::cuda::neighbors(
    const double (*points)[3],
    size_t n_points,
    const double box[3][3],
    const bool periodic[3],
    VesinOptions options,
    VesinNeighborList& neighbors
) {
    assert(neighbors.device.type == VesinCUDA);
    if (options.sorted) {
        throw std::runtime_error("CUDA implemented does not support sorted output yet");
    }

    // Check if CUDA is available
    checkCuda();

    // check that all pointers are are device pointers
    if (!is_device_ptr(getPtrAttributes(points), "points")) {
        throw std::runtime_error("`points` pointer is not allocated on a CUDA device");
    }

    if (!is_device_ptr(getPtrAttributes(box), "box")) {
        throw std::runtime_error("`box` pointer is not allocated on a CUDA device");
    }

    if (!is_device_ptr(getPtrAttributes(periodic), "periodic")) {
        throw std::runtime_error("`periodic` pointer is not allocated on a CUDA device");
    }

    auto device_id = get_device_id(points);
    if (device_id != get_device_id(box)) {
        throw std::runtime_error("`points` and `box` do not exist on the same device");
    }

    if (device_id != get_device_id(periodic)) {
        throw std::runtime_error("`points` and `periodic` do not exist on the same device");
    }

    if (device_id != neighbors.device.device_id) {
        throw std::runtime_error("`points`, `box` and `periodic` device differs from input neighbors device_id");
    }

    auto* extras = PLMD::metatomic::vesin::cuda::get_cuda_extras(&neighbors);
    size_t max_pairs_per_point = VESIN_DEFAULT_CUDA_MAX_PAIRS_PER_POINT;

    if (extras->allocated_device_id != device_id) {
        // first switch to previous device
        if (extras->allocated_device_id >= 0) {
            CUDART_SAFE_CALL(CUDART_INSTANCE.cudaSetDevice(extras->allocated_device_id));
        }
        // free any existing allocations
        reset(neighbors);
        // switch back to current device
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaSetDevice(device_id));
        extras->allocated_device_id = device_id;
    }

    if (extras->capacity >= n_points && (extras->length_ptr != nullptr)) {
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemset(extras->length_ptr, 0, sizeof(size_t)));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemset(extras->cell_check_ptr, 0, sizeof(int32_t)));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemset(extras->overflow_flag, 0, sizeof(int32_t)));
    } else {
        auto saved_device = extras->allocated_device_id;
        reset(neighbors);
        extras->allocated_device_id = saved_device;

        auto* env_max_pairs = std::getenv("VESIN_CUDA_MAX_PAIRS_PER_POINT");
        if (env_max_pairs != nullptr) {
            auto length = std::strlen(env_max_pairs);
            char* end = nullptr;
            errno = 0;
            auto parsed_max_pairs_per_point = std::strtoll(env_max_pairs, &end, 10);
            if (errno != 0 || end != env_max_pairs + length || parsed_max_pairs_per_point <= 0) {
                throw std::runtime_error(
                    "Invalid value for VESIN_CUDA_MAX_PAIRS_PER_POINT: '" +
                    std::string(env_max_pairs) + "'"
                );
            }
            max_pairs_per_point = static_cast<size_t>(parsed_max_pairs_per_point);
        }

        auto max_pairs = static_cast<size_t>(1.2 * static_cast<double>(n_points * max_pairs_per_point));
        extras->max_pairs = max_pairs;

        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&neighbors.pairs, sizeof(size_t) * max_pairs * 2));

        if (options.return_shifts) {
            CUDART_SAFE_CALL(
                CUDART_INSTANCE.cudaMalloc((void**)&neighbors.shifts, sizeof(int32_t) * max_pairs * 3)
            );
        }

        if (options.return_distances) {
            CUDART_SAFE_CALL(
                CUDART_INSTANCE.cudaMalloc((void**)&neighbors.distances, sizeof(double) * max_pairs)
            );
        }

        if (options.return_vectors) {
            CUDART_SAFE_CALL(
                CUDART_INSTANCE.cudaMalloc((void**)&neighbors.vectors, sizeof(double) * max_pairs * 3)
            );
        }

        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&extras->length_ptr, sizeof(size_t)));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemset(extras->length_ptr, 0, sizeof(size_t)));

        // Pinned host memory for async D2H copy
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaHostAlloc(
            (void**)&extras->pinned_length_ptr,
            sizeof(size_t),
            cudaHostAllocDefault
        ));

        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&extras->cell_check_ptr, sizeof(int32_t)));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemset(extras->cell_check_ptr, 0, sizeof(int32_t)));

        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&extras->overflow_flag, sizeof(int32_t)));
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemset(extras->overflow_flag, 0, sizeof(int32_t)));

        extras->capacity = static_cast<size_t>(1.2 * static_cast<double>(n_points));
    }

    const auto* d_positions = reinterpret_cast<const double*>(points);
    const auto* d_box = reinterpret_cast<const double*>(box);
    const auto* d_periodic = periodic;

    auto* d_pair_indices = reinterpret_cast<size_t*>(neighbors.pairs);
    auto* d_shifts = reinterpret_cast<int32_t*>(neighbors.shifts);
    auto* d_distances = neighbors.distances;
    auto* d_vectors = reinterpret_cast<double*>(neighbors.vectors);
    auto* d_pair_counter = extras->length_ptr;
    auto* d_cell_check = extras->cell_check_ptr;
    auto* d_overflow_flag = extras->overflow_flag;
    size_t max_pairs = extras->max_pairs;

    auto& factory = KernelFactory::instance(device_id);

    if (extras->box_diag == nullptr) {
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&extras->box_diag, sizeof(double) * 3));
    }
    if (extras->inv_box_brute == nullptr) {
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMalloc((void**)&extras->inv_box_brute, sizeof(double) * 9));
    }

    auto* box_check_kernel = factory.create(
        "mic_box_check",
        CUDA_BRUTEFORCE_CODE,
        "cuda_bruteforce.cu",
        {"-std=c++17"}
    );

    double* d_box_diag = extras->box_diag;
    double* d_inv_box_brute = extras->inv_box_brute;
    std::vector<void*> box_check_args = {
        static_cast<void*>(&d_box),
        static_cast<void*>(&d_periodic),
        static_cast<void*>(&options.cutoff),
        static_cast<void*>(&d_cell_check),
        static_cast<void*>(&d_box_diag),
        static_cast<void*>(&d_inv_box_brute),
    };

    box_check_kernel->launch(dim3(1), dim3(32), 0, nullptr, box_check_args, false);

    int32_t h_cell_check = 1;
    CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemcpy(&h_cell_check, d_cell_check, sizeof(int32_t), cudaMemcpyDeviceToHost));

    bool box_check_error = (h_cell_check & 1) != 0;
    bool is_orthogonal = (h_cell_check & 2) != 0;

    // Get box dimensions for auto algorithm selection
    double h_box_diag[3];
    CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemcpy(h_box_diag, d_box_diag, sizeof(double) * 3, cudaMemcpyDeviceToHost));
    double min_box_dim = std::min({h_box_diag[0], h_box_diag[1], h_box_diag[2]});
    bool cutoff_requires_cell_list = options.cutoff > min_box_dim / 2.0;

    bool use_cell_list;
    switch (options.algorithm) {
    case VesinBruteForce:
        if (box_check_error) {
            throw std::runtime_error("Invalid cutoff: too large for box dimensions");
        }
        use_cell_list = false;
        break;
    case VesinCellList:
        use_cell_list = true;
        break;
    case VesinAutoAlgorithm:
    default:
        // Use cell list if cutoff > half box size, or for large/non-orthogonal systems
        use_cell_list = cutoff_requires_cell_list || !is_orthogonal || n_points >= 5000;
        break;
    }

    if (use_cell_list) {
        NVTX_PUSH("cell_list_total");

        NVTX_PUSH("ensure_buffers");
        ensure_cell_list_buffers(extras->cell_list, n_points, MAX_CELLS);
        NVTX_POP();
        auto& cl = extras->cell_list;

        int32_t max_cells_int = static_cast<int32_t>(MAX_CELLS);
        int32_t min_particles_per_cell = MIN_PARTICLES_PER_CELL;

        size_t THREADS_PER_BLOCK = 256;
        size_t num_blocks_points = (n_points + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

        NVTX_PUSH("kernel0_grid_params");
        auto* grid_kernel = factory.create(
            "compute_cell_grid_params",
            CUDA_CELL_LIST_CODE,
            "cuda_cell_list.cu",
            {"-std=c++17"}
        );
        std::vector<void*> grid_args = {
            static_cast<void*>(&d_box),
            static_cast<void*>(&d_periodic),
            static_cast<void*>(&options.cutoff),
            static_cast<void*>(&max_cells_int),
            static_cast<void*>(&n_points),
            static_cast<void*>(&min_particles_per_cell),
            static_cast<void*>(&cl.inv_box),
            static_cast<void*>(&cl.n_cells),
            static_cast<void*>(&cl.n_search),
            static_cast<void*>(&cl.n_cells_total),
        };
        grid_kernel->launch(dim3(1), dim3(1), 0, nullptr, grid_args, false);
        NVTX_POP();

        NVTX_PUSH("memset_cell_counts");
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemset(cl.cell_counts, 0, sizeof(int32_t) * MAX_CELLS));
        NVTX_POP();

        NVTX_PUSH("memset_cell_starts");
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemset(cl.cell_starts, 0, sizeof(int32_t) * MAX_CELLS));
        NVTX_POP();

        NVTX_PUSH("kernel1_assign_cells");
        auto* assign_kernel = factory.create(
            "assign_cell_indices",
            CUDA_CELL_LIST_CODE,
            "cuda_cell_list.cu",
            {"-std=c++17"}
        );
        std::vector<void*> assign_args = {
            static_cast<void*>(&d_positions),
            static_cast<void*>(&cl.inv_box),
            static_cast<void*>(&d_periodic),
            static_cast<void*>(&cl.n_cells),
            static_cast<void*>(&n_points),
            static_cast<void*>(&cl.cell_indices),
            static_cast<void*>(&cl.particle_shifts),
        };
        assign_kernel->launch(
            dim3(num_blocks_points), dim3(THREADS_PER_BLOCK), 0, nullptr, assign_args, false
        );
        NVTX_POP();

        NVTX_PUSH("kernel2_count_particles");
        auto* count_kernel = factory.create(
            "count_particles_per_cell",
            CUDA_CELL_LIST_CODE,
            "cuda_cell_list.cu",
            {"-std=c++17"}
        );
        std::vector<void*> count_args = {
            static_cast<void*>(&cl.cell_indices),
            static_cast<void*>(&n_points),
            static_cast<void*>(&cl.cell_counts),
        };
        count_kernel->launch(
            dim3(num_blocks_points), dim3(THREADS_PER_BLOCK), 0, nullptr, count_args, false
        );
        NVTX_POP();

        NVTX_PUSH("kernel3_prefix_sum");
        auto* prefix_kernel = factory.create(
            "prefix_sum_cells",
            CUDA_CELL_LIST_CODE,
            "cuda_cell_list.cu",
            {"-std=c++17"}
        );
        std::vector<void*> prefix_args = {
            static_cast<void*>(&cl.cell_counts),
            static_cast<void*>(&cl.cell_starts),
            static_cast<void*>(&cl.n_cells_total),
        };
        size_t prefix_threads = 256;
        size_t shared_mem = sizeof(int32_t) * prefix_threads;
        prefix_kernel->launch(
            dim3(1), dim3(prefix_threads), shared_mem, nullptr, prefix_args, false
        );
        NVTX_POP();

        NVTX_PUSH("memcpy_cell_offsets");
        CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemcpy(
            cl.cell_offsets, cl.cell_starts, sizeof(int32_t) * MAX_CELLS, cudaMemcpyDeviceToDevice
        ));
        NVTX_POP();

        NVTX_PUSH("kernel4_scatter");
        auto* scatter_kernel = factory.create(
            "scatter_particles",
            CUDA_CELL_LIST_CODE,
            "cuda_cell_list.cu",
            {"-std=c++17"}
        );
        std::vector<void*> scatter_args = {
            static_cast<void*>(&d_positions),
            static_cast<void*>(&cl.cell_indices),
            static_cast<void*>(&cl.particle_shifts),
            static_cast<void*>(&cl.cell_offsets),
            static_cast<void*>(&n_points),
            static_cast<void*>(&cl.sorted_positions),
            static_cast<void*>(&cl.sorted_indices),
            static_cast<void*>(&cl.sorted_shifts),
            static_cast<void*>(&cl.sorted_cell_indices),
        };
        scatter_kernel->launch(
            dim3(num_blocks_points), dim3(THREADS_PER_BLOCK), 0, nullptr, scatter_args, false
        );
        NVTX_POP();

        NVTX_PUSH("kernel5_find_neighbors");
        auto* find_kernel = factory.create(
            "find_neighbors_optimized",
            CUDA_CELL_LIST_CODE,
            "cuda_cell_list.cu",
            {"-std=c++17"}
        );
        std::vector<void*> find_args = {
            static_cast<void*>(&cl.sorted_positions),
            static_cast<void*>(&cl.sorted_indices),
            static_cast<void*>(&cl.sorted_shifts),
            static_cast<void*>(&cl.sorted_cell_indices),
            static_cast<void*>(&cl.cell_starts),
            static_cast<void*>(&cl.cell_counts),
            static_cast<void*>(&d_box),
            static_cast<void*>(&d_periodic),
            static_cast<void*>(&cl.n_cells),
            static_cast<void*>(&cl.n_search),
            static_cast<void*>(&n_points),
            static_cast<void*>(&options.cutoff),
            static_cast<void*>(&options.full),
            static_cast<void*>(&d_pair_counter),
            static_cast<void*>(&d_pair_indices),
            static_cast<void*>(&d_shifts),
            static_cast<void*>(&d_distances),
            static_cast<void*>(&d_vectors),
            static_cast<void*>(&options.return_shifts),
            static_cast<void*>(&options.return_distances),
            static_cast<void*>(&options.return_vectors),
            static_cast<void*>(&max_pairs),
            static_cast<void*>(&d_overflow_flag)
        };
        size_t THREADS_PER_PARTICLE = 8;
        size_t particles_per_block = THREADS_PER_BLOCK / THREADS_PER_PARTICLE;
        size_t num_blocks_find = (n_points + particles_per_block - 1) / particles_per_block;
        find_kernel->launch(
            dim3(num_blocks_find), dim3(THREADS_PER_BLOCK), 0, nullptr, find_args, false
        );
        NVTX_POP();

        NVTX_POP(); // cell_list_total
    }

    if (!use_cell_list) {
        NVTX_PUSH("brute_force_total");

        size_t THREADS_PER_BLOCK = 128;
        double cutoff2 = options.cutoff * options.cutoff;

        size_t num_half_pairs = n_points * (n_points - 1) / 2;

        if (is_orthogonal) {
            if (options.full) {
                NVTX_PUSH("brute_force_full_orthogonal");
                auto* kernel = factory.create(
                    "brute_force_full_orthogonal",
                    CUDA_BRUTEFORCE_CODE,
                    "cuda_bruteforce.cu",
                    {"-std=c++17"}
                );

                std::vector<void*> args = {
                    static_cast<void*>(&d_positions),
                    static_cast<void*>(&d_box_diag),
                    static_cast<void*>(&d_periodic),
                    static_cast<void*>(&n_points),
                    static_cast<void*>(&cutoff2),
                    static_cast<void*>(&d_pair_counter),
                    static_cast<void*>(&d_pair_indices),
                    static_cast<void*>(&d_shifts),
                    static_cast<void*>(&d_distances),
                    static_cast<void*>(&d_vectors),
                    static_cast<void*>(&options.return_shifts),
                    static_cast<void*>(&options.return_distances),
                    static_cast<void*>(&options.return_vectors),
                    static_cast<void*>(&max_pairs),
                    static_cast<void*>(&d_overflow_flag)
                };

                size_t num_blocks = (num_half_pairs + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
                kernel->launch(
                    /*grid=*/dim3(std::max(num_blocks, static_cast<size_t>(1))),
                    /*block=*/dim3(THREADS_PER_BLOCK),
                    /*shared_mem_size=*/0,
                    /*cuda_stream=*/nullptr,
                    /*args=*/args,
                    /*synchronize=*/false
                );
                NVTX_POP();
            } else {
                NVTX_PUSH("brute_force_half_orthogonal");
                auto* kernel = factory.create(
                    "brute_force_half_orthogonal",
                    CUDA_BRUTEFORCE_CODE,
                    "cuda_bruteforce.cu",
                    {"-std=c++17"}
                );

                std::vector<void*> args = {
                    static_cast<void*>(&d_positions),
                    static_cast<void*>(&d_box_diag),
                    static_cast<void*>(&d_periodic),
                    static_cast<void*>(&n_points),
                    static_cast<void*>(&cutoff2),
                    static_cast<void*>(&d_pair_counter),
                    static_cast<void*>(&d_pair_indices),
                    static_cast<void*>(&d_shifts),
                    static_cast<void*>(&d_distances),
                    static_cast<void*>(&d_vectors),
                    static_cast<void*>(&options.return_shifts),
                    static_cast<void*>(&options.return_distances),
                    static_cast<void*>(&options.return_vectors),
                    static_cast<void*>(&max_pairs),
                    static_cast<void*>(&d_overflow_flag)
                };

                size_t num_blocks = (num_half_pairs + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
                kernel->launch(
                    /*grid=*/dim3(std::max(num_blocks, static_cast<size_t>(1))),
                    /*block=*/dim3(THREADS_PER_BLOCK),
                    /*shared_mem_size=*/0,
                    /*cuda_stream=*/nullptr,
                    /*args=*/args,
                    /*synchronize=*/false
                );
                NVTX_POP();
            }
        } else {
            if (options.full) {
                NVTX_PUSH("brute_force_full_general");
                auto* kernel = factory.create(
                    "brute_force_full_general",
                    CUDA_BRUTEFORCE_CODE,
                    "cuda_bruteforce.cu",
                    {"-std=c++17"}
                );

                std::vector<void*> args = {
                    static_cast<void*>(&d_positions),
                    static_cast<void*>(&d_box),
                    static_cast<void*>(&d_inv_box_brute),
                    static_cast<void*>(&d_periodic),
                    static_cast<void*>(&n_points),
                    static_cast<void*>(&cutoff2),
                    static_cast<void*>(&d_pair_counter),
                    static_cast<void*>(&d_pair_indices),
                    static_cast<void*>(&d_shifts),
                    static_cast<void*>(&d_distances),
                    static_cast<void*>(&d_vectors),
                    static_cast<void*>(&options.return_shifts),
                    static_cast<void*>(&options.return_distances),
                    static_cast<void*>(&options.return_vectors),
                    static_cast<void*>(&max_pairs),
                    static_cast<void*>(&d_overflow_flag)
                };

                size_t num_blocks = (num_half_pairs + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
                kernel->launch(
                    /*grid=*/dim3(std::max(num_blocks, static_cast<size_t>(1))),
                    /*block=*/dim3(THREADS_PER_BLOCK),
                    /*shared_mem_size=*/0,
                    /*cuda_stream=*/nullptr,
                    /*args=*/args,
                    /*synchronize=*/false
                );
                NVTX_POP();
            } else {
                NVTX_PUSH("brute_force_half_general");
                auto* kernel = factory.create(
                    "brute_force_half_general",
                    CUDA_BRUTEFORCE_CODE,
                    "cuda_bruteforce.cu",
                    {"-std=c++17"}
                );

                std::vector<void*> args = {
                    static_cast<void*>(&d_positions),
                    static_cast<void*>(&d_box),
                    static_cast<void*>(&d_inv_box_brute),
                    static_cast<void*>(&d_periodic),
                    static_cast<void*>(&n_points),
                    static_cast<void*>(&cutoff2),
                    static_cast<void*>(&d_pair_counter),
                    static_cast<void*>(&d_pair_indices),
                    static_cast<void*>(&d_shifts),
                    static_cast<void*>(&d_distances),
                    static_cast<void*>(&d_vectors),
                    static_cast<void*>(&options.return_shifts),
                    static_cast<void*>(&options.return_distances),
                    static_cast<void*>(&options.return_vectors),
                    static_cast<void*>(&max_pairs),
                    static_cast<void*>(&d_overflow_flag)
                };

                size_t num_blocks = (num_half_pairs + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
                kernel->launch(
                    /*grid=*/dim3(std::max(num_blocks, static_cast<size_t>(1))),
                    /*block=*/dim3(THREADS_PER_BLOCK),
                    /*shared_mem_size=*/0,
                    /*cuda_stream=*/nullptr,
                    /*args=*/args,
                    /*synchronize=*/false
                );
                NVTX_POP();
            }
        }

        NVTX_POP(); // brute_force_total
    }

    NVTX_PUSH("async_copy_and_sync");

    CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemcpyAsync(
        extras->pinned_length_ptr,
        d_pair_counter,
        sizeof(size_t),
        cudaMemcpyDeviceToHost,
        nullptr
    ));

    CUDART_SAFE_CALL(CUDART_INSTANCE.cudaDeviceSynchronize());

    // Check for overflow
    int h_overflow_flag = 0;
    CUDART_SAFE_CALL(CUDART_INSTANCE.cudaMemcpy(
        &h_overflow_flag,
        d_overflow_flag,
        sizeof(int),
        cudaMemcpyDeviceToHost
    ));

    if (h_overflow_flag != 0) {
        throw std::runtime_error(
            "The number of neighbor pairs exceeds the maximum capacity of " +
            std::to_string(max_pairs) + " (max_pairs_per_point=" +
            std::to_string(max_pairs_per_point) + "; n_points=" +
            std::to_string(n_points) + "). " +
            "Consider reducing the cutoff distance, or explicitly setting " +
            "VESIN_CUDA_MAX_PAIRS_PER_POINT as an environment variable."
        );
    }

    neighbors.length = *extras->pinned_length_ptr;

    NVTX_POP();
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

            PLMD::metatomic::vesin::cpu::neighbors(
                reinterpret_cast<const PLMD::metatomic::vesin::Vector*>(points),
                n_points,
                PLMD::metatomic::vesin::BoundingBox(matrix, periodic),
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
