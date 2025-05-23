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
#include <cassert>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <numeric>
#include <tuple>
#include <new>

#ifndef VESIN_CPU_CELL_LIST_HPP
#define VESIN_CPU_CELL_LIST_HPP

#include <vector>

#include "vesin.h"

#ifndef VESIN_TYPES_HPP
#define VESIN_TYPES_HPP

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
        return (*this)[0][0] * ((*this)[1][1] * (*this)[2][2] - (*this)[2][1] * (*this)[1][2])
             - (*this)[0][1] * ((*this)[1][0] * (*this)[2][2] - (*this)[1][2] * (*this)[2][0])
             + (*this)[0][2] * ((*this)[1][0] * (*this)[2][1] - (*this)[1][1] * (*this)[2][0]);
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
    BoundingBox(Matrix matrix, bool periodic): matrix_(matrix), periodic_(periodic) {
        if (periodic) {
            auto det = matrix_.determinant();
            if (std::abs(det) < 1e-30) {
                throw std::runtime_error("the box matrix is not invertible");
            }

            this->inverse_ = matrix_.inverse();
        } else {
            this->matrix_ = Matrix{{{
                {{1, 0, 0}},
                {{0, 1, 0}},
                {{0, 0, 1}}
            }}};
            this->inverse_ = matrix_;
        }
    }

    const Matrix& matrix() const {
        return this->matrix_;
    }

    bool periodic() const {
        return this->periodic_;
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
    bool periodic_;
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
namespace vesin { namespace cpu {

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

} // namespace vesin
} // namespace metatomic
} // namespace PLMD
} // namespace cpu

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
    auto cell_list = CellList(cell, options.cutoff);

    for (size_t i=0; i<n_points; i++) {
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

                if ((shift[0] + shift[1] + shift[2] == 0)
                    && (shift[2] < 0 || (shift[2] == 0 && shift[1] < 0))) {
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
    box_(box)
{
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

    for (size_t spatial=0; spatial<3; spatial++) {
        if (n_search_[spatial] < 1) {
            n_search_[spatial] = 1;
        }

        // don't look for neighboring cells if we have only one cell and no
        // periodic boundary condition
        if (n_cells[spatial] == 1 && !box.periodic()) {
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

    // deal with pbc by wrapping the atom inside if it was outside of the
    // cell
    CellShift shift;
    // auto (shift, cell_index) =
    if (box_.periodic()) {
        auto result = divmod(cell_index, cells_shape_);
        shift = CellShift{std::get<0>(result)};
        cell_index = std::get<1>(result);
    } else {
        shift = CellShift({0, 0, 0});
        cell_index = std::array<int32_t, 3>{
            std::clamp(cell_index[0], 0, static_cast<int32_t>(cells_shape_[0] - 1)),
            std::clamp(cell_index[1], 0, static_cast<int32_t>(cells_shape_[1] - 1)),
            std::clamp(cell_index[2], 0, static_cast<int32_t>(cells_shape_[2] - 1)),
        };
    }

    this->get_cell(cell_index).emplace_back(Point{index, shift});
}


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

                    if (!box_.periodic() && !shift_is_zero) {
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
    auto* new_ptr = reinterpret_cast<scalar_t (*)[N]>(std::realloc(ptr, new_size * sizeof(scalar_t[N])));

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

    double *new_distances = nullptr;
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
        compare_pairs(size_t (*pairs_)[2]): pairs(pairs_) {}

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
    assert(neighbors.device == VesinCPU);

    std::free(neighbors.pairs);
    std::free(neighbors.shifts);
    std::free(neighbors.vectors);
    std::free(neighbors.distances);
}
#include <cstring>
#include <string>



thread_local std::string LAST_ERROR;

extern "C" int vesin_neighbors(
    const double (*points)[3],
    size_t n_points,
    const double box[3][3],
    bool periodic,
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

    if (neighbors->device != VesinUnknownDevice && neighbors->device != device) {
        *error_message = "`neighbors` device and data `device` do not match, free the neighbors first";
        return EXIT_FAILURE;
    }

    if (device == VesinUnknownDevice) {
        *error_message = "got an unknown device to use when running simulation";
        return EXIT_FAILURE;
    }

    if (neighbors->device == VesinUnknownDevice) {
        // initialize the device
        neighbors->device = device;
    } else if (neighbors->device != device) {
        *error_message = "`neighbors.device` and `device` do not match, free the neighbors first";
        return EXIT_FAILURE;
    }

    try {
        if (device == VesinCPU) {
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
        } else {
            throw std::runtime_error("unknown device " + std::to_string(device));
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

    if (neighbors->device == VesinUnknownDevice) {
        // nothing to do
    } else if (neighbors->device == VesinCPU) {
        PLMD::metatomic::vesin::cpu::free_neighbors(*neighbors);
    }

    std::memset(neighbors, 0, sizeof(VesinNeighborList));
}
