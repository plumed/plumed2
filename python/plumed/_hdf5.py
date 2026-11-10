# /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   Copyright (c) 2026 The plumed team
#   (see the PEOPLE file at the root of the distribution for a list of names)
#
#   See http://www.plumed.org for more information.
#
#   This file is part of plumed, version 2.
#
#   plumed is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   plumed is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#
"""PLUMED HDF5 table version 1, independent of the native PLUMED kernel.

h5py and numpy are imported only when an HDF5 operation is requested.
"""

import operator
import os
from contextlib import ExitStack
from numbers import Integral

from ._io import FormatError


def is_hdf5(file, format):
    if format not in (None, "text", "hdf5"):
        raise ValueError("format must be None, 'text', or 'hdf5'")
    if format is not None:
        return format == "hdf5"
    if type(file).__module__.startswith("h5py."):
        return True
    name = (
        os.fspath(file)
        if isinstance(file, (str, os.PathLike))
        else getattr(file, "name", "")
    )
    return isinstance(name, str) and name.lower().endswith((".h5", ".hdf5"))


def _h5py():
    try:
        import h5py
    except ImportError as exc:
        raise ImportError(
            "PLUMED HDF5 support requires h5py; install it with 'pip install h5py'"
        ) from exc
    return h5py


def _open(stack, file, mode):
    h5py = _h5py()
    if isinstance(file, h5py.File):
        return file  # Borrowed: the caller owns its lifetime.
    if file is None:
        raise ValueError("HDF5 output requires a path, binary stream, or h5py.File")
    return stack.enter_context(h5py.File(file, mode))


def _text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _names(value):
    if not isinstance(value, (str, bytes)):
        raise FormatError("HDF5 field lists must be strings")
    names = _text(value).split()
    if len(names) != len(set(names)):
        raise FormatError("Duplicate HDF5 field names")
    return names


def _blocks(file):
    h5py = _h5py()
    if _text(file.attrs.get("plumed_format", "")) != "plumed-table":
        raise FormatError("Not a PLUMED HDF5 table")
    if file.attrs.get("plumed_format_version") != 1:
        raise FormatError("Unsupported PLUMED HDF5 table version")
    if set(file) != {"block_" + str(i) for i in range(len(file))}:
        raise FormatError(
            "PLUMED HDF5 blocks must be numbered consecutively from block_0"
        )
    result = []
    try:
        for i in range(len(file)):
            group = file["block_" + str(i)]
            if not isinstance(group, h5py.Group):
                raise FormatError("PLUMED HDF5 block is not a group")
            names = _names(group.attrs["plumed_columns"])
            constants = []
            for name in _names(group.attrs["plumed_constants"]):
                value = group.attrs[name]
                if getattr(value, "ndim", 0) != 0:
                    raise FormatError("PLUMED HDF5 constants must be scalars")
                if isinstance(value, bytes):
                    value = value.decode("utf-8")
                if hasattr(value, "item"):
                    value = value.item()
                constants.append((name, value, str(value)))
            lengths = set()
            for name in names:
                ds = group[name]
                if not isinstance(ds, h5py.Dataset) or ds.ndim != 1:
                    raise FormatError(
                        "PLUMED HDF5 columns must be one-dimensional datasets"
                    )
                if (
                    ds.dtype.kind not in "if"
                    and h5py.check_string_dtype(ds.dtype) is None
                ):
                    raise FormatError("Unsupported PLUMED HDF5 column type: " + name)
                lengths.add(len(ds))
            if len(lengths) > 1:
                raise FormatError("PLUMED HDF5 columns have inconsistent lengths")
            result.append((group, names, constants, next(iter(lengths), 0)))
    except (KeyError, TypeError, UnicodeError) as exc:
        raise FormatError("Malformed PLUMED HDF5 table: " + str(exc)) from exc
    return result


def _integer(value, name, positive=False):
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise ValueError(name + " must be an integer") from exc
    if isinstance(value, bool) or result < (1 if positive else 0):
        raise ValueError(
            name + (" must be positive" if positive else " must be non-negative")
        )
    return result


def _columns(names, usecols):
    if usecols is None:
        return names
    if callable(usecols):
        return [name for name in names if usecols(name)]
    requested = list(usecols)
    if all(isinstance(item, str) for item in requested):
        if not set(requested).issubset(names):
            raise ValueError("usecols contains unknown HDF5 column names")
        return [name for name in names if name in requested]
    if all(
        isinstance(item, Integral) and not isinstance(item, bool) for item in requested
    ):
        if any(item < 0 or item >= len(names) for item in requested):
            raise ValueError("usecols contains out-of-range HDF5 column indices")
        return [name for i, name in enumerate(names) if i in requested]
    raise ValueError("usecols must contain only names or only integer positions")


class _Reader:
    def __init__(
        self,
        file,
        *,
        chunksize,
        usecols,
        skiprows,
        nrows,
        index_col,
        enable_constants,
        process,
    ):
        self.stack = ExitStack()
        self._closed = False
        try:
            self.file = _open(self.stack, file, "r")
            self.blocks = _blocks(self.file)
            if usecols is not None and not callable(usecols):
                usecols = list(usecols)
            self.selected = [_columns(names, usecols) for _, names, _, _ in self.blocks]
            self.size = (
                _integer(chunksize, "chunksize", True)
                if chunksize is not None
                else 65536
            )
            self.limit = _integer(nrows, "nrows") if nrows is not None else None
            self.skip_count, self.skip_set, self.skip_func = 0, set(), None
            if callable(skiprows):
                self.skip_func = skiprows
            elif skiprows is not None:
                try:
                    operator.index(skiprows)
                except TypeError:
                    self.skip_set = {_integer(i, "skiprows entry") for i in skiprows}
                else:
                    self.skip_count = _integer(skiprows, "skiprows")
            if index_col is True:
                raise ValueError("index_col=True is not a valid column selection")
            if (
                index_col is not None
                and index_col is not False
                and not isinstance(index_col, (str, Integral))
            ):
                try:
                    index_col = list(index_col)
                except TypeError as exc:
                    raise ValueError(
                        "index_col must be a name, position, or sequence"
                    ) from exc
                if not index_col:
                    index_col = None
            self.index_col = index_col
            self.mode, self.process = enable_constants, process
            self.offset = self.emitted = 0
            self._iterator = self._frames()
        except BaseException:
            self.close()
            raise

    def _frames(self):
        import numpy as np
        import pandas as pd

        h5py = _h5py()
        simple_skip = self.skip_func is None and not self.skip_set
        for (group, names, constants, length), selected in zip(
            self.blocks, self.selected
        ):

            def finish(df, constants=constants):
                if self.index_col is not None and self.index_col is not False:
                    columns = list(df.columns)
                    keys = (
                        self.index_col
                        if isinstance(self.index_col, (list, tuple))
                        else [self.index_col]
                    )
                    keys = [columns[k] if isinstance(k, Integral) else k for k in keys]
                    df = df.set_index(keys)
                return self.process(df, constants if self.mode != "no" else [])

            pending = None
            pos = 0
            if self.limit != 0 and simple_skip:
                pos = min(length, max(0, self.skip_count - self.offset))
            while pos < length or (length == 0 and pos == 0):
                if (
                    self.limit is not None
                    and self.emitted >= self.limit
                    and self.limit != 0
                ):
                    break
                end = pos if self.limit == 0 else min(pos + self.size, length)
                mask = None
                if simple_skip:
                    if self.limit is not None:
                        end = min(end, pos + self.limit - self.emitted)
                else:
                    indices = range(self.offset + pos, self.offset + end)
                    mask = np.fromiter(
                        (
                            i >= self.skip_count
                            and i not in self.skip_set
                            and (self.skip_func is None or not self.skip_func(i))
                            for i in indices
                        ),
                        bool,
                    )
                    if self.limit is not None and self.limit != 0:
                        remaining = self.limit - self.emitted
                        positions = np.flatnonzero(mask)
                        if len(positions) >= remaining:
                            end = pos + int(positions[remaining - 1]) + 1
                            mask = mask[: end - pos]
                    if not mask.any() and length and self.limit != 0:
                        pos = end
                        continue
                data = {}
                for name in selected:
                    ds = group[name]
                    data[name] = (
                        ds.asstr(encoding="utf-8")[pos:end]
                        if h5py.check_string_dtype(ds.dtype)
                        else ds[pos:end]
                    )
                df = pd.DataFrame(data, index=pd.RangeIndex(end - pos))
                if mask is not None:
                    df = df.loc[mask]
                df.index = pd.RangeIndex(self.emitted, self.emitted + len(df))
                self.emitted += len(df)
                if len(df):
                    pending = df if pending is None else pd.concat([pending, df])
                    if len(pending) >= self.size:
                        yield finish(pending.iloc[: self.size])
                        pending = pending.iloc[self.size :]
                elif length == 0 or self.limit == 0:
                    yield finish(df)
                pos = end
                if length == 0 or self.limit == 0:
                    break
            if pending is not None and len(pending):
                yield finish(pending)
            self.offset += length
            if self.limit is not None and self.emitted >= self.limit:
                return

    def __iter__(self):
        return self

    def __next__(self):
        if self._closed:
            raise StopIteration
        try:
            return next(self._iterator)
        except BaseException:
            self.close()
            raise

    def close(self):
        if self._closed:
            return
        self._closed = True
        try:
            iterator = getattr(self, "_iterator", None)
            if iterator is not None:
                iterator.close()
        finally:
            self.stack.close()

    def __del__(self):
        if not getattr(self, "_closed", True):
            self.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def read(file, **kwargs):
    import pandas as pd

    reader = _Reader(file, **kwargs)
    if kwargs["chunksize"] is not None:
        return reader
    with reader:
        frames = []
        constants = None
        for frame in reader:
            if kwargs["enable_constants"] == "metadata":
                current = frame.plumed_constants
                if constants is not None and current != constants:
                    raise FormatError(
                        "HDF5 constants change between blocks; use enable_constants='columns', "
                        "enable_constants=False, or chunksize"
                    )
                constants = current
            frames.append(frame)
        if not frames:
            columns = reader.selected[0] if reader.selected else []
            raw = pd.DataFrame(columns=columns)
            if reader.index_col is not None and reader.index_col is not False:
                keys = (
                    reader.index_col
                    if isinstance(reader.index_col, (list, tuple))
                    else [reader.index_col]
                )
                raw = raw.set_index(
                    [columns[k] if isinstance(k, Integral) else k for k in keys]
                )
            values = (
                reader.blocks[0][2] if reader.blocks and reader.mode != "no" else []
            )
            return reader.process(raw, values)
        # pandas normally inserts float NaNs for absent integer columns, losing
        # values above 2**53. Keep integers nullable when schemas change, and
        # use objects if a field's type changes between integer and float.
        columns = dict.fromkeys(name for frame in frames for name in frame.columns)
        for name in columns:
            present = [frame for frame in frames if name in frame.columns]
            kinds = {frame[name].dtype.kind for frame in present}
            if kinds <= {"i"} and len(present) != len(frames):
                for frame in present:
                    frame[name] = frame[name].astype("Int64")
            elif "i" in kinds and "f" in kinds:
                for frame in present:
                    frame[name] = frame[name].astype(object)
        index_parts = [
            [frame.index.get_level_values(level) for frame in frames]
            for level in range(frames[0].index.nlevels)
        ]
        mixed_levels = [
            {index.dtype.kind for index in indexes} >= {"i", "f"}
            for indexes in index_parts
        ]
        result = pd.concat(frames)
        if any(mixed_levels):
            index_arrays = [
                pd.Index([v for index in indexes for v in index], dtype=object)
                if mixed
                else indexes[0].append(indexes[1:])
                for indexes, mixed in zip(index_parts, mixed_levels)
            ]
            names = [
                parts[0].name
                if all(index.name == parts[0].name for index in parts)
                else None
                for parts in index_parts
            ]
            if len(index_arrays) == 1:
                result.index = index_arrays[0].rename(names[0])
            else:
                result.index = pd.MultiIndex.from_arrays(index_arrays, names=names)
        if constants is not None:
            result.plumed_constants = constants
        return result


def _field_name(value):
    if (
        not isinstance(value, str)
        or not value
        or any(c.isspace() for c in value)
        or "/" in value
        or "\x00" in value
        or value in (".", "..")
    ):
        raise ValueError(
            "HDF5 field names must be nonempty strings without whitespace, '/' or NUL"
        )
    value.encode("utf-8")  # Validate before an output file can be truncated.
    return value


def write(df, file):
    import numpy as np
    import pandas as pd

    h5py = _h5py()
    if isinstance(df.index, pd.MultiIndex):
        raise TypeError("Writing dataframes with MultiIndexes is not supported")
    names = [_field_name(name) for name in df.columns]
    columns = [df.iloc[:, i] for i in range(len(names))]
    if df.index.name is not None:
        names.insert(0, _field_name(df.index.name))
        columns.insert(0, df.index)
    if len(set(names)) != len(names) or not names:
        raise ValueError(
            "HDF5 output requires unique column names and at least one column"
        )
    # Validate everything before opening (and potentially replacing) the destination.
    datasets = []
    for name, column in zip(names, columns):
        kind = getattr(column.dtype, "kind", "O")
        if kind in "biu":
            if pd.isna(column).any():
                raise ValueError(
                    "Missing integer values cannot be stored in PLUMED HDF5 tables"
                )
            if (
                kind == "u"
                and len(column)
                and int(column.max()) > np.iinfo(np.int64).max
            ):
                raise OverflowError("PLUMED HDF5 integers must fit in signed 64 bits")
            data, dtype = np.asarray(column, dtype="<i8"), "<i8"
        elif kind == "f":
            data, dtype = np.asarray(column, dtype="<f8"), "<f8"
        elif kind in "OUS" and all(
            isinstance(value, str) and "\x00" not in value for value in column
        ):
            data = np.asarray([value.encode("utf-8") for value in column], dtype=object)
            dtype = h5py.string_dtype("ascii")
        else:
            raise TypeError(
                "Unsupported PLUMED HDF5 column type for "
                + name
                + ": "
                + str(column.dtype)
            )
        datasets.append((name, data, dtype))
    constants = {}
    for constant in getattr(df, "plumed_constants", []):
        if len(constant) not in (2, 3):
            raise ValueError("Constants must have two or three entries")
        name = _field_name(constant[0])
        if name.startswith("plumed_") or name in constants:
            raise ValueError("Duplicate or reserved HDF5 constant name: " + name)
        value = constant[-1]
        # A numeric HDF5 constant has its normal representation in entry 2.
        # Retain that type on round trips; expressions/other original spellings
        # still take precedence over converted numeric values in entry 1.
        if (
            len(constant) == 3
            and type(constant[1]) in (int, float)
            and str(constant[1]) == value
        ):
            value = constant[1]
        if not isinstance(value, (str, int, float, np.integer, np.floating)):
            raise TypeError("Unsupported HDF5 constant type: " + name)
        if isinstance(value, str):
            if "\x00" in value:
                raise ValueError("HDF5 string constants cannot contain NUL")
            value.encode("utf-8")
        elif isinstance(value, (int, np.integer)):
            value = int(value)
            if not np.iinfo(np.int64).min <= value <= np.iinfo(np.int64).max:
                raise OverflowError(
                    "PLUMED HDF5 integer constants must fit in signed 64 bits"
                )
        else:
            value = float(value)
        constants[name] = value
    with ExitStack() as stack:
        handle = _open(stack, file, "w")
        if len(handle) or len(handle.attrs):
            raise ValueError("A supplied h5py.File must be empty")

        def attribute(obj, name, value):
            # Native H5Table uses variable-length H5T_C_S1 (ASCII-tagged).
            # Store UTF-8 bytes in the same type: HDF5 cannot convert a UTF-8
            # string type to the ASCII memory type requested by native readers.
            if isinstance(value, str):
                obj.attrs.create(
                    name, value.encode("utf-8"), dtype=h5py.string_dtype("ascii")
                )
            else:
                obj.attrs[name] = value

        attribute(handle, "plumed_format", "plumed-table")
        handle.attrs["plumed_format_version"] = np.int64(1)
        block = handle.create_group("block_0")
        attribute(block, "plumed_columns", " ".join(names))
        attribute(block, "plumed_constants", " ".join(constants))
        for name, value in constants.items():
            attribute(block, name, value)
        for name, data, dtype in datasets:
            block.create_dataset(
                name, data=data, dtype=dtype, maxshape=(None,), chunks=True
            )
        handle.flush()
