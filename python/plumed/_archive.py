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
"""The uncompressed PLUMED HDF5 USTAR batch profile (no extraction)."""

import io
import os
import tarfile
import warnings
from contextlib import contextmanager

from ._errors import FormatError


def recognizes(file):
    if not isinstance(file, (str, os.PathLike)):
        return False
    try:
        file = os.path.realpath(os.fspath(file))
        with open(file, "rb") as stream:
            header = stream.read(1024)
            if header[257:262] == b"ustar" or header == bytes(1024):
                return True
        with open(os.fspath(file) + ".plumed-trim", "rb") as recovery:
            return recovery.read(17) == b"plumed-h5-trim-v1 "
    except OSError:
        return False


def _number(value):
    if not value or not value.isdecimal():
        raise FormatError("Invalid PLUMED archive block/row metadata")
    result = int(value)
    if result > 2**64 - 1:
        raise FormatError("PLUMED archive counter overflow")
    return result


@contextmanager
def members(path):
    """Snapshot the headers, then yield one closed HDF5 image at a time.

    The optional durable trim record overlays only the replacement boundary
    batch; this also permits reading after a process dies during truncation.
    """
    path = os.path.realpath(os.fspath(path))
    with open(path, "rb") as stream:
        recovery = None
        try:
            try:
                recovery = open(os.fspath(path) + ".plumed-trim", "rb")
            except FileNotFoundError:
                pass
            length = os.fstat(stream.fileno()).st_size
            replacement = None
            if recovery is not None:
                record = recovery.read(512).split(b"\0", 1)[0].decode("ascii").split()
                if len(record) != 6 or record[0] != "plumed-h5-trim-v1":
                    raise FormatError("Invalid trim recovery record")
                prefix, count, size, block, rows = map(_number, record[1:])
                if (
                    prefix > length
                    or prefix % 4096
                    or size != os.fstat(recovery.fileno()).st_size - 512
                ):
                    raise FormatError("Inconsistent trim recovery record")
                length = prefix
                replacement = (count, size, block, rows)
            index = []
            offset = 0
            previous = 0
            while length - offset >= 512:
                stream.seek(offset)
                header = stream.read(512)
                if header == bytes(512):
                    break
                try:
                    member = tarfile.TarInfo.frombuf(header, "utf-8", "strict")
                except (tarfile.HeaderError, UnicodeError) as exc:
                    if index:
                        warnings.warn(
                            "Ignoring damaged/unpublished HDF5 archive tail",
                            RuntimeWarning,
                        )
                        break
                    raise FormatError("Invalid first TAR batch header") from exc
                if (
                    header[257:262] != b"ustar"
                    or not member.isfile()
                    or member.uname != "plumed-h5-v1"
                ):
                    raise FormatError("Unsupported PLUMED archive version/profile")
                if member.name != "batch_" + str(len(index) + 1) + ".h5":
                    raise FormatError(
                        "Archive members must be consecutively numbered batch_N.h5"
                    )
                values = member.gname.split(":")
                if len(values) != 2:
                    raise FormatError("Invalid batch block/row metadata")
                block, rows = map(_number, values)
                if (
                    not rows
                    or block not in (previous, previous + 1)
                    or (not index and block)
                ):
                    raise FormatError("Invalid logical block sequence")
                rounded = (member.size + 511) // 512 * 512
                if member.size <= 0:
                    raise FormatError("Invalid TAR member size")
                if rounded > length - offset - 512:
                    warnings.warn(
                        "Ignoring incomplete HDF5 archive tail", RuntimeWarning
                    )
                    break
                index.append((stream, offset + 512, member.size, block, rows))
                previous = block
                offset += 512 + rounded
            if 0 < length - offset < 512:
                warnings.warn("Ignoring incomplete HDF5 archive tail", RuntimeWarning)
            if replacement is not None:
                count, size, block, rows = replacement
                if offset != length or len(index) != count:
                    raise FormatError("Trim prefix does not match recovery record")
                if size:
                    index.append((recovery, 512, size, block, rows))

            def images():
                for source, offset, size, block, rows in index:
                    source.seek(offset)
                    image = source.read(size)
                    if len(image) != size:
                        raise FormatError("Truncated HDF5 batch")
                    yield image, block, rows

            yield images()
        finally:
            if recovery is not None:
                recovery.close()


def blocks(path):
    from ._hdf5 import _blocks, _h5py

    previous = None
    previous_block = None
    with members(path) as images:
        for image, block, rows in images:
            with _h5py().File(io.BytesIO(image), "r") as handle:
                descriptions = _blocks(handle)
                if len(descriptions) != 1 or descriptions[0][3] != rows:
                    raise FormatError("HDF5 batch disagrees with its TAR header")
                group, names, constants, _ = descriptions[0]
                signature = (
                    tuple((name, group[name].dtype.kind) for name in names),
                    tuple(
                        (name, type(value).__name__, repr(value))
                        for name, value, _ in constants
                    ),
                )
                if block == previous_block and signature != previous:
                    raise FormatError(
                        "Schema/constants change within a logical archive block"
                    )
                previous, previous_block = signature, block
                yield descriptions[0], block


@contextmanager
def _writer_lock(path):
    descriptor = os.open(path + ".plumed-lock", os.O_CREAT | os.O_RDWR, 0o666)
    try:
        if os.name == "nt":
            import msvcrt

            msvcrt.locking(descriptor, msvcrt.LK_NBLCK, 1)
        else:
            import fcntl

            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield
    finally:
        os.close(descriptor)


def write(df, path):
    """Atomically export a dataframe as ordinary HDF5 batches in a USTAR file."""
    import stat
    import tempfile

    from . import _hdf5

    destination = os.path.realpath(os.fspath(path))
    parent = os.path.dirname(os.path.abspath(destination))
    with (
        _writer_lock(destination),
        tempfile.TemporaryDirectory(prefix="plumed-h5-", dir=parent) as directory,
    ):
        temporary = os.path.join(directory, "table.tar")
        with open(temporary, "wb") as stream:
            with tarfile.open(
                fileobj=stream, mode="w", format=tarfile.USTAR_FORMAT
            ) as archive:
                for first in range(0, len(df), 4096):
                    frame = df.iloc[first : first + 4096]
                    if hasattr(df, "plumed_constants"):
                        frame.plumed_constants = df.plumed_constants
                    image = io.BytesIO()
                    _hdf5.write(frame, image)
                    data = image.getvalue()
                    size = (len(data) + 512 + 4095) // 4096 * 4096 - 512
                    member = tarfile.TarInfo("batch_" + str(first // 4096 + 1) + ".h5")
                    member.uname = "plumed-h5-v1"
                    member.gname = "0:" + str(len(frame))
                    member.size = size
                    archive.addfile(member, io.BytesIO(data + bytes(size - len(data))))
            stream.flush()
            try:
                os.chmod(temporary, stat.S_IMODE(os.stat(destination).st_mode))
            except FileNotFoundError:
                pass
            os.fsync(stream.fileno())
            if os.sys.platform == "darwin":
                import fcntl

                fcntl.fcntl(stream.fileno(), fcntl.F_FULLFSYNC)
        os.replace(temporary, destination)
        if os.name != "nt":
            descriptor = os.open(parent, os.O_RDONLY)
            try:
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
