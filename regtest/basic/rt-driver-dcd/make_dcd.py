"""Write the small dcd files of this regtest, one for each flavour of the format read by the dcd reader of PLUMED
(they are in the regtest: this script is kept to document how they were made).
usage: python3 make_dcd.py OUTDIR"""
import struct, sys, os, math
out = sys.argv[1]

def write(name, natoms=4, nframes=3, endian='<', marker='i', charmm=True, cell='cos',
          fixed=None, fourdims=False, truncate=0, xplor_delta=False):
    m = endian + marker
    def rec(payload):
        return struct.pack(m, len(payload)) + payload + struct.pack(m, len(payload))
    hdr = bytearray(80)
    def put_i(off, v): hdr[off:off+4] = struct.pack(endian+'i', v)
    put_i(0, nframes); put_i(4, 10); put_i(8, 5)
    nfixed = 0 if fixed is None else natoms - len(fixed)
    put_i(32, nfixed)
    if charmm:
        hdr[36:40] = struct.pack(endian+'f', 0.02)
        put_i(40, 1 if cell else 0)
        put_i(44, 1 if fourdims else 0)
        put_i(76, 24)
    else:
        hdr[36:44] = struct.pack(endian+'d', 0.02)
    data = rec(b'CORD' + bytes(hdr))
    titles = [b'REMARKS test file'.ljust(80), b'REMARKS for plumed'.ljust(80)]
    data += rec(struct.pack(endian+'i', len(titles)) + b''.join(titles))
    data += rec(struct.pack(endian+'i', natoms))
    if fixed is not None:
        data += rec(struct.pack(endian+'%di' % len(fixed), *[i+1 for i in fixed]))
    for f in range(nframes):
        if charmm and cell:
            a, b, c = 30.0+f, 32.5+0.1*f, 35.25
            if cell == 'cos':
                # gamma, beta, alpha as cosines (CHARMM, NAMD>2.5)
                ab, ac, bc = math.cos(math.radians(80.0+f)), math.cos(math.radians(95.5)), math.cos(math.radians(88.0))
            elif cell == 'matrix':
                # CHARMM (c36b2 and later) triclinic cell: the symmetric shape matrix h11, h21, h22, h31, h32, h33
                a, ab, b, ac, bc, c = 30.0+f, 14.5, 31.75, 9.625, -2.5, 32.5
            elif cell == 'ortho':
                ab = ac = bc = 0.0
            else:
                # NAMD 2.5: degrees
                ab, ac, bc = 80.0+f, 95.5, 88.0
            data += rec(struct.pack(endian+'6d', a, ab, b, ac, bc, c))
        n = natoms if (fixed is None or f == 0) else len(fixed)
        idx = list(range(natoms)) if (fixed is None or f == 0) else fixed
        for d in range(3):
            vals = [1.25*i + 10.0*d + 0.123456*f + 0.001*i*i for i in idx]
            data += rec(struct.pack(endian+'%df' % n, *vals))
        if fourdims:
            data += rec(struct.pack(endian+'%df' % n, *[0.5]*n))
    if truncate:
        data = data[:-truncate]
    with open(os.path.join(out, name), 'wb') as fh:
        fh.write(data)

write('charmm_le.dcd')
write('charmm_be.dcd', endian='>')
write('charmm_ortho.dcd', cell='ortho')
write('namd25_degrees.dcd', cell='deg')
write('nocell.dcd', cell=None)
write('charmm_matrix.dcd', cell='matrix')
write('xplor.dcd', charmm=False, cell=None)
write('xplor_be.dcd', charmm=False, cell=None, endian='>')
write('fixed_le.dcd', natoms=6, nframes=4, fixed=[1, 3, 4])
write('fixed_be.dcd', natoms=6, nframes=4, fixed=[1, 3, 4], endian='>')
write('fourdims.dcd', fourdims=True)
write('charmm_i8.dcd', marker='q')
write('charmm_i8_be.dcd', marker='q', endian='>')
write('fixed_i8.dcd', natoms=6, nframes=4, fixed=[0, 5], marker='q')
write('truncated.dcd', truncate=7)
write('truncated_fixed.dcd', natoms=6, nframes=4, fixed=[1, 3, 4], truncate=3)
