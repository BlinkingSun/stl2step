#!/usr/bin/env python3
"""Write a single-body A9 fixture: hex tube (2 cylinders) + sloped top (oblique)."""
import math
import os
import struct
import sys

N = 12
R_OUT = 10.0
R_IN = 4.0
Z0 = 0.0


def z_top(x, _y):
    # Slope in X: |n·â| = |n·z| is neither ~0 nor ~1.
    return 8.0 + 0.25 * x


def nrm(a, b, c):
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    L = math.sqrt(nx * nx + ny * ny + nz * nz)
    if L <= 0.0:
        return (0.0, 0.0, 1.0)
    return (nx / L, ny / L, nz / L)


def add_tri(tris, a, b, c):
    tris.append((nrm(a, b, c), a, b, c))


def ring(r, zfun):
    out = []
    for i in range(N):
        ang = 2.0 * math.pi * i / N
        x, y = r * math.cos(ang), r * math.sin(ang)
        out.append((x, y, zfun(x, y)))
    return out


def write_stl(path, tris):
    hdr = b"nonprismatic-control A9 hex-tube sloped-top"
    with open(path, "wb") as f:
        f.write(hdr + b"\0" * (80 - len(hdr)))
        f.write(struct.pack("<I", len(tris)))
        for n, a, b, c in tris:
            f.write(struct.pack("<12fH", n[0], n[1], n[2],
                                a[0], a[1], a[2], b[0], b[1], b[2],
                                c[0], c[1], c[2], 0))


def main():
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    out = os.path.join(root, "tests", "corpus", "nonprismatic-control.stl")
    if len(sys.argv) > 1:
        out = sys.argv[1]
    tris = []
    bot_o = ring(R_OUT, lambda x, y: Z0)
    bot_i = ring(R_IN, lambda x, y: Z0)
    top_o = ring(R_OUT, z_top)
    top_i = ring(R_IN, z_top)

    for i in range(N):
        j = (i + 1) % N
        # Bottom cap (outward -Z): outer[i] -> inner[i] -> inner[j] -> outer[j]
        add_tri(tris, bot_o[i], bot_i[i], bot_i[j])
        add_tri(tris, bot_o[i], bot_i[j], bot_o[j])
        # Sloped top (outward roughly +Z)
        add_tri(tris, top_o[i], top_o[j], top_i[j])
        add_tri(tris, top_o[i], top_i[j], top_i[i])
        # Outer wall
        add_tri(tris, bot_o[i], bot_o[j], top_o[j])
        add_tri(tris, bot_o[i], top_o[j], top_o[i])
        # Inner wall (inward)
        add_tri(tris, bot_i[i], top_i[i], top_i[j])
        add_tri(tris, bot_i[i], top_i[j], bot_i[j])

    os.makedirs(os.path.dirname(out), exist_ok=True)
    write_stl(out, tris)
    print("wrote", out, "tris", len(tris))


if __name__ == "__main__":
    main()
