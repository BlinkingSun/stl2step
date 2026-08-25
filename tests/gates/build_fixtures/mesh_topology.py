"""Mesh topology matching tests/harness/mesh_harness.cpp (I5 local ids)."""
from __future__ import annotations

import math
import struct
from collections import deque
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

Vec3 = Tuple[float, float, float]
Tri = Tuple[int, int, int]


def dbl_bits(v: float) -> int:
    if v == 0.0:
        v = 0.0
    return struct.unpack("<Q", struct.pack("<d", v))[0]


def weld_stl_triangles(
    raw_tris: Sequence[Tuple[Vec3, Vec3, Vec3]],
) -> Tuple[List[Vec3], List[Tri]]:
    """Exact-bit weld matching mesh_harness.cpp (weldTol=0)."""
    pts: List[Vec3] = []
    remap: Dict[Tuple[int, int, int], int] = {}
    tris: List[Tri] = []

    def local_vid(p: Vec3) -> int:
        k = (dbl_bits(p[0]), dbl_bits(p[1]), dbl_bits(p[2]))
        if k not in remap:
            remap[k] = len(pts)
            pts.append(p)
        return remap[k]

    for a, b, c in raw_tris:
        ia, ib, ic = local_vid(a), local_vid(b), local_vid(c)
        if ia == ib or ib == ic or ia == ic:
            continue
        tris.append((ia, ib, ic))
    return pts, tris


def read_binary_stl(path: str) -> Tuple[List[Vec3], List[Tri]]:
    raw: List[Tuple[Vec3, Vec3, Vec3]] = []
    with open(path, "rb") as f:
        f.read(80)
        (n,) = struct.unpack("<I", f.read(4))
        for _ in range(n):
            f.read(12)
            a = struct.unpack("<3f", f.read(12))
            b = struct.unpack("<3f", f.read(12))
            c = struct.unpack("<3f", f.read(12))
            f.read(2)
            raw.append((a, b, c))
    return weld_stl_triangles(raw)


def vadd(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vsub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vscale(a: Vec3, s: float) -> Vec3:
    return (a[0] * s, a[1] * s, a[2] * s)


def vdot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def vcross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def vnorm(a: Vec3) -> float:
    return math.sqrt(vdot(a, a))


def vunit(a: Vec3) -> Vec3:
    n = vnorm(a)
    if n < 1e-30:
        return (0.0, 0.0, 0.0)
    return (a[0] / n, a[1] / n, a[2] / n)


def tri_area(pts: Sequence[Vec3], tri: Tri) -> float:
    a, b, c = (pts[tri[0]], pts[tri[1]], pts[tri[2]])
    return 0.5 * vnorm(vcross(vsub(b, a), vsub(c, a)))


def tri_normal(pts: Sequence[Vec3], tri: Tri) -> Vec3:
    a, b, c = (pts[tri[0]], pts[tri[1]], pts[tri[2]])
    return vunit(vcross(vsub(b, a), vsub(c, a)))


def ekey(a: int, b: int) -> Tuple[int, int]:
    return (a, b) if a < b else (b, a)


@dataclass
class Component:
    pts: List[Vec3]
    tris: List[Tri]
    comp_edges: List[Tuple[int, int]] = field(default_factory=list)
    tri_edges: List[Tuple[int, int, int]] = field(default_factory=list)
    tri_dirs: List[int] = field(default_factory=list)
    root: int = 0
    vol: float = 0.0
    open: int = 0
    conflict: int = 0
    non_manifold: int = 0

    @property
    def n_tri(self) -> int:
        return len(self.tris)

    @property
    def n_edge(self) -> int:
        return len(self.comp_edges)

    def build_topology(self) -> None:
        edge_map: Dict[Tuple[int, int], int] = {}
        edge_recs: List[Dict[str, int]] = []
        self.comp_edges = []
        self.tri_edges = []
        self.tri_dirs = []
        self.vol = 0.0

        for tri in self.tris:
            v = [tri[0], tri[1], tri[2], tri[0]]
            te: List[int] = []
            dirs = 0
            for s in range(3):
                gu, gv = v[s], v[s + 1]
                k = ekey(gu, gv)
                if k not in edge_map:
                    eid = len(self.comp_edges)
                    edge_map[k] = eid
                    self.comp_edges.append(k)
                    edge_recs.append({"cnt": 0, "fwd": 0})
                else:
                    eid = edge_map[k]
                edge_recs[eid]["cnt"] += 1
                if gu < gv:
                    edge_recs[eid]["fwd"] += 1
                    dirs |= 1 << s
                te.append(eid)
            self.tri_edges.append((te[0], te[1], te[2]))
            self.tri_dirs.append(dirs)
            a, b, c = (self.pts[tri[0]], self.pts[tri[1]], self.pts[tri[2]])
            self.vol += (
                a[0] * (b[1] * c[2] - b[2] * c[1])
                - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])
            ) / 6.0

        self.open = self.conflict = self.non_manifold = 0
        for rec in edge_recs:
            if rec["cnt"] == 1:
                self.open += 1
            elif rec["cnt"] > 2:
                self.non_manifold += 1
            elif rec["fwd"] != 1:
                self.conflict += 1

    def diag(self) -> float:
        if not self.pts:
            return 0.0
        xs = [p[0] for p in self.pts]
        ys = [p[1] for p in self.pts]
        zs = [p[2] for p in self.pts]
        mn = (min(xs), min(ys), min(zs))
        mx = (max(xs), max(ys), max(zs))
        return vnorm(vsub(mx, mn))


def write_binary_stl(path: str, pts: Sequence[Vec3], tris: Sequence[Tri]) -> None:
    with open(path, "wb") as f:
        header = b"stl2step build_fixtures" + b"\0" * (80 - 23)
        f.write(header[:80])
        f.write(struct.pack("<I", len(tris)))
        for tri in tris:
            a, b, c = (pts[tri[0]], pts[tri[1]], pts[tri[2]])
            n = vcross(vsub(b, a), vsub(c, a))
            n = vunit(n)
            f.write(struct.pack("<3f", *n))
            f.write(struct.pack("<3f", *a))
            f.write(struct.pack("<3f", *b))
            f.write(struct.pack("<3f", *c))
            f.write(struct.pack("<H", 0))


def cylinder_verts(r: float, z0: float, z1: float, n: int, phase: float = 0.0) -> List[Vec3]:
    out: List[Vec3] = []
    for i in range(n):
        ang = phase + 2.0 * math.pi * i / n
        out.append((r * math.cos(ang), r * math.sin(ang), z0))
    for i in range(n):
        ang = phase + 2.0 * math.pi * i / n
        out.append((r * math.cos(ang), r * math.sin(ang), z1))
    return out


def cylinder_tris(n: int, closed: bool = True) -> List[Tri]:
    tris: List[Tri] = []
    if closed:
        for i in range(n):
            j = (i + 1) % n
            tris.append((i, j, i + n))
            tris.append((j, j + n, i + n))
    else:
        for i in range(n - 1):
            j = i + 1
            tris.append((i, j, i + n))
            tris.append((j, j + n, i + n))
    return tris


def rect_plane(z: float, x0: float, x1: float, y0: float, y1: float, base: int = 0) -> Tuple[List[Vec3], List[Tri]]:
    pts = [(x0, y0, z), (x1, y0, z), (x1, y1, z), (x0, y1, z)]
    tris = [(base, base + 1, base + 2), (base, base + 2, base + 3)]
    return pts, tris


def merge_components(parts: List[Tuple[List[Vec3], List[Tri]]]) -> Component:
    pts: List[Vec3] = []
    tris: List[Tri] = []
    for p, t in parts:
        off = len(pts)
        pts.extend(p)
        tris.extend((a + off, b + off, c + off) for a, b, c in t)
    c = Component(pts=pts, tris=tris)
    c.build_topology()
    return c


def edge_between(comp: Component, va: int, vb: int) -> int:
    return {ekey(comp.comp_edges[i][0], comp.comp_edges[i][1]): i for i in range(comp.n_edge)}[
        ekey(va, vb)
    ]


def walk_ring(comp: Component, verts: Sequence[int]) -> Tuple[List[int], List[int]]:
    """Return (meshEdges, meshVerts) for a closed vertex ring."""
    edges: List[int] = []
    for i in range(len(verts)):
        a, b = verts[i], verts[(i + 1) % len(verts)]
        edges.append(edge_between(comp, a, b))
    return edges, list(verts)


def walk_open(comp: Component, verts: Sequence[int]) -> Tuple[List[int], List[int]]:
    edges = [edge_between(comp, verts[i], verts[i + 1]) for i in range(len(verts) - 1)]
    return edges, list(verts)


def cyl_dvol(r: float, n_sides: int, axial: float, outward: bool, n_bands: Optional[int] = None) -> float:
    """D4 circular-segment prism. n_bands defaults to n_sides (full 360°)."""
    if n_sides < 3 or r <= 0:
        return 0.0
    gamma = 2.0 * math.pi / n_sides
    if gamma >= math.pi:
        return 0.0
    sigma = 1.0 if outward else -1.0
    bands = n_sides if n_bands is None else n_bands
    # A_reg = bands * 2 * R * H * sin(gamma/2)
    # dVol  = sigma * A_reg * R * (gamma - sin gamma) / (4 * sin(gamma/2))
    #       = sigma * bands * H * (R^2/2) * (gamma - sin gamma)
    return sigma * bands * axial * (r * r / 2.0) * (gamma - math.sin(gamma))


def plane_dvol(pts: Sequence[Vec3], tris: Sequence[Tri], ax_dir: Vec3, ax_loc: Vec3) -> float:
    d = vunit(ax_dir)
    total = 0.0
    for tri in tris:
        a = tri_area(pts, tri)
        h = sum(vdot(d, vsub(pts[v], ax_loc)) for v in tri) / 3.0
        total -= a * h
    return total


def fillet_dvol(r: float, length: float, outward: bool, span: float = math.pi / 2) -> float:
    """D4 circular-segment prism for a fillet strip of arc `span` and axial `length`."""
    if r <= 0 or length <= 0 or span <= 0:
        return 0.0
    sigma = 1.0 if outward else -1.0
    # exact circular segment area * length: (R^2/2)(span - sin span) * L
    return sigma * (r * r / 2.0) * (span - math.sin(span)) * length


class MeshBuilder:
    """Exact-bit welded vertex table; triangle order is the STL / local-id order."""

    def __init__(self) -> None:
        self.pts: List[Vec3] = []
        self.tris: List[Tri] = []
        self._key: Dict[Tuple[int, int, int], int] = {}

    def v(self, x: float, y: float, z: float) -> int:
        p: Vec3 = (float(x), float(y), float(z))
        k = (dbl_bits(p[0]), dbl_bits(p[1]), dbl_bits(p[2]))
        i = self._key.get(k)
        if i is None:
            i = len(self.pts)
            self._key[k] = i
            self.pts.append(p)
        return i

    def tri(self, a: int, b: int, c: int) -> int:
        if a == b or b == c or a == c:
            raise ValueError(f"degenerate triangle {(a, b, c)}")
        self.tris.append((a, b, c))
        return len(self.tris) - 1

    def quad(self, a: int, b: int, c: int, d: int) -> List[int]:
        return [self.tri(a, b, c), self.tri(a, c, d)]

    def component(self) -> Component:
        c = Component(pts=list(self.pts), tris=list(self.tris))
        c.build_topology()
        return c

    def orient_manifold(self) -> None:
        """Flood-fill so every cnt==2 edge is walked opposite ways (conflict==0)."""
        n = len(self.tris)
        if n == 0:
            return

        def edge_map() -> Dict[Tuple[int, int], List[Tuple[int, int, int]]]:
            em: Dict[Tuple[int, int], List[Tuple[int, int, int]]] = {}
            for t, tri in enumerate(self.tris):
                v = (tri[0], tri[1], tri[2], tri[0])
                for s in range(3):
                    a, b = v[s], v[s + 1]
                    em.setdefault(ekey(a, b), []).append((t, a, b))
            return em

        oriented = [False] * n
        for seed in range(n):
            if oriented[seed]:
                continue
            q: deque = deque([seed])
            oriented[seed] = True
            while q:
                t = q.popleft()
                em = edge_map()
                tri = self.tris[t]
                vv = (tri[0], tri[1], tri[2], tri[0])
                for s in range(3):
                    a, b = vv[s], vv[s + 1]
                    uses = em.get(ekey(a, b), [])
                    if len(uses) != 2:
                        continue
                    for u, ua, ub in uses:
                        if u == t or oriented[u]:
                            continue
                        if (ua, ub) == (a, b):
                            x, y, z = self.tris[u]
                            self.tris[u] = (x, z, y)
                        oriented[u] = True
                        q.append(u)


def count_manifold_components(comp: Component) -> int:
    """Harness loadMesh split: unite triangles across edges used exactly twice."""
    n = len(comp.tris)
    if n == 0:
        return 0
    parent = list(range(n))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def unite(a: int, b: int) -> None:
        a, b = find(a), find(b)
        if a != b:
            parent[a] = b

    e2t: Dict[int, List[int]] = {}
    for t in range(n):
        for s in range(3):
            eid = comp.tri_edges[t][s]
            e2t.setdefault(eid, []).append(t)
    for ts in e2t.values():
        if len(ts) == 2:
            unite(ts[0], ts[1])
    return len({find(i) for i in range(n)})


def poly_area2(poly: Sequence[Tuple[float, float]]) -> float:
    n = len(poly)
    a = 0.0
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        a += x1 * y2 - x2 * y1
    return 0.5 * a


def earclip_xy(poly: Sequence[Tuple[float, float]]) -> List[Tuple[int, int, int]]:
    """Ear-clip a simple polygon in 2D. `poly` CCW or CW. Returns index triples."""
    n0 = len(poly)
    if n0 < 3:
        return []
    idx = list(range(n0))
    area = poly_area2(poly)
    if area < 0:
        idx.reverse()
        area = -area

    def at(i: int) -> Tuple[float, float]:
        return poly[idx[i]]

    def is_ear(i: int) -> bool:
        n = len(idx)
        a, b, c = at((i - 1) % n), at(i), at((i + 1) % n)
        cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        if cross <= 1e-18:
            return False
        for j in range(n):
            if j in ((i - 1) % n, i, (i + 1) % n):
                continue
            p = at(j)
            def sign(u, v, w):
                return (v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0])
            s1, s2, s3 = sign(a, b, p), sign(b, c, p), sign(c, a, p)
            if (s1 >= -1e-18 and s2 >= -1e-18 and s3 >= -1e-18) or (
                s1 <= 1e-18 and s2 <= 1e-18 and s3 <= 1e-18
            ):
                if abs(s1) + abs(s2) + abs(s3) > 1e-18:
                    return False
        return True

    tris: List[Tuple[int, int, int]] = []
    guard = 0
    while len(idx) > 3 and guard < n0 * n0:
        guard += 1
        n = len(idx)
        found = False
        for i in range(n):
            if is_ear(i):
                tris.append((idx[(i - 1) % n], idx[i], idx[(i + 1) % n]))
                del idx[i]
                found = True
                break
        if not found:
            break
    if len(idx) == 3:
        tris.append((idx[0], idx[1], idx[2]))
    return tris


def _parse_owner(o: Tuple[str, int]) -> Tuple[int, int]:
    if o[0] == "r":
        return o[1], -1
    if o[0] == "i":
        return -1, o[1]
    return -1, -1


def extract_chains(
    comp: Component,
    tri_region: Sequence[int],
    tri_island: Sequence[int],
) -> List[dict]:
    """Boundary chains with I2/I8 splits. Owner-pair walks; split at deg!=2 or >=3 owners."""
    n_tri = len(comp.tris)

    def own(t: int) -> Tuple[str, int]:
        if tri_region[t] >= 0:
            return ("r", tri_region[t])
        if tri_island[t] >= 0:
            return ("i", tri_island[t])
        return ("x", -1)

    uses: Dict[int, List[Tuple[int, int, int, Tuple[str, int]]]] = {}
    for t in range(n_tri):
        tri = comp.tris[t]
        te = comp.tri_edges[t]
        v = (tri[0], tri[1], tri[2], tri[0])
        for s in range(3):
            eid = te[s]
            uses.setdefault(eid, []).append((t, v[s], v[s + 1], own(t)))

    boundary: Dict[int, List[Tuple[int, int, int, Tuple[str, int]]]] = {}
    nm_verts: Set[int] = set()
    for eid, u in uses.items():
        owners = {rec[3] for rec in u}
        if len(u) > 2:
            a, b = comp.comp_edges[eid]
            nm_verts.add(a)
            nm_verts.add(b)
        if len(u) == 1:
            boundary[eid] = u
        elif len(owners) >= 2:
            boundary[eid] = u

    v_owners: Dict[int, Set[Tuple[str, int]]] = {}
    for eid, u in boundary.items():
        ow = {rec[3] for rec in u}
        if len(u) == 1:
            ow.add(("x", -1))
        a, b = comp.comp_edges[eid]
        v_owners.setdefault(a, set()).update(ow)
        v_owners.setdefault(b, set()).update(ow)
    split_verts = {v for v, ow in v_owners.items() if len(ow) >= 3}
    split_verts |= nm_verts

    def pair_of(eid: int) -> frozenset:
        u = boundary[eid]
        ows = {rec[3] for rec in u}
        if len(u) == 1:
            ows.add(("x", -1))
        return frozenset(ows)

    groups: Dict[frozenset, List[int]] = {}
    for eid in boundary:
        groups.setdefault(pair_of(eid), []).append(eid)

    chains: List[dict] = []
    for pair, eids in groups.items():
        adj: Dict[int, List[Tuple[int, int]]] = {}
        for eid in eids:
            va, vb = comp.comp_edges[eid]
            adj.setdefault(va, []).append((eid, vb))
            adj.setdefault(vb, []).append((eid, va))
        terminals = {v for v, nbrs in adj.items() if len(nbrs) != 2 or v in split_verts}
        used: Set[int] = set()

        def walk(start_v: int, start_eid: int) -> Tuple[List[int], List[int], bool]:
            va, vb = comp.comp_edges[start_eid]
            cur = vb if va == start_v else va
            edges = [start_eid]
            verts = [start_v, cur]
            used.add(start_eid)
            while True:
                if cur == start_v and len(edges) >= 2:
                    verts.pop()
                    return edges, verts, True
                if cur in terminals and cur != start_v:
                    return edges, verts, False
                nxts = [(eid, oth) for eid, oth in adj[cur] if eid not in used]
                if not nxts:
                    return edges, verts, False
                eid, oth = min(nxts, key=lambda x: x[0])
                used.add(eid)
                edges.append(eid)
                cur = oth
                verts.append(cur)

        for v in sorted(terminals):
            for eid, _oth in sorted(adj[v], key=lambda x: x[0]):
                if eid in used:
                    continue
                edges, verts, closed = walk(v, eid)
                chains.append(_orient_chain(comp, uses, pair, edges, verts, closed))
        for eid in sorted(eids):
            if eid in used:
                continue
            va, vb = comp.comp_edges[eid]
            edges, verts, closed = walk(va, eid)
            chains.append(_orient_chain(comp, uses, pair, edges, verts, closed))

    chains.sort(key=lambda c: min(c["meshEdges"]) if c["meshEdges"] else 0)
    return chains


def _orient_chain(
    comp: Component,
    uses: Dict[int, List[Tuple[int, int, int, Tuple[str, int]]]],
    pair: frozenset,
    edges: List[int],
    verts: List[int],
    closed: bool,
) -> dict:
    if not edges:
        return {
            "regA": -1, "regB": -1, "islandA": -1, "islandB": -1,
            "tangent": False, "closedLoop": closed, "meshEdges": [], "meshVerts": verts,
        }
    eid0 = edges[0]
    a, b = verts[0], verts[1] if len(verts) > 1 else verts[0]
    left: Optional[Tuple[str, int]] = None
    for _t, fr, to, ow in uses[eid0]:
        if fr == a and to == b:
            left = ow
            break
    if left is None:
        edges = list(reversed(edges))
        if closed:
            verts = [verts[0]] + list(reversed(verts[1:]))
        else:
            verts = list(reversed(verts))
        a, b = verts[0], verts[1] if len(verts) > 1 else verts[0]
        for _t, fr, to, ow in uses[eid0]:
            if fr == a and to == b:
                left = ow
                break
    if left is None:
        left = next(iter(pair))
    others = [o for o in pair if o != left]
    right = others[0] if others else ("x", -1)
    ra, ia = _parse_owner(left)
    rb, ib = _parse_owner(right)
    return {
        "regA": ra,
        "regB": rb,
        "islandA": ia,
        "islandB": ib,
        "tangent": False,
        "closedLoop": closed,
        "meshEdges": edges,
        "meshVerts": verts,
    }


def chain_g1(chain: dict, regions: Sequence[dict], eps_plane: float = 0.02) -> bool:
    """refit.hpp BoundaryChain::tangent geometric test."""
    ra, rb = chain["regA"], chain["regB"]
    if ra < 0 or rb < 0:
        return False
    a, b = regions[ra], regions[rb]
    sin3 = math.sin(math.radians(3.0))
    cos3 = math.cos(math.radians(3.0))

    def is_pl(r: dict) -> bool:
        return r["type"] == "plane"

    def is_cy(r: dict) -> bool:
        return r["type"] == "cylinder"

    def ax(r: dict) -> Tuple[Vec3, Vec3]:
        return tuple(r["ax"]["loc"]), tuple(r["ax"]["dir"])  # type: ignore

    if is_pl(a) and is_cy(b):
        pl, cy = a, b
    elif is_pl(b) and is_cy(a):
        pl, cy = b, a
    elif is_cy(a) and is_cy(b):
        loc1, d1 = ax(a)
        loc2, d2 = ax(b)
        d1u, d2u = vunit(d1), vunit(d2)
        if abs(vdot(d1u, d2u)) < cos3:
            return False
        w = vsub(loc2, loc1)
        cross = vcross(d1u, w)
        dist = vnorm(cross)
        r1, r2 = a["radius"], b["radius"]
        return min(abs(dist - abs(r1 + r2)), abs(dist - abs(r1 - r2))) <= eps_plane
    elif is_pl(a) and is_pl(b):
        n1 = vunit(tuple(a["ax"]["dir"]))
        n2 = vunit(tuple(b["ax"]["dir"]))
        return abs(vdot(n1, n2)) >= cos3
    else:
        return False

    loc_p, n = ax(pl)
    loc_c, d = ax(cy)
    n = vunit(n)
    d = vunit(d)
    if abs(vdot(d, n)) > sin3:
        return False
    dist = abs(vdot(vsub(loc_c, loc_p), n))
    return abs(dist - cy["radius"]) <= eps_plane


def region_vertices(comp: Component, tris: Sequence[int]) -> List[int]:
    seen: Set[int] = set()
    out: List[int] = []
    for t in tris:
        for v in comp.tris[t]:
            if v not in seen:
                seen.add(v)
                out.append(v)
    return out


def lowest_id_xdir(pts: Sequence[Vec3], verts: Sequence[int], loc: Vec3, axis: Vec3) -> Vec3:
    if not verts:
        return (1.0, 0.0, 0.0)
    vid = min(verts)
    radial = vsub(pts[vid], loc)
    axial = vunit(axis)
    radial = vsub(radial, vscale(axial, vdot(radial, axial)))
    u = vunit(radial)
    if vnorm(u) < 1e-12:
        # pick a perpendicular
        if abs(axial[0]) < 0.9:
            u = vunit(vcross(axial, (1.0, 0.0, 0.0)))
        else:
            u = vunit(vcross(axial, (0.0, 1.0, 0.0)))
    return u
