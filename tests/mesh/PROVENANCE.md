# seamed_cylinder.step

Plain right circular cylinder authored with OpenCASCADE for `mesh_edges_seam`.

| Parameter | Value |
|---|---|
| Radius | 10 mm |
| Height | 30 mm |
| Generator | `BRepPrimAPI_MakeCylinder(10.0, 30.0)` + `STEPControl_Writer` |
| Source | `gen_seamed_cylinder.cpp` (one-shot; not built by default) |

## sha256

```
e1fcd16158f8d2a787b96e72acdf22b2f8321d23c5ad97873cb80394d94c05ed
```

## Expected drawable edge topology

After `--mesh --edges`:

- **3 drawable edges** in `MESH_RESULT.edges`: 2 circular cap edges + 1 axial seam
- **91 segments** in the Format A buffer (45 tessellated circle segments per cap,
  1 straight seam segment) — far below facet-boundary count (528 for 176 triangles)
