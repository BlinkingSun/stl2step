// One-shot generator for tests/mesh/seamed_cylinder.step (R=10 mm, H=30 mm).
// SPDX-License-Identifier: MIT

#include <BRepPrimAPI_MakeCylinder.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Writer.hxx>

#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <out.step>\n", argv[0]);
        return 1;
    }
    BRepPrimAPI_MakeCylinder cyl(10.0, 30.0);
    STEPControl_Writer w;
    if (w.Transfer(cyl.Shape(), STEPControl_AsIs) != IFSelect_RetDone) return 1;
    if (w.Write(argv[1]) != IFSelect_RetDone) return 1;
    return 0;
}
