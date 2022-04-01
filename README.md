# meshPython2Foam

Write OpenFOAM blockMeshDict from python

Inspired by pypi.org/project/ofblockmeshdicthelper/

Examples illustrate how to use this utility. The underlying idea is
that blocks in blockMesh become cells in the mesh, i.e. we write
every cell from python as a block and blockMesh then takes care
of the indexing and so on. Of course blocks may still be split into
many cells as originally envisioned by blockMesh.

The strength of this utility is that it takes care of doublicate vertices. 
In python, we add a cell (a block actually) and it's neighbor, which uses 
some of the same vertices, but in the final blockMeshDict, the overlapping 
vertices only show up once.

If we create very many blocks, blockMesh will take a long time
(many hours) to generate the mesh as it can not run in parallel.
A mesh with 860850 vertices took 3.5h to run blockMesh (used 2.8GB of RAM)
A mesh with 5.9e6 vertices did not finish in 19h (used 17GB of RAM so far)
