# VTK-Output
## The UnstructuredGrid format (.vtu)
The general format for UnstructuredGrid is:
```
<VTKFile type=”UnstructuredGrid” ...>
  <UnstructuredGrid>
    <Piece NumberOfPoints=”#” NumberOfCells=”#”>
      <PointData>...</PointData>
      <CellData>...</CellData>
      <Points>...</Points>
      <Cells>...</Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
```


The following composes the .vtu files for convenient read. It is supposed to have .pvd ending.
```
<VTKFile type="Collection">
<Collection>
<DataSet timestep=".00000" file="./vtk/particles_000000.000000.vtu"/>
<DataSet timestep="16.00000" file="./vtk/particles_000001.000000.vtu"/>
<DataSet timestep="32.00000" file="./vtk/particles_000002.000000.vtu"/>
<DataSet timestep="48.00000" file="./vtk/particles_000003.000000.vtu"/>
<DataSet timestep="64.00000" file="./vtk/particles_000004.000000.vtu"/>
<DataSet timestep="80.00000" file="./vtk/particles_000005.000000.vtu"/>
<DataSet timestep="96.00000" file="./vtk/particles_000006.000000.vtu"/>
<DataSet timestep="112.00000" file="./vtk/particles_000007.000000.vtu"/>
<DataSet timestep="128.00000" file="./vtk/particles_000008.000000.vtu"/>
<DataSet timestep="144.00000" file="./vtk/particles_000009.000000.vtu"/>
<DataSet timestep="160.00000" file="./vtk/particles_000010.000000.vtu"/>
</Collection>
</VTKFile>
```
## Implementation:

All Paraview/VTU output functionality is implemented within the
MPMOutputVTK class.
