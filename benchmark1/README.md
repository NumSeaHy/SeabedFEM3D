# Solving 3-D pulsant sphere problem by using a Perfectly Matched Layer inside Gridap.

The problem is fully parameterised, so use the `Configuration.jl` file to introduce the desire parameters.

The `Mesh.jl` file has inside a `FunMesh` function the let generate meshes with different number of elements `N` per wavelength and then save the meshes into the `data` folder. 

The `Run.jl` file has inside the `RunFEM` function that takes as input the mesh name and run the finite element simulation by using Gridap. Inside the file a convergence test cis performed by using the meshes generated with the `FunMesh` function inside the `Mesh.jl` file.


### How to run
```julia
include("src/Run.jl")
```

### Authors
This work has been carried out by Andres Prieto Aneiros (andres.prieto@udc.es) and Pablo Rubial Yáñez (p.rubialy@udc.es) during the work developed in the [NumSeaHy](https://dm.udc.es/m2nica/en/node/157) project.

### Interesting links

Editing the mesh size of a 3D mesh build with the OCC kernel instead of the Geo one.

[GmshPythonExample](https://jsdokken.com/src/tutorial_gmsh.html)

Generating a volume mesh from a CAD stl file:

[GmshCAD1](https://project.inria.fr/softrobot/documentation/from-design-to-mesh-generation-using-freecad-and-gmsh/)

[GmshCAD2](https://www.pygimli.org/_examples_auto/1_meshing/plot_cad_tutorial.html)

### License
 <p xmlns:cc="http://creativecommons.org/ns#" >This work is licensed under <a href="http://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1"></a></p> 

