# PotConstTri3D
This is an implementation of a Laplace equation BEM solver to use alongside FreeCAD and Gmsh. Some examples are solved.

This program is written in Julia 0.6.2.

To run the solver, simply write in the Julia shell:

`include("PotConstTri3D")`

Be sure that you're in the same directory of the program. 

If you want to run a simulation of your own model, change the dad.jl file and put the name of the mesh file and the boundary conditions on the faces.
