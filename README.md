# vor-pack - A fast particle packing algorythim for arbitrary shaped regions.

This is a program for packing spheres into arbitrary user defined regions, controling the distribution of sphere radii. It is inteded to be used as a pre processing procedure for generating clumps of bonded particles that can be used for descrete element modeling for granular materials.

The code takes advantage of the Voronoi diagram as well as Lloyds algorythm for voronoi cell relaxation to determine particle positions. One of the advantages of this method is that it provides a fast huristic solution to the sphere packing problem without any dynamic simulation, therfore avoiding unnecesary computational overhead.

It is capable of constructing clumps of highly polydisperse particles, where the largest particle is up to ten times larger than the smallest 

# Dependancies
The code is writen in python and relies on several third party depedancies.

To istall all dependacies run the following in your terminal
'''pip install pyvoro
pip install trimesh
pip install cvxpy
pip install numpy
pip install pandas
pip install matplotlib'''

or equivelent in your chosen package manager

#Input
the code takes two input files:

- A watertight mesh in stl file format
- A tab seperated PSD distribution in csv or txt file format representing sived particle data. Format: two tab seperated columns: Diamater and volume/mass percentage passing

#Output
The code writes a single output file specifying particle possitions and radii in the format: x y z r

This format is lightweight and can be recognised by many popular DEM and visualisation software packages.

#Running the code

to run the code use:

``` python voro_pack.py input_geom input_PSD output_file```

input_geom
-This must be the path to your mesh.

input_PSD
-This must be the path to your particle size distribution.

output_file
-This must be the path to and name you choose for your output file.

#To be added in the future
- 2D Jupyter notebook version of the code for explanation
- some test PSD's and example clumps
