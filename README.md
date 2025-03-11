# vor-pack - A fast particle packing algorithm for arbitrary shaped regions.

This is a program for packing spheres into arbitrary user defined regions, controlling the distribution of sphere radii. It is intended to be used as a pre-processing procedure for generating clumps of bonded particles that can be used for discrete element modeling for granular materials.

The code takes advantage of the Voronoi diagram as well as Lloyds algorythm for voronoi cell relaxation to determine particle positions. One of the advantages of this method is that it provides a fast heuristic solution to the sphere packing problem without any dynamic simulation, therefore avoiding unnecessary computational overhead.

It is capable of constructing clumps of highly polydisperse particles, where the largest particle is up to ten times larger than the smallest 

## Dependencies
The code is writen in python and relies on several third party dependencies.

To install all dependencies run the following in your terminal
```
pip install pip install -r requirements.txt 
```

## Input
the code takes two input files:

- A watertight mesh in stl file format
- A tab seperated PSD distribution in csv or txt file format representing sived particle data. Format: two tab seperated columns: Diameter and volume/mass percentage passing

## Output
The code writes a single output file specifying particle positions and radii in the format: x y z r

This format is lightweight and can be recognised by many popular DEM and visualisation software packages.

## Running the code

to run the code use:

``` python voro_pack.py input_geom input_PSD output_file```

input_geom
-This must be the path to your mesh.

input_PSD
-This must be the path to your particle size distribution.

output_file
-This must be the path to and name you choose for your output file.

### Acknowledging Voro-pack
The voro-pack code can be cited with this jornal article, which provides a description and a number of illustrative examples

https://doi.org/10.1016/j.partic.2025.01.014
Chao Zhang, Connor O'Shaughnessy, Sadaf Maramizonouz, Vasileios Angelidakis, Sadegh Nadimi,
Controlling fragment size distribution for modelling the breakage of multi-sphere particles,
Particuology