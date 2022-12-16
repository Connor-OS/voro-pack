# vor-pack - A fast particle packing algorythim for arbitrary shaped regions.

to run the code use:

``` python voro_pack.py input_geom input_PSD output_file```

input_geom
-This must be the path to a watertight mesh in stl file format

input_PSD
-This must be the path to a tab seperated PSD distribution in csv or txt file format representing sived particle data.

format: two tab seperated columns: Diamater and volume/mass percentage passing

output_file
-This must be the path to the file where you want to write the final particle postions to

format: x y z r

To be added:
- 2D Jupyter notebook with examples
- some test PSD
