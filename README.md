# vor-pack  !In Progress!
fast particle packing algorythim for arbitrary shaped regions using voro++ and trimesh.

to run the code use

``` python voro_pack.py input_geom.stl PSD.csv output_file.txt ```

included are a number of example stl files and potential particle size distributions to test.

Known bugs
code may intermittently fail in the disk sampling function when trying to place paricle outside of the initlaised grid. This happens approx 1/20 times you run the code. I plan to completely overhall this function in the future.

Voronoi digram might intermittently fail if useing an stl much larger than the ones included in the repository. This is because I was too lazy to make this function generalisable. Will fix this in the near future.
