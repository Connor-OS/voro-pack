import pyvoro as voro
import trimesh
from scipy import ndimage
import sys
import math
import random
import time
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

##TODO
#add buffer to bounding box
#raise more exceptions if part of the code fails due to inputs
#determine the difference between input PSD and final particle PSD 
#rearange particles to improve packing density??
#split a few more operations into individual functions

if(len(sys.argv) != 4):
    print("Error: Input arguments are of the form: \nvor_pack.py input_geom input_PSD output_file\n PLease check inputs and try again")
    exit()

if(sys.argv[1][-4:].lower() != '.stl'):
    print('Error: Input geometry must be an stl file')
    exit()

#Load mesh and calculate bounding box
mesh = trimesh.load(sys.argv[1])
input_PSD = sys.argv[2]
output_file = sys.argv[3]

bb = mesh.bounds
xmin,ymin,zmin = bb[0]
xmax,ymax,zmax = bb[1]
width, depth, height = bb[1]-bb[0]

#user defined size distribution Size distribution
size_dist = pd.read_csv(f'./{input_PSD}',delimiter='\t', names=['diam','vol']) #users sould uses vol percent but scaling will be the same
size_dist['diam'] = size_dist['diam']/1000 #convert from m to μm
size_dist['rad'] = size_dist['diam']/2
size_dist['freq'] = size_dist['vol']/(size_dist['vol'] * (size_dist['diam'])**3)
size_dist['freq'] = size_dist['freq']/(size_dist['freq'].sum())
size_dist['freq'] = size_dist['freq'].fillna(0.)
cfd = [size_dist['freq'][0]]
for i,v in enumerate(size_dist['freq'][1:]):
    cfd.append(v+cfd[i])
size_dist['cfd'] = cfd
print('User defined PSD\n')
print(size_dist)
dmin = size_dist['diam'].min()
dmax = size_dist['diam'].max()

#=============================================
def main():
    start = time.perf_counter()
    ##DISK SAMPLING
    print('disk sampling...')
    tic = time.perf_counter()
    points = disk_sampling()
    toc = time.perf_counter()
    print('COMPLETE\n')
    print(f'time elapsed {toc - tic}\n')

    ##BUILD VORONOI
    print('building voronoi...')
    tic = time.perf_counter()
    #use voro++ to create a weighted voronoi diagram from the point sampling positions
    vor = voro.compute_voronoi(
        points[:,:3],
        [[xmin-0.1,xmax+0.1], [ymin-0.1,ymax+0.1], [zmin-0.1,zmax+0.1]], # limits
        2.0, # block size
        periodic=[True]*3,
        radii=points[:,3]
    )
    toc = time.perf_counter()
    print('COMPLETE')
    print(f'time elapsed {toc - tic}\n')

    ##PARTICLE PLACEMENT
    print('Placing particles...')
    tic = time.perf_counter()
    particles = []
    for i,region in enumerate(vor):
        print(f'placing particle {i}/{len(vor)}',end='\r')
        particles.append(circ_inscribed([face['vertices'] for face in region['faces']],np.array(region['vertices'])))
    particles = np.array([i for i in particles if i[0]])
    toc = time.perf_counter()
    print('COMPLETE')
    print(f'time elapsed {toc - tic}\n')
    
    ##TRIM TO MESH
    print('Trimming to mesh...')
    tic = time.perf_counter()
    exists = []
    exists_ = []
    chunk = 500
    i=1
    while chunk*i <len(particles):
        exists += mesh.contains(particles[chunk*(i-1):chunk*i,:3]).tolist()
        exists_ += mesh.contains(points[chunk*(i-1):chunk*i,:3]).tolist()
        i+=1
    exists += mesh.contains(particles[chunk*(i-1):len(particles),:3]).tolist()
    exists_ += mesh.contains(points[chunk*(i-1):len(particles),:3]).tolist()
    particles = [v for i,v in enumerate(particles) if exists[i] == True]
    points = [v for i,v in enumerate(points) if exists_[i] == True]
    toc = time.perf_counter()
    print('COMPLETE')
    print(f'time elapsed {toc - tic}\n')

    ##OUTPUTS
    print('Packing COMPLETE, Outputing results..\n')
    stop = time.perf_counter()
    print(f'total time take {stop - start}\n')
    output(particles,points,output_file)


#=============================================
def generate_radius():
    rnd =  np.random.random()
    i = 0
    while rnd > size_dist['cfd'][i] and i<len(size_dist['cfd']): i+=1
    return (size_dist['rad'][i]-((size_dist['cfd'][i]-rnd)*(size_dist['rad'][i]-size_dist['rad'][i-1]))/(size_dist['cfd'][i]-size_dist['cfd'][i-1]))


#=============================================
def disk_sampling():
    #(Bridson,2007) Fast Poisson Disk Sampling in Arbitrary Dimensions
    #constants
    n = 3 #Dimensions
    attempts = 30 
    #initialise grid
    w = dmin
    search_range = math.ceil(dmax/dmin)
    cols = int(width/w)
    rows = int(depth/w)
    plains = int(height/w)
    grid = np.array([[None,None,None,None]]*cols*rows*plains) #empty grid to stor points and radii
    active = []
    n_part = 1

    #step 1: pick a random starting point
    pos = np.concatenate((np.random.rand(n)*[width,depth,height], [generate_radius()]), axis=0)
    col = int(pos[0]/w)
    row = int(pos[1]/w)
    plain = int(pos[2]/w)
    grid[col+row*cols+plain*rows*cols] = pos
    active.append(pos)

    #step 2: Poission disc sampling algorythm
    r_new = generate_radius()
    while len(active) > 0:
        randID = random.randint(0,len(active)-1)
        pos = active[randID]
        found = False
        for i in range(attempts):
            valid_sample = True
            a1 = random.random()*2*np.pi
            a2 = random.random()*2*np.pi
            r = pos[3]+r_new
            sample = pos +[np.sin(a1)*r,np.cos(a1)*np.cos(a2)*r,np.sin(a2)*r,r_new-pos[3]]
            col =int(sample[0]/w)
            row = int(sample[1]/w)
            plain = int(sample[2]/w)
            #check sample validity
            if(col>-1 and row>-1 and plain>-1 and col<cols and row<rows and plain<plains and not grid[col+row*cols+plain*rows*cols][0]):
                for i in range(-search_range,search_range+1):
                    for j in range(-search_range,search_range+1):
                        for k in range(-search_range,search_range+1):
                            if((col+i)>-1 and (row+j)>-1 and (plain+k)>-1 and (col+i)<cols and (row+j)<rows and (plain+k)<plains):
                                neighbour = grid[(col+i)+(row+j)*cols+(plain+k)*rows*cols]    
                                if(neighbour[0] and ((sample[0]-neighbour[0])**2+(sample[1]-neighbour[1])**2+(sample[2]-neighbour[2])**2)**0.5 < neighbour[3]+r_new):
                                    valid_sample = False
                #accept sample and amend active list
                if(valid_sample == True):
                    found = True
                    n_part += 1
                    grid[col+row*cols+plain*rows*cols] = sample
                    active.append(sample)
                    r_new = generate_radius()
                    break
        if(not found):
            active.pop(randID)
    grid = np.array([i+[xmin,ymin,zmin,0] for i in grid if i[0]])
    print(f'sampled {len(grid)} points')
    return grid


#=============================================
def circ_inscribed(faces,vertices):
    #this function finds the largest inscribed circle inside voronoi region and returns it as a particle
    GRID_DIVS = 25 #This constant has a massive affect on performance, it needs to be tuned, increase in this constant will increase the prcesion of the inscribed circle
    centroid = np.array([0.0,0.0,0.0])
    for v in vertices:
        centroid += v
    centroid /= len(vertices)
    xr,yr,zr = max(vertices[:,0])-min(vertices[:,0]),max(vertices[:,1])-min(vertices[:,1]),max(vertices[:,2])-min(vertices[:,2])
    div = min(xr,yr,zr)/GRID_DIVS #divide cell into a uniform grid for distance transform later
    #construct grid
    x_axis = np.arange(min(vertices[:,0]-div), max(vertices[:,0]+div), div)
    y_axis = np.arange(min(vertices[:,1]-div), max(vertices[:,1]+div), div)
    z_axis = np.arange(min(vertices[:,2]-div), max(vertices[:,2]+div), div)
    grid_points = [[x,y,z] for x in x_axis for y in y_axis for z in z_axis]

    nx,ny,nz = len(x_axis),len(y_axis),len(z_axis)
    grid = np.ones([nx,ny,nz]) #ones mask
    #apply mask for each face to remove points outside cell
    for face in faces:
        #find equation of plane
        p0, p1, p2 = np.array(vertices[face[0]]), np.array(vertices[face[1]]), np.array(vertices[face[2]])
        v1,v2 = p2-p0, p1-p0
        cp = np.cross(v1,v2)
        a,b,c = cp
        d = np.dot(cp,p2)
        #flip plane if contains centroid
        if(a*centroid[0] + b*centroid[1] + c*centroid[2] < d):
            a, b, c, d = -a,-b, -c, -d 
        is_outside = lambda p: a*p[0] + b*p[1] + c*p[2] < d
        grid = grid - np.array([is_outside(i) for i in grid_points]).reshape([nx,ny,nz])

    grid = grid.clip(min=0) #grid points inside the cell are 1 outside the cell are 0

    dt_image = ndimage.distance_transform_edt(grid) #create euclidian distance transform matrix

    #initalise particle with max inscribed radius
    r = max(dt_image.flatten())
    cent_inscribed = np.where(dt_image == r)
    r = r*div-div/2 #convert to the simulation space and take one div from diam to ensure no particle overlap

    cent_inscribed = list(zip(cent_inscribed[0],cent_inscribed[1],cent_inscribed[2]))[0]
    inscribed_circle = [x_axis[cent_inscribed[0]], y_axis[cent_inscribed[1]], z_axis[cent_inscribed[2]], r]
    return inscribed_circle


#=============================================
def output(particles,points,output):
    # Write output files
    with open(output,'w') as file:
        # file.writelines([f'{len(particles)}\n','\n']) #print number of particles first for LAMMPS Style output
        for i,v in enumerate(particles):
            file.write(f'{v[0]} {v[1]} {v[2]} {v[3]}\n')

    #uncomment for points output    
    # with open('points.txt','w') as file:
    #     # file.writelines([f'{len(particles)}\n','\n']) #print number of particles first for LAMMPS Style output
    #     for i in points:
    #         file.write(f'{i[0]} {i[1]} {i[2]} {i[3]}\n')

    data2 = pd.DataFrame(points)
    data2.columns = ['x', 'y', 'z','r']
    data2['r'].plot.hist(bins=40)

    #Show size dist
    data = pd.DataFrame(particles)
    data.columns = ['x', 'y', 'z','r']
    data['r'].plot.hist(bins=40)

    plt.legend(['points','particles'])
    plt.show()


if __name__ == '__main__':
    main()