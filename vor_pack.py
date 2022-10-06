import numpy as np
from matplotlib import pyplot as plt
import random
import pandas as pd
import time
import trimesh 
from scipy.spatial import Voronoi
from scipy import ndimage
import sys

if(len(sys.argv) != 2):
    print("Error: Input arguments are of the form: input geometry, e.g:\nvor_pack.py geom.stl\nPLease check inputs and try again")
    exit()

if(sys.argv[1][-3:] != '.stl'):
    print('Error: Input geometry must be an stl file')
    exit()

#Load mesh and calculate bounding box
mesh = trimesh.load(sys.argv[1])
bb = mesh.bounds
xmin,ymin,zmin = bb[0]
xmax,ymax,zmax = bb[1]
width, depth, height = bb[1]-bb[0]

#user defined size distribution Size distribution
# size_dist = pd.read_csv('./CrushedSandPSD.csv')

#Parameters for normal distribution
rmean = 0.1
rsd = 0.05

##TODO
#Add user defined particle size distribution
#Configure min_dist and div
#Add buffer to bounding box
#Handle extreme cases, very large diffrence between largest, to smallest particles
#Raise more exceptions if part of the code fails due to inputs

#=============================================
def main():
    # creat bounding box around mesh ... do disk sampling inside this bounding box (or mesh with buffer?)
    print('disk sampling...')
    tic = time.perf_counter()
    points = disk_sampling()
    toc = time.perf_counter()
    print('COMPLETE\n')
    print(f'time elapsed {toc - tic}\n')


    print('building voronoi...')
    tic = time.perf_counter()
    vor  = Voronoi(points)
    toc = time.perf_counter()
    print('COMPLETE')
    print(f'time elapsed {toc - tic}\n')


    # print('annoying unpack...')
    # tic = time.perf_counter()
    vertices = vor.vertices
    #unpack finite regions
    regions = [i for i in vor.regions if i and -1 not in i]

    faces = []
    for i,v in enumerate(regions):
        faces.append([w for w in vor.ridge_vertices if set(w) <= set(v)]) #Becomes very costly with highnumbers of particle and is 
    toc = time.perf_counter()
    # print('COMPLETE')
    # print(f'time elapsed {toc - tic}\n')

    print('Placing particles...')
    particles = []
    tic = time.perf_counter()
    for i,region in enumerate(regions):
        print(f'placing particle {i}/{len(regions)}',end='\r')
        particles.append(circ_inscribed(region,faces[i],vertices))

    particles = np.array([i for i in particles if i[0]])
    toc = time.perf_counter()
    print('COMPLETE')
    print(f'time elapsed {toc - tic}\n')
    
    #Trim to stl mesh
    print('Trimming to mesh...')
    tic = time.perf_counter()

    exists = []
    chunk = 500
    i=1
    while chunk*i <len(particles):
        exists += mesh.contains(particles[chunk*(i-1):chunk*i,:3]).tolist()
        i+=1
    exists += mesh.contains(particles[chunk*(i-1):len(particles),:3]).tolist()
    particles = [v for i,v in enumerate(particles) if exists[i] == True]
    toc = time.perf_counter()
    print('COMPLETE')
    print(f'time elapsed {toc - tic}\n')

    print('Packing COMPLETE, Outputing results..')
    output(points,vertices,particles)


#=============================================
def disk_sampling():
    #(Bridson,2007) Fast Poisson Disk Sampling in Arbitrary Dimensions
    #constants
    n = 3
    attempts = 30
    #initialise grid
    w = rmean/n**0.5
    cols = int(width/w)
    rows = int(depth/w)
    plains = int(height/w)
    grid = np.array([[None,None,None,None]]*cols*rows*plains)
    active = []
    n_part = 1

    #step 1: pick a random starting point
    pos = np.random.rand(n+1)*[width,depth,height,rmean]
    col = int(pos[0]/w)
    row = int(pos[1]/w)
    plain = int(pos[2]/w)
    grid[col+row*cols+plain*rows*cols] = pos
    active.append(pos)

    r_new =  np.random.normal(rmean,rsd) #replace with sampling from user generated distribution

    #step 2: Poission disc sampling algorythm
    while len(active) > 0:
        randID = random.randint(0,len(active)-1)
        pos = active[randID]
        found = False
        for i in range(attempts):
            valid_sample = True
            a1 = random.random()*2*np.pi
            a2 = random.random()*2*np.pi
            r = pos[3]+r_new
            # r = random.random()*(pos[3]+r_new)+(pos[3]+r_new) #old method, adds unecesary space to packing
            sample = pos +[np.sin(a1)*r,np.cos(a1)*np.cos(a2)*r,np.sin(a2)*r,r_new-pos[3]]
            col =int(sample[0]/w)
            row = int(sample[1]/w)
            plain = int(sample[2]/w)
            #check sample validity
            if(col>-1 and row>-1 and plain>-1 and col<cols and row<rows and plain<plains and not grid[col+row*cols+plain*rows*cols][0]):
                for i in range(-2,3):
                    for j in range(-2,3):
                        for k in range(-2,3):
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
                    r_new =  np.random.normal(rmean,rsd)
                    break
        if(not found):
            active.pop(randID)
    grid = np.array([i[:3]+[xmin,ymin,zmin] for i in grid if i[0]])

    return grid


#=============================================
def circ_inscribed(region,faces,vertices):
    #find define region and initalise grid of points over it
    GRID_DIVS = 15 #This constant has a massive affect on performance, it needs to be tuned, any increase in this constant will increase the prcesion of the inscribed circle
    reg_verts = []
    centroid = np.array([0.0,0.0,0.0])
    for i,v in enumerate(vertices):
        if i in region:
            reg_verts.append(v)
            centroid += v
    reg_verts = np.array(reg_verts)
    centroid /= len(reg_verts)
    #break out if centroi outside bounding box
    if centroid[0] < xmin or centroid[0] > xmax or centroid[1] < ymin or centroid[1] > ymax or centroid[2] < zmin or centroid[2] > zmax:
        return [None,None,None,None]
    xr,yr,zr = max(reg_verts[:,0])-min(reg_verts[:,0]),max(reg_verts[:,1])-min(reg_verts[:,1]),max(reg_verts[:,2])-min(reg_verts[:,2])
    div = min(xr,yr,zr)/GRID_DIVS
    #point grid
    x = np.arange(min(reg_verts[:,0]), max(reg_verts[:,0]), div)
    y = np.arange(min(reg_verts[:,1]), max(reg_verts[:,1]), div)
    z = np.arange(min(reg_verts[:,2]), max(reg_verts[:,2]), div)
    nx,ny,nz = len(x),len(y),len(z)
    grid = np.ones([nx,ny,nz]) #ones mask
    #apply mask for each face
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
        #slice plain out of mask
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(a*x[i] + b*y[j] + c*z[k] < d):
                        grid[i,j,k] = 0

    points = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if(grid[i,j,k] == 1):
                    points.append([x[i],y[j],z[k]])


    dt_image = ndimage.distance_transform_edt(grid)
    r = max(dt_image.flatten())
    cent_inscribed = np.where(dt_image == r)
    r = r*div

    cent_inscribed = list(zip(cent_inscribed[0],cent_inscribed[1],cent_inscribed[2]))[0]
    inscribed_circle = [x[cent_inscribed[0]], y[cent_inscribed[1]], z[cent_inscribed[2]], r]
    return inscribed_circle


#=============================================
def output(points,vertices,particles):
    # Write dump files

    ## Uncomment to print voronoi seeds and vertices 
    # with open('test_points.txt','w') as file:
    #     file.writelines([f'{len(points)}\n','\n'])
    #     for i in points:
    #         file.write(f'{i[0]} {i[1]} {i[2]}\n')

    # with open('test_vertices.txt','w') as file:
    #     file.writelines([f'{len(vertices)}\n','\n'])
    #     for i in vertices:
    #         file.write(f'{i[0]} {i[1]} {i[2]}\n')

    with open('clump.txt','w') as file:
        # file.writelines([f'{len(particles)}\n','\n']) print number of particles first for LAMMPS Style output
        for i in particles:
            file.write(f'{i[0]} {i[1]} {i[2]} {i[3]}\n')

    #Show size dist
    data = pd.DataFrame(particles)
    data.columns = ['x', 'y', 'z','r']
    data['r'].plot.hist(bins=40)
    plt.show()

if __name__ == '__main__':
    main()