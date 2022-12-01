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
import cvxpy as cp

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

width, depth, height = bb[1]-bb[0]
width = max(width, depth, height)
depth = max(width, depth, height)
height = max(width, depth, height)
xmin,ymin,zmin = bb[0]
# xmax,ymax,zmax = bb[1]
xmax,ymax,zmax = xmin+width,ymin+depth,zmin+height
# width, depth, height = bb[1]-bb[0]

print(xmin,xmax)
print(ymin,ymax)
print(zmin,zmax)

#user defined size distribution Size distribution
#read from csv file format: two tab seperated columns: Diamater  volume/mass percentage
size_dist = pd.read_csv(f'./{input_PSD}',delimiter='\t', names=['diam','vol']) #users sould uses vol percent but scaling will be the same
size_dist['diam'] = size_dist['diam']/1000 #convert from m to μm
size_dist['rad'] = size_dist['diam']/2
size_dist['freq'] = size_dist['vol']/(size_dist['vol'] * (size_dist['rad'])**3)
size_dist['freq'] = size_dist['freq']/(size_dist['freq'].sum())
size_dist['freq'] = size_dist['freq'].fillna(0.)
# cfd = [size_dist['freq'][0]]
cfd = [0]
for i,v in enumerate(size_dist['freq'][:-1]):
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
    particles = disk_sampling()
    # output([1,1,1,1],points,output_file)
    toc = time.perf_counter()
    print('COMPLETE\n')
    print(f'time elapsed {toc - tic}\n')

    ##LLOYDS
    print('Voronoi relaxation...')
    tic = time.perf_counter()
    vor,points = lloyds(particles,100)
    toc = time.perf_counter()
    print('COMPLETE')
    print(f'time elapsed {toc - tic}\n')

    ##PARTICLE PLACEMENT
    print('Placing particles...')
    tic = time.perf_counter()
    particles = []
    for i,region in enumerate(vor):
        print(f'placing particle {i}/{len(vor)}',end='\r')
        particles.append(circ_inscribe([face['vertices'] for face in region['faces']],np.array(region['vertices'])))
    particles = np.array(particles)
    toc = time.perf_counter()
    print('COMPLETE')
    print(f'time elapsed {toc - tic}\n')
    
    ##TRIM TO MESH
    # print('Trimming to mesh...')
    # tic = time.perf_counter()
    # exists = []
    # exists_ = []
    # chunk = 500
    # i=1
    # while chunk*i <len(particles):
    #     exists += mesh.contains(particles[chunk*(i-1):chunk*i,:3]).tolist()
    #     exists_ += mesh.contains(points[chunk*(i-1):chunk*i,:3]).tolist()
    #     i+=1
    # exists += mesh.contains(particles[chunk*(i-1):len(particles),:3]).tolist()
    # exists_ += mesh.contains(points[chunk*(i-1):len(points),:3]).tolist()
    # particles = [v for i,v in enumerate(particles) if exists[i] == True]
    # points = [v for i,v in enumerate(points) if exists[i] == True]
    # toc = time.perf_counter()
    # print('COMPLETE')
    # print(f'time elapsed {toc - tic}\n')

    ##OUTPUTS
    print('Packing COMPLETE, Outputing results..\n')
    stop = time.perf_counter()
    print(f'total time take {stop - start}\n')
    output(particles,points)


#=============================================
def generate_radius():
    #This function generates a new particle radius sampled from the user defined PSD
    rnd =  np.random.random()
    i = 0
    while rnd > size_dist['cfd'][i+1] and i<len(size_dist['cfd']): i+=1
    return (size_dist['rad'][i] + ((rnd - size_dist['cfd'][i])*(size_dist['rad'][i+1]-size_dist['rad'][i]))/(size_dist['cfd'][i+1]-size_dist['cfd'][i]))


#=============================================
def disk_sampling():
    volume = width*height*depth
    TARGET_DEN = 0.64
    pack_den = 0.0
    radii = []
    while(pack_den < TARGET_DEN):
        radius = generate_radius()
        pack_den += (4/3*np.pi * radius**3)/volume
        radii.append(radius)
    
    pos = []
    for i in range(len(radii)):
        pos.append(np.random.random(3)*[width,height,depth]+[xmin,ymin,zmin])
    pos = np.array(pos)

    particles = np.array((pos[:,0],pos[:,1],pos[:,2],radii)).T
    return particles


#=============================================
def lloyds(particles,steps):
    scale = 0
    with open(output_file,'a') as file:
        for step in range(steps):
            vor = voro.compute_voronoi(
                particles[:,:3],
                [[xmin,xmax], [ymin,ymax], [zmin,zmax]], # limits
                2.0, # block size
                periodic=[True]*3,
                radii=particles[:,3]*scale
            )
            for i,cell in enumerate(vor):
                centroid = [0,0,0]
                for j in cell['vertices']:
                    centroid[0] += j[0]
                    centroid[1] += j[1]
                    centroid[2] += j[2]
                centroid = [j/len(cell['vertices']) for j in centroid]
                particles[i][0] = centroid[0] +0.01*(random.random()-0.5)
                particles[i][1] = centroid[1] +0.01*(random.random()-0.5)
                particles[i][2] = centroid[2] +0.01*(random.random()-0.5)

            if(scale < 1):
                scale += 0.01
            else:
                scale = 1

            # Write output files
            file.writelines([f'{len(particles)}\n','\n']) #print number of particles first for LAMMPS Style output
            for i,v in enumerate(particles):
                file.write(f'{v[0]} {v[1]} {v[2]} {v[3]*scale}\n')
            print(f'iteration {step}/{steps}',end='\r')

    return vor,particles

#=============================================
def circ_inscribe(faces,vertices):
    centroid = np.array([vertices[:,0].sum(),vertices[:,1].sum(),vertices[:,2].sum()])/len(vertices)
    B = np.ones(len(vertices))
    A = []
    count = 0
    for face in faces:
        #find equation of plane
        p0, p1, p2 = np.array(vertices[face[0]]), np.array(vertices[face[1]]), np.array(vertices[face[2]])
        v1,v2 = p2-p0, p1-p0
        cross_product = np.cross(v1,v2)
        a,b,c = cross_product
        d = np.dot(cross_product,p2)
        a,b,c = a/d,b/d,c/d
        #flip plane if contains centroid
        if(a*centroid[0] + b*centroid[1] + c*centroid[2] > 1):
            a, b, c, d = -a,-b, -c, -d 
            B[count] = - B[count]
        count += 1
        A.append([a,b,c])

    A = np.array(A)
    # # variables
    radius = cp.Variable(1)
    center = cp.Variable(3)

    constraints = [A[i].T@center + np.linalg.norm(A[i], 2)*radius <= B[i] for i in range(len(A))]
    objective = cp.Maximize(radius)
    p = cp.Problem(objective, constraints)
    # The optimal objective is returned by p.solve().
    result = p.solve()
    # The optimal value
    return([center.value[0],center.value[1],center.value[2],radius.value[0]])



#=============================================
def output(particles,points):
    # # uncomment for points output    
    with open('points.txt','w') as file:
        file.writelines([f'{len(particles)}\n','\n']) #print number of particles first for LAMMPS Style output
        for i in points:
            file.write(f'{i[0]} {i[1]} {i[2]} {i[3]}\n')

    # Write output files
    with open(output_file,'a') as file:
        file.writelines([f'{len(particles)}\n','\n']) #print number of particles first for LAMMPS Style output
        for i,v in enumerate(particles):
            file.write(f'{v[0]} {v[1]} {v[2]} {v[3]}\n')

    #Show size dist
    data = pd.DataFrame(points)
    data.columns = ['x', 'y', 'z','r']
    data['r'].plot.hist(bins=40)

    #Show size dist
    data = pd.DataFrame(particles)
    data.columns = ['x', 'y', 'z','r']
    data['r'].plot.hist(bins=40)
    
    # plt.yscale('log')
    plt.plot(size_dist['rad'],size_dist['freq']*len(particles)*0.5)

    plt.title('Particle radii resulting from different packing methods')
    plt.ylabel('frequency')
    plt.xlabel('radius (mm)')
    plt.legend(['Exact packing','Dense packing','Particle size distribution'])
    plt.show()


if __name__ == '__main__':
    main()