import pyvoro as voro
import trimesh
import numpy as np
import pandas as pd
import cvxpy as cp
from matplotlib import pyplot as plt
import sys
import time

if(len(sys.argv) != 4):
    print("Error: Input arguments are of the form: \nvor_pack.py input_geom input_PSD output_file\n PLease check inputs and try again")
    exit()

if(sys.argv[1][-4:].lower() != '.stl'):
    print('Error: Input geometry must be an stl file')
    exit()

#load mesh
mesh = trimesh.load(sys.argv[1])
input_PSD = sys.argv[2]
output_file = sys.argv[3]

bb = mesh.bounds
width, depth, height = bb[1]-bb[0]
#scale bounding box to longest side
#produces better results in thin particles, no effect for more cubic bb's
width = max(width, depth, height)
depth = max(width, depth, height)
height = max(width, depth, height)
xmin,ymin,zmin = bb[0]
xmax,ymax,zmax = xmin+width,ymin+depth,zmin+height

#user defined size distribution
try:
    size_dist = pd.read_csv(f'./{input_PSD}',delimiter='\t', names=['diam','vol_passing'])
except:
    'failed to read input PSD, format: two tab seperated columns: diam  percent_passing'
size_dist['diam'] = size_dist['diam']/1000 #convert from μm to m
size_dist['rad'] = size_dist['diam']/2
#create frequency column
size_dist['freq'] = size_dist['vol_passing']/(size_dist['vol_passing'] * (size_dist['rad'])**3)
size_dist['freq'] = size_dist['freq']/(size_dist['freq'].sum())
size_dist['freq'] = size_dist['freq'].fillna(0.)
#create cumulative frequency column
cfd = [0]
for i,v in enumerate(size_dist['freq'][1::]):
    cfd.append(v+cfd[i])
size_dist['cfd'] = cfd

#print PSD to output as a sanity check
print('User defined PSD\n')
print(size_dist,'\n')

#check particle size ratio
dmin = size_dist['diam'].min()
dmax = size_dist['diam'].max()
lambda_ratio = dmax/dmin
if lambda_ratio > 10:
    print('\n**Warning** particle size ratio Lambda is greater than 10\nClump generation is not guaranteed accurate at this stage, consider reducing the range of your PSD\n')


#=============================================
def main():
    start = time.perf_counter()
    ##PARTICLE SCATTERING
    print('scattering particles...')
    tic = time.perf_counter()
    particles = scatter_particles()
    toc = time.perf_counter()
    print('COMPLETE')
    print(f'generated {len(particles)}')
    print(f'time elapsed {toc - tic}\n')

    ##LLOYDS
    print('Voronoi relaxation...')
    tic = time.perf_counter()
    vor = lloyds(particles,100)
    toc = time.perf_counter()
    print('COMPLETE')
    print(f'time elapsed {toc - tic}\n')

    ##INSCRIBE PARTICLES
    print('Placing particles...')
    tic = time.perf_counter()
    particles = place_particles(vor)
    toc = time.perf_counter()
    print('COMPLETE')
    print(f'time elapsed {toc - tic}\n')
    
    ##TRIM TO MESH
    print('Trimming to mesh...')
    tic = time.perf_counter()
    particles = trim_mesh(particles)
    toc = time.perf_counter()
    print('COMPLETE')
    print(f'time elapsed {toc - tic}\n')

    ##OUTPUTS
    print('Packing COMPLETE, Outputing results..')
    stop = time.perf_counter()
    print(f'total time take {stop - start}\n')
    output(particles)


#=============================================
def generate_radius() -> float:
    '''Generates a new particle radius sampled from the user defined PSD using a linear interpolation between sive sizes'''
    rnd =  np.random.random()
    i = 0
    while rnd > size_dist['cfd'][i] and i<len(size_dist['cfd']): i+=1
    return (size_dist['rad'][i-1] + ((rnd - size_dist['cfd'][i-1])*(size_dist['rad'][i]-size_dist['rad'][i-1]))/(size_dist['cfd'][i]-size_dist['cfd'][i-1]))


#=============================================
def scatter_particles() -> np.ndarray:
    '''Scatter particles randomly accross domain up to a prescribed packing density'''
    volume = width*height*depth
    TARGET_DEN = 0.6
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
def lloyds(particles: np.ndarray,steps: np.ndarray) -> np.ndarray:
    '''Implementation of Lloyds algorythm, with growing cell weights'''
    scale = 0.0 #cell weight scaling [0:1]
    growth_rate = 0.1 #controls rate of increase of scale
    energy_func = 1 #Energy function stopping criteria
    step = 0
    while(scale<1 or energy_func > 0.3):
        try:
            #Build voronoi based on current particle positions and scale
            vor = voro.compute_voronoi(
                particles[:,:3],
                [[xmin,xmax], [ymin,ymax], [zmin,zmax]],
                2.0,
                periodic=[True]*3,
                radii=particles[:,3]*scale
            )
        except:
            #The voronoi can fail when a small particle get's "swallowed" by a larger
            #reducing growth rate will help stop this but make the process slower
            print('growth rate too high, reducing and continuing..')
            scale -= growth_rate
            growth_rate /= 2
            if growth_rate < 0.0005:
                print('Clump generation falied lambda too high, exiting')
                exit()
        else:
            energy_func = 0
            #compute centroid for each cell
            for i,cell in enumerate(vor):
                centroid = [0,0,0]
                for j in cell['vertices']:
                    centroid[0] += j[0]
                    centroid[1] += j[1]
                    centroid[2] += j[2]
                centroid = [j/len(cell['vertices']) for j in centroid]
                #energy function proportional to particle movement
                energy_func += ((particles[i][0]-centroid[0])**2
                                +(particles[i][1]-centroid[1])**2
                                +(particles[i][2]-centroid[2])**2)**0.5
                #move particles
                particles[i][0] = centroid[0]
                particles[i][1] = centroid[1]
                particles[i][2] = centroid[2]
            #grow scale factor
            if(scale < 1):
                scale += growth_rate
            else:
                scale = 1
            energy_func /= len(particles)
            print(f'Voronoi relaxation step {step}',end='\r')
            step += 1
    print(f'Tolerance reached after step {step}')
    return vor


#=============================================
def inscribe_particle(faces: np.ndarray, vertices: np.ndarray) -> list:
    '''Draws inscribed sphere inside convex polyhedra, useing cvxpy module and convex optimisation'''
    centroid = np.array([vertices[:,0].sum(),vertices[:,1].sum(),vertices[:,2].sum()])/len(vertices)
    B = np.ones(len(vertices))
    A = []
    count = 0
    #find constraint equations describing cell faces
    for face in faces:
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
    #find max radius and Chebeshev center
    radius = cp.Variable(1)
    center = cp.Variable(3)
    constraints = [A[i].T@center + np.linalg.norm(A[i], 2)*radius <= B[i] for i in range(len(A))]
    objective = cp.Maximize(radius)
    p = cp.Problem(objective, constraints)
    result = p.solve()
    return([center.value[0],center.value[1],center.value[2],radius.value[0]])


#=============================================
def place_particles(vor: np.ndarray) -> np.ndarray:
    '''places a particle inside each voronoi cell'''
    particles = []
    for i,region in enumerate(vor):
        print(f'placing particle {i}/{len(vor)}',end='\r')
        particles.append(inscribe_particle([face['vertices'] for face in region['faces']],np.array(region['vertices'])))
    particles = np.array(particles)
    print(f'placing particle {len(vor)}/{len(vor)}')
    return particles


#=============================================
def trim_mesh(particles: np.ndarray) -> np.ndarray:
    '''trims particles that have centers outside the mesh'''
    exists = []
    chunk = 500
    i=1
    while chunk*i < len(particles):
        exists += mesh.contains(particles[chunk*(i-1):chunk*i,:3]).tolist()
        i+=1
    exists += mesh.contains(particles[chunk*(i-1):len(particles),:3]).tolist()
    particles = np.array([v for i,v in enumerate(particles) if exists[i] == True])
    return particles


#=============================================
def output(particles: np.ndarray):
    '''writes particle positions to output_file, optionally uncomment some of the code to produce grpahs'''
    # Write output files
    with open(output_file,'w') as file:
        for i,v in enumerate(particles):
            file.write(f'{v[0]} {v[1]} {v[2]} {v[3]}\n')

    ##Uncomment section below to show PSD from the clump generation

    # cvd = [0]
    # for i,v in enumerate(size_dist['vol_passing'][1:]):
    #     cvd.append(v+cvd[i])
    # vol_passing = []
    # for sive_size in size_dist['rad']:
    #     vol_passing.append(sum([(4/3)*np.pi*radius**3 for radius in particles[:,3] if radius < sive_size]))
    # vol_passing = np.array(vol_passing)*100/vol_passing[-1]
    # plt.plot(size_dist['diam'],cvd)
    # plt.plot(size_dist['diam'],vol_passing)
    # plt.title('real vs sampled PSD siving')
    # plt.xscale('log')
    # plt.xlabel('size (μm)')
    # plt.ylabel('% volume passing')
    # plt.legend(['Real PSD', 'Sampled PSD'])
    # plt.show()

if __name__ == '__main__':
    main()