import numpy as np
import pandas as p
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D, art3d
from matplotlib import rc, colors
from itertools import izip


def generate_layer_index_set(layers, num_points):
    # Assumptions: - 8 interpolation points generated
    #              - potential_error and point_positions are in order from
    #              inside to outside, in order per target point, per patch
    it = np.arange(num_points/8)
    point_index = np.array([
            np.array(idx) for idx in izip(
                    8*it + i for i in layers
                )
        ])
    point_index = point_index.reshape(num_points/(8/len(layers)))
    return point_index
def plot_all_interior_points(point_positions, error):
    fig = pl.figure()
    #for i in xrange(4):
    #ax = fig.add_subplot(221+i, projection='3d')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')
    #ith_layer_position = po
    ax.scatter(point_positions[dim[0]],
        point_positions[dim[1]],
        point_positions[dim[2]],
            #c=np.log(error),
            c = 'b',
            marker='o')
    pl.show()

def plot_each_layer_interior_points():
    fig = pl.figure()
    #for i in xrange(8):
    for i in xrange(1):
        '''
        if i% 4 is 0:
            fig = pl.figure()
        '''
        ax = fig.add_subplot(111, projection='3d')
        #ax = fig.add_subplot(221+i%4, projection='3d')
        ax.set_aspect('equal')
        point_positions = point_positions_init.take(generate_layer_index_set([i],
            num_points))
        potential_error = potential_error_init.take(generate_layer_index_set([i],
            num_points))
        print generate_layer_index_set([i],  num_points)
        '''
        error = potential_error[dim[0]]**2 + \
                potential_error[dim[1]]**2+\
                potential_error[dim[2]]**2
        '''
        error = potential_error[dim[0]]
        print np.log10(error)
        sc = ax.scatter(point_positions[dim[0]],
            point_positions[dim[1]],
            point_positions[dim[2]],
                c=np.log10(error),
                marker='o'
                )
        '''
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        x = .952*np.outer(np.cos(u), np.sin(v))
        y = .952*np.outer(np.sin(u), np.sin(v))
        z = .952*np.outer(np.ones(np.size(u)), np.cos(v))
        s = ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b')
        s.set_facecolor((0,0,1, .1))
        s.set_edgecolor((0,0,1, .1))
        '''
        pl.colorbar(sc);
    pl.show()

def plot_distance_to_targets(targets_init):
    fig = pl.figure()
    for i in xrange(1):
        ax = fig.add_subplot(111+i, projection='3d')
        #ax = fig.add_subplot(221+i, projection='3d')
        ax.set_aspect('equal')
        point_positions = point_positions_init.take(generate_layer_index_set([i],
            num_points))
        point_positions = point_positions.set_index(np.arange(num_points/8))
        '''
        repeated_target_index = np.array([
               np.array([j for j in xrange(1)]) for i in xrange(num_points/8)
            ]).reshape((num_points/8,))
        '''
        repeated_target_index = np.array([i for i in xrange(num_points/8)])
        targets = targets_init.take(repeated_target_index)
        '''
        distance_to_target = \
                (targets[dim[0]] - point_positions[dim[0]])**2 + \
                (targets[dim[1]] - point_positions[dim[1]])**2+\
                (targets[dim[2]] - point_positions[dim[2]])**2
        '''
        distance_to_target = \
                (point_positions[dim[0]])**2 + \
                (point_positions[dim[1]])**2+\
                (point_positions[dim[2]])**2
        distance_to_target = np.sqrt(distance_to_target)
        p.set_option('display.max_rows', len(distance_to_target))
        print p.concat([distance_to_target, point_positions], axis=1)
        quit()
        sc = ax.scatter(point_positions[dim[0]],
            point_positions[dim[1]],
            point_positions[dim[2]],
                c=distance_to_target,
                marker='o'
                )
        pl.colorbar(sc);
    pl.show()



rc('text', usetex=True)
#num_points = 38400*2
num_points = 12305
potential_error_init = p.read_csv('interpolation_point_errors.txt',
        delim_whitespace=True, nrows=num_points)
point_positions_init = p.read_csv('interpolation_point_positions.txt',
        delim_whitespace=True, nrows=num_points)
point_positions_init = p.read_csv('interpolation_points_true_potential.txt',
        delim_whitespace=True, nrows=num_points)
targets_init = p.read_csv('target_point_positions.txt', delim_whitespace=True,
        nrows = num_points/8)
#potential_error = potential_error[::]
#point_positions = point_positions[::]

# Number of on-surface target points loaded
#interior_point_index =  interior_point_index.reshape(num_points/2) 
#interior_point_index = generate_layer_index_set([0,1,2,3], num_points)
interior_point_index = generate_layer_index_set([0], num_points)
potential_error = potential_error_init.take(interior_point_index)
point_positions = point_positions_init.take(interior_point_index)

dim = ['m_0', 'm_1','m_2']
dim = ['pos_1', 'pos_2','pos_3']
plot_all_interior_points(point_positions, None)
#plot_each_layer_interior_points()
#plot_distance_to_targets(targets_init)
