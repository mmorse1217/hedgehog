import pandas as pd
import numpy as np
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D

def pid_to_color(number):
    colors = ['b', 'r', 'g', 'm', 'c', 'y', 'k']
    return colors[number % len(colors)]

sample_positions = pd.read_csv('subsampled_patches.txt', delim_whitespace=True)
position_index = ['m_0', 'm_1','m_2']
#sample_positions= sample_positions[position_index]
f = pl.figure()
ax = f.add_subplot(111, projection='3d')
ax.scatter(
        sample_positions[position_index[0]],
        sample_positions[position_index[1]],
        sample_positions[position_index[2]],
        marker='o',
        c=map(pid_to_color, sample_positions['v'])
        )
pl.show()
