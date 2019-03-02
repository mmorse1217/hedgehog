from read import load_pydict,load_vtk_point_data, generate_latex_table
import argparse as ap
import os
import re
import glob
import pylab as pl
import numpy as np
import itertools as it
import pandas as pd

def check_for_results_folder(folder, file_prefix):
    results_for_any_prefix = \
            [len(glob.glob(os.path.join(folder, prefix + '*.pydict'))) > 1 for prefix in file_prefix]
    return any(results_for_any_prefix)

def greens_identity_table(folder):
    pass

def build_table(data_dir, kernels, geometries, test_name, prop):
    table_rows = []
    row_names = []
    # For each test that was run...
    for kernel_name in kernels:
        for geom_name in geometries:
            test_dir = os.path.join(data_dir, kernel_name + test_name + geom_name)
            # there is only one directory in test_dir/output/ that has test results saved.
            # find it with os.walk
            for folder, __, files in os.walk(test_dir):
                if check_for_results_folder(folder, 'qbx_greens_id'):
                    # once we found the results, load the python dicts and 
                    # aggregate the values 'prop' from each result dict into a list
                    result_filenames = [f for f in files if '.pydict' in f]
                    results = [load_pydict(os.path.join(folder,f)) for f in result_filenames]
                    table_row = np.array([r[prop] for r in results])
                    table_rows.append( pd.Series(table_row))
                    row_names.append(kernel_name + ', ' + geom_name)
    table = pd.DataFrame.from_items(zip(row_names, table_rows)).T
    return table

data_dir = 'data/results2/'

# order by geometry, then kernel
kernels = ['laplace', 'navier', 'stokes']
geometries = ['cube', 'newtorus', 'pipe']
prop = 'max absolute error'

print build_table(data_dir, kernels, geometries, '_greens_identity_', prop)
print build_table(data_dir, kernels, geometries, '_qbx_vs_singular_quad_', prop)
#geometries = ['cube',  'pipe', 'ttorus']
#print build_table(data_dir, kernels, geometries, '_qbx_eval_conv_', prop)
