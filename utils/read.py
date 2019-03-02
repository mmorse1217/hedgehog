from ast import literal_eval
from vtk.util.numpy_support import vtk_to_numpy
import vtk
import glob
import os

def load_pydict(filepath):
    #filename = 'test_greens_identity/24_patches_11-15-2017_16-40-53.pydict'
    data = {}
    with open(filepath,'r') as f:
        data.update(literal_eval(f.read()))
    return data

def load_vtk_point_data(filepath):
    reader = vtk.vtkXMLPolyDataReader() 
    reader.SetFileName(filepath) 
    reader.Update()
    polyDataOutput = reader.GetOutput() 
    points = vtk_to_numpy(polyDataOutput.GetPoints().GetData())
    error = vtk_to_numpy(polyDataOutput.GetPointData().GetArray(0))
    return points, error

def collect_max_error(test_path):
    print test_path
    print glob.glob(test_path + '/*.pydict')
    results = []
    for test_case in glob.glob(test_path + '/*.pydict'):
        result = load_pydict(test_case)
        # read in structured data with the same name
        vtk_filename = test_case[:-7]+'.vtp'
        points, error = load_vtk_point_data(vtk_filename)
        max_error = error.max()
        result['max error'] = max_error
        results.append(result)

    return results

def generate_latex_table(test_path): 
    domains = next(os.walk(test_path))[1]
    results = []
    for domain in domains:
        equation_types = next(os.walk(test_path + domain))[1]
        print domain, equation_types
        for equation_type in equation_types:
            test_results = collect_max_error(test_path + domain + '/'+equation_type)
            results += test_results
    print results
    #for i in os.walk(test_path):
    #print i[0]

