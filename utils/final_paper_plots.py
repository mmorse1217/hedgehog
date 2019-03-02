from read import load_pydict,load_vtk_point_data, generate_latex_table
import numpy as np
import pylab as pl
import vtk
import glob
import os


def get_prop_list(results, s):
    return np.array([r[s] for r in results])

def load_test_results(path_to_results, file_prefix):
    results = []
    for result_file in glob.glob(path_to_results+'*.pydict'):
        result = load_pydict(result_file)
        if file_prefix in result_file.split('/')[-1]:
            results.append(result)
    return results

def fit_line(x,y):
    x = np.log10(x)
    y = np.log10(y)
    n = x.shape[0]
    assert x.shape == y.shape
    A = np.hstack([np.ones((n,1)),x.reshape(n,1)])
    return np.linalg.lstsq(A,y, rcond=1e-13)

def compute_rate(error,spacing):
    assert error.shape == spacing.shape
    n = error.shape[0]
    rate = np.zeros((n-1,))
    for i in xrange(n-1):
        #print spacing[i+1],spacing[i], error[i+1], error[i]
        rate[i] = (spacing[i+1] -spacing[i])/(error[i+1] - error[i])
    return rate

def extract_props_from_results(properties, results):
    return { prop: get_prop_list(results, prop) for prop in properties}

def plot_greens_identity(path_to_tests, out_dir):
    domain_names = {'cube': 'Cube', 'newtorus':'Torus'}
    for domain in domain_names.keys():
        pl.figure()
        for kernel in ['laplace', 'stokes', 'navier']:
            test_name = kernel + '_greens_identity_' + domain
            test_output_path = path_to_tests + test_name + '/output/test_greens_identity/'
            test_output_path = test_output_path + '/' + domain + '/' + kernel +'/'
            
            results = load_test_results(test_output_path, 'qbx')
            
            props_dep = ['max relative error']
            props_indep = ['face-map max patch size', 'coarse spacing']
            prop_values_dep = extract_props_from_results(props_dep, results)
            prop_values_indep = extract_props_from_results(props_indep, results)

            x = prop_values_indep['face-map max patch size']#/\
                    #np.floor(1./prop_values_indep['coarse spacing'])+1.
            y = prop_values_dep['max relative error']
            print 'convergence rate', kernel, domain,':',fit_line(x,y)[0]
            #print y
            print 'convergence rate',compute_rate(np.log10(x), np.log10(y))
            pl.loglog(x,y, 'o-',label=kernel.capitalize())
        pl.legend(loc='best')#, bbox_to_anchor=(1, .5))
        pl.title('Greens Identity Error Convergence: '+domain_names[domain])
        pl.ylabel(r'$\log\left(\|u- \tilde{u}\|_\infty/\|u\|_\infty\right)$')
        pl.xlabel(r'$L_{max}$')
        pl.savefig(out_dir+'greens_identity_'+domain,
                bbox_inches='tight',dpi=150)

def plot_solver_convergence(path_to_tests, out_dir):
    domain_names = {'cube': 'Cube', 'pipe':'Pipe', 'ttorus2': 'Genus 2 surface'}
    for domain in domain_names.keys():
        pl.figure()
        for kernel in ['laplace', 'stokes', 'navier']:
            test_name = kernel + '_qbx_eval_conv_' + domain
            test_name = test_name[:-1] if domain is 'ttorus2' else test_name
            test_output_path = path_to_tests + test_name + '/output/test_qbkix_eval_convergence/'
            test_output_path = test_output_path + '/' + domain + '/' + kernel +'/'
            
            results = load_test_results(test_output_path, 'qbx')
            
            props_dep = ['max relative error']
            props_indep = ['face-map max patch size', 'coarse spacing']
            prop_values_dep = extract_props_from_results(props_dep, results)
            prop_values_indep = extract_props_from_results(props_indep, results)

            x = prop_values_indep['face-map max patch size']#/\
                    #np.floor(1./prop_values_indep['coarse spacing'])+1.
            y = prop_values_dep['max relative error']
            print 'convergence rate', kernel, domain,':',fit_line(x,y)[0]
            pl.loglog(x,y,'o-', label=kernel.capitalize())
        
        pl.legend(loc='best')#, bbox_to_anchor=(1, .5))
        pl.title('Solver Convergence: '+domain_names[domain])
        pl.ylabel(r'$\log\left\|u- \tilde{u}\|_\infty\|u\|_\infty\right)$')
        pl.xlabel(r'$L_{max}$')
        pl.savefig(out_dir+'solver_convergence_'+domain,
                bbox_inches='tight',dpi=150)

def plot_comparison_error(path_to_tests, out_dir):
    test_paths = {
            'laplace_qbx_vs_singular_quad': {
                'domain' : 'cube',
                'kernel' : 'laplace',
                'dir'    : 'cube/laplace/'
                }, \
            'navier_qbx_vs_singular_quad': {
                'domain' : 'cube',
                'kernel' : 'navier',
                'dir'    : 'cube/navier/'
                }, \
            'sph_harm_qbx_vs_singular_quad': {
                'domain' : 'sphere',
                'kernel' : 'laplace',
                'dir'    : 'sphere/laplace/'
                }, \
            'navier_qbx_vs_singular_quad_const_den': {
                'domain' : 'cube_const_den',
                'kernel' : 'navier',
                'dir'    : 'const_density/cube/navier/'
                }, \
            'laplace_qbx_vs_singular_quad_solve': {
                'domain' : 'cube_solve',
                'kernel' : 'laplace',
                'dir'    : 'solve/cube/laplace/'
                }, \
            'face_map_qbx_vs_singular_quad': {
                'domain' : 'cube',
                'kernel' : 'laplace',
                'dir'    : 'face_map_vs_blended/cube/laplace/'
                },\
            'face_map_qbx_vs_singular_quad_navier': {
                'domain' : 'cube',
                'kernel' : 'navier',
                'dir'    : 'face_map_vs_blended/cube/navier/'
                }
            }
    def plot_error_vs_coarse_spacing(results, *args, **kwargs):
        props_dep = ['max absolute error']
        props_indep = ['coarse spacing']
        prop_values_dep = extract_props_from_results(props_dep, results)
        prop_values_indep = extract_props_from_results(props_indep, results)
        x = prop_values_indep['coarse spacing']
        y = prop_values_dep['max absolute error']
        pl.loglog(x,y, *args, **kwargs)

    def plot_error_vs_face_map_size(results, *args, **kwargs):
        props_dep = ['max absolute error']
        props_indep = ['coarse spacing','face-map max patch size']
        prop_values_dep = extract_props_from_results(props_dep, results)
        prop_values_indep = extract_props_from_results(props_indep, results)
        x = prop_values_indep['face-map max patch size']/\
               (np.floor(1./prop_values_indep['coarse spacing'])+1.)
        y = prop_values_dep['max absolute error']
        pl.loglog(x,y, *args, **kwargs)

    # plot error comparision
    '''
    for test_name, test_data in test_paths.items()[:-2]:


        pl.figure()
        test_output_path = path_to_tests + test_name + '/output/test_qbkix_vs_singular_quad/' + test_data['dir']
        #test_output_path = test_output_path + '/' + domain + '/' + kernel +'/'
        qbx_results = load_test_results(test_output_path, 'qbx')
        singular_results = load_test_results(test_output_path, 'singular_eval')
        plot_error_vs_coarse_spacing(qbx_results, 'o-', label='QBX')
        plot_error_vs_coarse_spacing(singular_results, 'o-', label='Singular Quad.')
        
        #pl.legend(loc='center left', bbox_to_anchor=(1, .5))
        pl.legend(loc='best')#, bbox_to_anchor=(1, .5))
        pl.title('QBX vs. Singular Quadrature: ' + test_data['domain'].capitalize()\
                +', ' + test_data['kernel'].capitalize())
        pl.ylabel(r'$\log\|u- \tilde{u}\|_\infty$')
        pl.xlabel(r'$h$')
        pl.savefig(out_dir+'comparison_'+test_data['domain']+'_'+test_data['kernel']+'_error',  bbox_inches='tight',dpi=150)

    '''

    #pl.figure()
    f, ax = pl.subplots(ncols=2, nrows=1, sharex=True, sharey=True,
            figsize=(6,3))
    test_name = 'face_map_qbx_vs_singular_quad'
    test_data = test_paths['face_map_qbx_vs_singular_quad']
    test_output_path = path_to_tests + test_name + '/output/test_qbkix_vs_singular_quad/' + test_data['dir']
    #test_output_path = test_output_path + '/' + domain + '/' + kernel +'/'
    qbx_results = load_test_results(test_output_path, 'qbx')
    singular_results = load_test_results(test_output_path, 'singular_eval')
    pl.subplot(121)
    plot_error_vs_face_map_size(qbx_results, 'bo-', label='QBX')
    pl.legend()
    pl.ylabel(r'$\log\|u- \tilde{u}\|_\infty$')
    pl.xlabel(r'$L_{max}$')
    pl.subplot(122)
    plot_error_vs_coarse_spacing(singular_results, 'ro-', label='Singular Quad.')
    pl.legend()
    pl.xlabel(r'$h$')
    
    #pl.legend(loc='best')#, bbox_to_anchor=(1, .5))
    #pl.title('QBX + Face-map vs. Singular Quadrature + Blended: ' + test_data['domain'].capitalize()\
            #+', ' + test_data['kernel'].capitalize())
    #pl.ylabel(r'$\log\|u- \tilde{u}\|_\infty$')
    #f.suptitle('QBX vs. Singular Quadrature: ' + test_data['domain'].capitalize()\
            #+', ' + test_data['kernel'].capitalize())
    f.suptitle('QBX + Face-map vs. Singular Quadrature + Blended:')
    pl.savefig(out_dir+'comparison_face_map_vs_blendsurf_'+test_data['domain']+'_'+test_data['kernel']+'_error',  bbox_inches='tight',dpi=150)
    pl.figure()
    
    test_name = 'face_map_qbx_vs_singular_quad_navier'
    test_data = test_paths['face_map_qbx_vs_singular_quad_navier']
    test_output_path = path_to_tests + test_name + '/output/test_qbkix_vs_singular_quad/' + test_data['dir']
    #test_output_path = test_output_path + '/' + domain + '/' + kernel +'/'
    qbx_results = load_test_results(test_output_path, 'qbx')
    singular_results = load_test_results(test_output_path, 'singular_eval')
    plot_error_vs_face_map_size(qbx_results, 'o-', label='QBX')
    plot_error_vs_coarse_spacing(singular_results, 'o-', label='Singular Quad.')
    
    pl.legend(loc='best')#, bbox_to_anchor=(1, .5))
    pl.title('QBX + Face-map vs. Singular Quadrature + Blended: ' + test_data['domain'].capitalize()\
            +', ' + test_data['kernel'].capitalize())
    pl.ylabel(r'$\log\|u- \tilde{u}\|_\infty$')
    pl.xlabel(r'$L_{max}$')
    pl.savefig(out_dir+'comparison_face_map_vs_blendsurf_'+test_data['domain']+'_'+test_data['kernel']+'_error',  bbox_inches='tight',dpi=150)

def plot_comparison_timing(path_to_tests, out_dir):
    const_den = {
            'laplace_qbx_vs_singular_quad': {
                'domain' : 'cube',
                'kernel' : 'laplace',
                'dir'    : 'cube/laplace/'
                }, \
            'navier_qbx_vs_singular_quad_const_den': {
                'domain' : 'cube_const_den',
                'kernel' : 'navier',
                'dir'    : 'const_density/cube/navier/'
                }
            }
    test_paths = {
            'navier_qbx_vs_singular_quad': {
                'domain' : 'cube',
                'kernel' : 'navier',
                'dir'    : 'cube/navier/'
                }, \
            'sph_harm_qbx_vs_singular_quad': {
                'domain' : 'sphere',
                'kernel' : 'laplace',
                'dir'    : 'sphere/laplace/'
                }, \
            'laplace_qbx_vs_singular_quad_solve': {
                'domain' : 'cube_solve',
                'kernel' : 'laplace',
                'dir'    : 'solve/cube/laplace/'
                }, \
            'face_map_qbx_vs_singular_quad': {
                'domain' : 'cube',
                'kernel' : 'laplace',
                'dir'    : 'face_map_vs_blended/cube/laplace/'
                },\
            'face_map_qbx_vs_singular_quad_navier': {
                'domain' : 'cube',
                'kernel' : 'navier',
                'dir'    : 'face_map_vs_blended/cube/navier/'
                }
            }
    def plot_timing_vs_error_qbx(results, ax, *args, **kwargs):
        props_indep = ['total fmm time',  'total qbx time', 'total matvec time', 
                'total density interp time']
        
        props_indep_names = ['FMM',  'extrap', 'total', 
                'density interp']
        props_indep_names = [ 'QBX ' + p for p in props_indep_names]
        
        props_markers = ['*',  '^', 'o', 
                '+']
        props_markers = ['b-'+ m for m in props_markers]

        props_dep = ['max absolute error', 'FMM benchmark time', 'number of GMRES iterations']
        prop_values_dep = extract_props_from_results(props_dep, results)
        prop_values_indep = extract_props_from_results(props_indep, results)
        
        x = prop_values_dep['max absolute error']
        fmm = prop_values_dep['FMM benchmark time']
        gmres_iters = prop_values_dep['number of GMRES iterations']
        plots = []
        for prop, name, marker in zip(props_indep, props_indep_names, props_markers):
            y = prop_values_indep[prop]/fmm/gmres_iters
            if prop is 'total matvec time':
                a = ax[0]
                #pl.subplot(131)
                #pl.loglog(x,y, marker, label=name,*args, **kwargs)
                p, = a.loglog(x,y, marker, label=name,*args, **kwargs)
                #coeffs,_,__,___ = fit_line(np.log10(x),np.log10(y))
                #a.loglog(x, (10**coeffs[0] *x**coeffs[1]), 'k--')
                #plots.append(p)
            #pl.subplot(132)
            #pl.loglog(x,y, marker, label=name,*args, **kwargs)
            a = ax[1]
            p, = a.loglog(x,y, marker, label=name,*args, **kwargs)
            plots.append(p)
        return plots, props_indep_names
    def plot_timing_vs_error_singular(results, ax, *args, **kwargs):
        props_indep = ['total fmm time', 'total inaccurate subtract time', 'total fft upsample time', 'total polar quad time', 'total matvec time']
        
        props_indep_names = ['FMM',  'subtract', 'FFT', 'polar quad', 'total']
        props_indep_names = [ 'Singular ' + p for p in props_indep_names]

        props_markers = ['*',  'D', '+', '^', 'o']
        props_markers = ['r-'+ m for m in props_markers]

        props_dep = ['max absolute error', 'FMM benchmark time', 'number of GMRES iterations']

        prop_values_dep = extract_props_from_results(props_dep, results)
        prop_values_indep = extract_props_from_results(props_indep, results)
        x = prop_values_dep['max absolute error']
        fmm = prop_values_dep['FMM benchmark time']
        gmres_iters = prop_values_dep['number of GMRES iterations']
        plots = []
        for prop, name, marker in zip(props_indep, props_indep_names,
                props_markers):
            y = prop_values_indep[prop]/fmm/gmres_iters
            if prop is 'total matvec time':
                #pl.subplot(131)
                a = ax[0]
                #pl.loglog(x,y, marker, label=name,*args, **kwargs)
                p, = a.loglog(x,y, marker, label=name,*args, **kwargs)
                #plots.append(p)
            #pl.subplot(133)
            a = ax[2]
            #pl.loglog(x,y, marker, label=name,*args, **kwargs)
            p, = a.loglog(x,y, marker, label=name,*args, **kwargs)
            plots.append(p)
        return plots, props_indep_names
    figure_size = (12,4) 
    # plot timing comparision
    for test_name, test_data in test_paths.items()[:-2]:

        f, ax = pl.subplots(ncols=3, nrows=1, sharex=True, sharey=True, figsize=figure_size)
        test_output_path = path_to_tests + test_name + '/output/test_qbkix_vs_singular_quad/' + test_data['dir']
        #test_output_path = test_output_path + '/' + domain + '/' + kernel +'/'
        qbx_results = load_test_results(test_output_path, 'qbx')
        singular_results = load_test_results(test_output_path, 'singular_eval')
        print test_name
        qbx_plots, qbx_names = plot_timing_vs_error_qbx(qbx_results,ax)
        singular_plots,singular_names = plot_timing_vs_error_singular(singular_results,ax)
        
        pl.legend(qbx_plots+ singular_plots,qbx_names+singular_names,
                loc='center left', bbox_to_anchor=(1, .5))
        f.suptitle('QBX vs. Singular Quadrature: ' + test_data['domain'].capitalize()\
                +', ' + test_data['kernel'].capitalize())
        ax[1].set_xlabel(r'$\log\|u- \tilde{u}\|_\infty$')
        ax[0].set_ylabel(r'$T/T_{FMM}$')
        f.savefig(out_dir+'comparison_'+test_data['domain']+'_'+test_data['kernel']+'_timing',  bbox_inches='tight',dpi=150)


    f, ax = pl.subplots(ncols=3, nrows=1, sharex=True, sharey=True, figsize=figure_size)
    test_name = 'face_map_qbx_vs_singular_quad'
    test_data = test_paths['face_map_qbx_vs_singular_quad']
    test_output_path = path_to_tests + test_name + '/output/test_qbkix_vs_singular_quad/' + test_data['dir']
    #test_output_path = test_output_path + '/' + domain + '/' + kernel +'/'
    qbx_results = load_test_results(test_output_path, 'qbx')
    singular_results = load_test_results(test_output_path, 'singular_eval')
    qbx_plots, qbx_names = plot_timing_vs_error_qbx(qbx_results,ax)
    singular_plots, singular_names = plot_timing_vs_error_singular(singular_results,ax)
    
    pl.legend(qbx_plots+ singular_plots, qbx_names + singular_names,
            loc='center left', bbox_to_anchor=(1, .5))
    f.suptitle('QBX + Face-map vs. Singular Quadrature + Blended: ' + test_data['domain'].capitalize()\
            +', ' + test_data['kernel'].capitalize())
    ax[1].set_xlabel(r'$\log\|u- \tilde{u}\|_\infty$')
    ax[0].set_ylabel(r'$T/T_{FMM}$')
    f.savefig(out_dir+'comparison_face_map_vs_blendsurf_'+test_data['domain']+'_'+test_data['kernel']+'_timing',  bbox_inches='tight',dpi=150)

    f, ax = pl.subplots(ncols=3, nrows=1, sharex=True, sharey=True, figsize=figure_size)
    test_name = 'face_map_qbx_vs_singular_quad_navier'
    test_data = test_paths['face_map_qbx_vs_singular_quad_navier']
    test_output_path = path_to_tests + test_name + '/output/test_qbkix_vs_singular_quad/' + test_data['dir']
    #test_output_path = test_output_path + '/' + domain + '/' + kernel +'/'
    qbx_results = load_test_results(test_output_path, 'qbx')
    singular_results = load_test_results(test_output_path, 'singular_eval')
    qbx_plots, qbx_names = plot_timing_vs_error_qbx(qbx_results,ax)
    singular_plots, singular_names = plot_timing_vs_error_singular(singular_results,ax)
    
    pl.legend(qbx_plots+ singular_plots, qbx_names + singular_names,
            loc='center left', bbox_to_anchor=(1, .5))
    f.suptitle('QBX + Face-map vs. Singular Quadrature + Blended: ' + test_data['domain'].capitalize()\
            +', ' + test_data['kernel'].capitalize())
    ax[1].set_xlabel(r'$\log\|u- \tilde{u}\|_\infty$')
    ax[0].set_ylabel(r'$T/T_{FMM}$')
    f.savefig(out_dir+'comparison_face_map_vs_blendsurf_'+test_data['domain']+'_'+test_data['kernel']+'_timing',  bbox_inches='tight',dpi=150)


test_results_dir = 'data/results14/'
out_dir = 'output/plots_final12/'
#plot_greens_identity(test_results_dir, out_dir)
plot_solver_convergence(test_results_dir, out_dir)
#plot_comparison_error(test_results_dir, out_dir)
#plot_comparison_timing(test_results_dir, out_dir)

