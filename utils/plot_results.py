from read import load_pydict,load_vtk_point_data, generate_latex_table
import argparse as ap
import os
import re
import glob
import pylab as pl
import numpy as np
import itertools as it
from matplotlib.font_manager import FontProperties


class CustomArgFormatter(ap.ArgumentDefaultsHelpFormatter,
                         ap.RawDescriptionHelpFormatter):
    pass

class Plot(object):
    def __init__(self,**kwargs):
        self._opts = self._parse_cl(**kwargs)

        self.set_defaults()

    def __getattr__(self,name):
        return self._opts[name]

    def _parse_cl(self,callback=None,predoc=None,postdoc=None):
        # commandline parser
        clp = ap.ArgumentParser(
            formatter_class=CustomArgFormatter,
            description=predoc,epilog=postdoc,
            )
        clp.add_argument('--semilogx', help='plot x-axis in semilog', action='store_true')
        clp.add_argument('--semilogy', help='plot y-axis in semilog', action='store_true')
        clp.add_argument('--loglog', help='x and y-axes in loglog', action='store_true')
        clp.add_argument('--xvar', help='specifies which property in the results pydict to extract as the x-axis of the plot',action='append')
        clp.add_argument('--yvar', help='specifies which property in the results pydict to extract as the y-axis of the plot. Multiple values correspond to distinct calls to pl.plot()', action='append')
        clp.add_argument('--prefix', help='categories of files with prefixes', action='append')
        clp.add_argument('--label', help='legend labels for each yvar', action='append')
        clp.add_argument('--legend', help='show a legend on the plot', action='store_true')
        clp.add_argument('--name', help='name of each prefix category to put in a legend', action='append')
        clp.add_argument('--scale-by', help='property to scale all other by (y-vars only)',default='')
        clp.add_argument('--output-file', help='name of final plot to be save')
        clp.add_argument('--x-label', help='text to insert on x-axes (copied verbatim). please escape $ for latex (i.e. pass \$)')
        clp.add_argument('--y-label', help='text to insert on y-axes (copied verbatim). please escape $ for latex (i.e. pass \$)')
        clp.add_argument('--filename', help='filename of final plotted output. by default this is prepended by problem relevant data (geometry, kernel, test type)')
        clp.add_argument('--outdir', help='path to output directory for plots. current directory by default', default='')
        clp.add_argument('--indir', help='path to input data for plots. current directory by default', default='data/')
        
        if callback is not None: 
            callback(clp)
        return vars(clp.parse_args())

    def set_defaults(self):
        pass

    def aggregate_results(self, folder):
        results = {} 
        for file_prefix in self.prefix:
            category_results = []
            
            # collect results dictionaries into a list and associate it with the prefix
            pattern = os.path.join(folder, file_prefix + '*.pydict')
            
            for result_file in glob.glob(pattern):
                result = load_pydict(result_file)
                category_results.append(result)

            results[file_prefix] = category_results
        return results

    # folder - directory of interest
    # returns - boolean indicating whether any test results live in this folder
    def check_for_results_folder(self, folder):
        results_for_any_prefix = \
                [len(glob.glob(os.path.join(folder, prefix + '*.pydict'))) > 1 for prefix in self.prefix]
        return any(results_for_any_prefix)

    
    # results - list of python dicts containing results
    # key - string of interest
    # returns - list of values corr to key in each element of results
    def extract(self, results, key):
        try:
            return np.array([r[key] for r in results])
        except KeyError:
            # if either of the requested values aren't present, skip plotting them
            print 'failed to extract key ', key
            raise KeyError

    def filename_from_folder(self, folder):
        '''
        filename = re.sub('output\/', '', folder)
        filename = re.sub('\/', '_', filename)
        filename = os.path.join(self.outdir, filename)
        filename = filename + '_' + self.output_file
        '''
        filename = folder.split('/')[2]

        filename = os.path.join(self.outdir, filename)
        filename = filename + '_' + self.output_file
        return filename

    def plot(self, folder, filename=None, *args, **kwargs):
        if not filename:
            filename = self.filename_from_folder(folder)
        
        # need an equal number of plotting labels and  y variables to plot on  ..
        assert len(self.yvar) is len(self.label)

        # need an equal number of x- and y- vars, or alternatively, if only one 
        # x-var is passed, it is assumed to plot several y-vars with the same 
        # x-var value. only 1 or n x-vars allowed, then
        assert len(self.xvar) is len(self.yvar) or len(self.xvar) is 1
        
        # accumulate results in the categories that we specified
        test_results= self.aggregate_results(folder)
        assert len(test_results) is len(self.prefix)
      
        colormap = pl.cm.viridis
        k = len(self.name)
        color = it.cycle((colormap(float(i)/float(k)) for i in xrange(k)))
        marker_list = ('o', '*', '^', '+', 'x','D')
        marker = it.cycle(marker_list)
        # make a new plot
        pl.figure()
        # iterate over categories
        for name, results_per_category in it.izip(self.name,test_results.itervalues()):
            plot_color = color.next()
            marker = it.cycle(marker_list)
            x_var_iter = self.xvar if len(self.xvar) is len(self.yvar) \
                    else [self.xvar[0] for i in xrange(self.yvar)]
            # for each specified y-var, plot a curve
            print x_var_iter, self.yvar, self.label
            for x_var, y_var, label in it.izip(x_var_iter, self.yvar, self.label):
                # get x- and y-values to plot over
                y_values = None
                x_values = None
                print 'getting', x_var, y_var
                try:
                    y_values = self.extract(results_per_category, y_var)
                    x_values = self.extract(results_per_category, x_var)
                except KeyError:
                    # if either of the requested values aren't present, skip plotting them
                    continue
                
                plot_marker = marker.next()

                print x_values
                print y_values
                print ('%s %s' % (name, label))
                plot_props = {
                        'label': ('%s %s' % (name, label)), 
                        'color': plot_color, 
                        'marker': plot_marker
                        }
                if self.semilogx:
                    pl.semilogx(x_values, y_values, **plot_props)
                elif self.semilogy:
                    pl.semilogy(x_values, y_values, **plot_props)
                elif self.loglog:
                    pl.loglog(x_values, y_values, **plot_props)
                else:
                    pl.plot(x_values, y_values, **plot_props)

        if(self.legend):
            pl.legend(loc='center left', bbox_to_anchor=(1, .5))
                              #fancybox=True, shadow=True)
        pl.xlabel(self.x_label)
        pl.ylabel(self.y_label)
        pl.savefig(filename,  bbox_inches='tight',dpi=300)

def plot_results_in_dir():

    p = Plot()
    #for folder, __, ___ in os.walk('data/results3/uniform_vs_adaptive_squished_cube/'):
    #for folder, __, ___ in os.walk('output/test_adaptive/'):
    for folder, __, ___ in os.walk(p.indir):
        if p.check_for_results_folder(folder):
            p.plot(folder)
plot_results_in_dir()
#python utils/plot_results.py --xvar="coarse spacing" --yvar "max absolute error" --yvar "total fmm time time (s)" --prefix "qbx__ref_lvl_" --prefix "singular_eval__ref_lvl_" --label "error" --label "FMM time"
