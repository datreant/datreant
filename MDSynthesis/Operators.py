"""
Operator objects that perform some kind of operation using container objects
as input.

"""
import numpy as np
import os
import sys
import cPickle
from multiprocessing import Process
import matplotlib.pyplot as plt

import MDAnalysis
from MDAnalysis.core.log import ProgressMeter

class Analysis(object):
    """Base class for analysis on individual Sim objects.
        
    """
    
    def __init__(self, *args, **kwargs):
        """Initialize analysis object.

        :Arguments:
            *args
                Base-derived objects to analyze 

        """
        self.systems = list(args)

    def run(self, **kwargs):
        """Obtain timeseries data.

        :Keywords:
            *force*
                If True, force recollection of data; default False

            **kwargs passed to `:meth:self._run_system()`
        """
        joblist = []
        force = kwargs.pop('force', False)

        for system in self.systems:
            if (not self._datacheck(system)) or force:
                p = (Process(target=self._run_system, args=(system,), kwargs=kwargs))
                p.start()
                joblist.append(p)
            else:
                system.logger.info('{} data already present; skipping data collection.'.format(self.__name__))

            # update analysis list in each object
            if not self.__name__ in system.metadata['analysis_list']:
                system.metadata['analysis_list'].append(self.__name__)
                system.save()

        for p in joblist:
            p.join()

    def _run_system(self, system, **kwargs):
        """Run timeseries collection for single system.

        """
        system.logger.info("Running {} analysis on '{}'...".format(self.__name__, system.metadata['name']))

        # set up data storage structure
        sys_results = {'time': np.zeros((len(system.universe.trajectory),), dtype=float),
                      }

        # iterate through trajectory; collect raw data
        system.logger.info("Collecting timeseries...")
        pm = ProgressMeter(system.universe.trajectory.numframes, interval=100)
        system.universe.trajectory[0]
        for ts in system.universe.trajectory:
            pm.echo(ts.frame)
            sys_results['time'][system.universe.trajectory.frame - 1] = system.universe.trajectory.time
            # stuff to do
        
        self.save(system, sys_results)


    def analyze(self, **kwargs):
        """Perform analysis of timeseries.

        """
        # make sure data loaded into each system; should use try/catch here
        self.load()

    def save(self, system, sys_results):
        """Save results to main data file.

        :Arguments:
            *system*
                system to save data for
            *sys_results*
                results for system
        """
        analysis_dir = os.path.join(system.metadata['basedir'], self.__name__)
        system._makedirs(analysis_dir)
        main_file = os.path.join(analysis_dir, '{}.pkl'.format(self.__name__))

        with open(main_file, 'wb') as f:
            cPickle.dump(sys_results, f)

    def load(self, **kwargs):
        """Load data for each system if not already loaded.

        :Keywords:
            *force*
                If True, force reload of data; default False
        """
        force = kwargs.pop('force', False)

        # make sure data loaded into each system; should use try/catch here
        for system in self.systems:
            if (not self.__name__ in system.analysis.keys()) or force:
                system.load(self.__name__)
    
    def _datacheck(self, system):
        """Check if data file already present.

        :Arguments:
            *system*
                Base-derived object

        :Returns:
            *present*
                True if data is already present; False otherwise
        """
        analysis_dir = os.path.join(system.metadata['basedir'], self.__name__)
        main_file = os.path.join(analysis_dir, '{}.pkl'.format(self.__name__))
        return os.path.isfile(main_file)

class AnalysisSet(object):
    """Base class for analysis on SimSet objects.

    """
