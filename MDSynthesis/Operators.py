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
        """Generate analysis object.

        Analysis objects perform an analysis operation on individual Sims,
        looping through the trajectory to look at each timestep.

        When multiple Sims are loaded, each Analysis method operates on them
        in parallel using separate processes.

        Important methods::

        :method:`run()` - perform this time-consuming trajectory-loop data
                          collection

        :method:`analyze()` - post-collection operations that do not directly
                              require trajectory information

        This class is meant to be a parent class to a specific analysis for
        specific Sim-derived objects. When inheriting, one need only replace::
        
        :method:`_run_system_pre()`
        :method:`_run_system_loop()`
        :method:`_run_system_post()`
        :method:`analyze()`

        at minimum to implement specific results.

        :Arguments:
            *args
                Sim (or Sim-derived) objects to analyze

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
                system.logger.info('{} data already present; skipping data collection.'.format(self.__class__.__name__))

            # update analysis list in each object
            if not self.__class__.__name__ in system.metadata['analysis_list']:
                system.metadata['analysis_list'].append(self.__class__.__name__)
                system.save()

        for p in joblist:
            p.join()

    def _run_system(self, system, **kwargs):
        """Run timeseries collection for single system.

        :Arguments:
            *system*
                Sim (or Sim-derived) object

        :Keywords:
            **kwargs passed to :method:`_run_system_pre()`
                               :method:`_run_system_loop()`
                               :method:`_run_system_post()`

        """
        system.logger.info("Running {} analysis on '{}'...".format(self.__class__.__name__, system.metadata['name']))

        # set up data storage structure
        sys_results = {'time': np.zeros((len(system.universe.trajectory),), dtype=float),
                      }

        self._run_system_pre(system, sys_results, **kwargs)

        # iterate through trajectory; collect raw data
        system.logger.info("Collecting timeseries...")
        pm = ProgressMeter(system.universe.trajectory.numframes, interval=100)
        system.universe.trajectory[0]
        for ts in system.universe.trajectory:
            pm.echo(ts.frame)
            sys_results['time'][system.universe.trajectory.frame - 1] = system.universe.trajectory.time
            self._run_system_loop(system, sys_results, **kwargs)

        self._run_system_post(system, sys_results, **kwargs)
        self._save(system, sys_results)
 
    def _run_system_pre(system, sys_results, **kwargs):
        """Operations to be performed before run loop.

        :Arguments:
            *system*
                Sim (or Sim-derived) object
            *sys_results*
                dict storing results for *system*

        :Keywords:

        """
        return

    def _run_system_loop(system, sys_results, **kwargs):
        """Operations to be performed inside of run loop.

        :Arguments:
            *system*
                Sim (or Sim-derived) object
            *sys_results*
                dict storing results for *system*

        :Keywords:
    
        """
        return

    def _run_system_post(system, sys_results, **kwargs):
        """Operations to be performed after run loop.

        :Arguments:
            *system*
                Sim (or Sim-derived) object
            *sys_results*
                dict storing results for *system*

        :Keywords:

        """
        return

    def analyze(self, **kwargs):
        """Perform analysis of timeseries.

        Does not require stepping through trajectory.

        """
        # make sure data loaded into each system; should use try/catch here
        self._load()

    def _save(self, system, sys_results):
        """Save results to main data file.

        :Arguments:
            *system*
                system to save data for
            *sys_results*
                results for system
        """
        analysis_dir = os.path.join(system.metadata['basedir'], self.__class__.__name__)
        system._makedirs(analysis_dir)
        main_file = os.path.join(analysis_dir, '{}.pkl'.format(self.__class__.__name__))

        with open(main_file, 'wb') as f:
            cPickle.dump(sys_results, f)

    def _load(self, **kwargs):
        """Load data for each system if not already loaded.

        :Keywords:
            *force*
                If True, force reload of data; default False
        """
        force = kwargs.pop('force', False)

        # make sure data loaded into each system; should use try/catch here
        for system in self.systems:
            if (not self.__class__.__name__ in system.analysis.keys()) or force:
                system.load(self.__class__.__name__)
    
    def _datacheck(self, system):
        """Check if data file already present.

        :Arguments:
            *system*
                Base-derived object

        :Returns:
            *present*
                True if data is already present; False otherwise
        """
        analysis_dir = os.path.join(system.metadata['basedir'], self.__class__.__name__)
        main_file = os.path.join(analysis_dir, '{}.pkl'.format(self.__class__.__name__))
        return os.path.isfile(main_file)

class AnalysisSet(object):
    """Base class for analysis on SimSet objects.

    """
