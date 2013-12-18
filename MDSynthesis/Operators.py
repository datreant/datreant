"""
Operator objects that perform some kind of operation using Container objects
as input.

"""
import numpy as np

import MDAnalysis
from MDAnalysis.core.log import ProgressMeter

from Core import *

class Analysis(OperatorCore):
    """Base class for analysis on individual Sim objects.
        
    """
    
    def __init__(self, *args, **kwargs):
        """Generate analysis object.

        Analysis objects perform an analysis operation on individual Sims,
        looping through the trajectory to look at each timestep.

        When multiple Sims are loaded, each Analysis method operates on them
        in parallel using separate processes.

        Important methods::

        :method:`run()` - perform time-consuming trajectory-loop data
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
        super(Analysis, self).__init__(*args, **kwargs)

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
        system._logger.info("Running {} analysis on '{}'...".format(self.__class__.__name__, system.metadata['name']))

        # set up data storage structure
        sys_results = {'time': np.zeros((len(system.universe.trajectory),), dtype=float),
                      }

        self._run_system_pre(system, sys_results, **kwargs)

        # iterate through trajectory; collect raw data
        system._logger.info("Collecting timeseries...")
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

class MetaAnalysis(OperatorCore):
    """Base class for analysis on Group objects.

    """
    def __init__(self, *args, **kwargs):
        """Generate meta-analysis object.

        MetaAnalysis objects perform analysis operations on Groups, allowing
        for direct comparison of trajectories to distill larger patterns.

        When multiple Groups are loaded, each MetaAnalysis method operates on
        them in parallel using separate processes.

        Important methods::
    
        :method:`run()` - perform time-consuming trajectory-loop data
                          collection

        :method:`analyze()` - post-collection operations that do not directly
                              require trajectory information

        This class is meant to be a parent class to a specific analysis for
        specific Group-derived objects. When inheriting, one need only replace::

        :method:`_run_system_pre()`
        :method:`_run_system_loop()`
        :method:`_run_system_post()`
        :method:`analyze()`

        at minimum to implement specific results.

        :Arguments:
            *args
                Group (or Group-derived) objects to analyze
        """
        super(MetaAnalysis, self).__init__(*args, **kwargs)
