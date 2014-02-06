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
        
        :method:`_run_container_pre()`
        :method:`_run_container_loop()`
        :method:`_run_container_post()`
        :method:`analyze()`

        at minimum to implement specific results.

        :Arguments:
            *args
                Sim (or Sim-derived) objects to analyze

        """
        super(Analysis, self).__init__(*args, **kwargs)

    def _run_container(self, container, **kwargs):
        """Run timeseries collection for single container.

        :Arguments:
            *container*
                Sim (or Sim-derived) object

        :Keywords:
            *universe*
                key of universe to use for analysis [``'main'``]

            **kwargs passed to :method:`_run_container_pre()`
                               :method:`_run_container_loop()`
                               :method:`_run_container_post()`

        """
        container._logger.info("Running {} analysis on '{}'...".format(self.__class__.__name__, container.metadata['name']))

        # set up data storage structure
        con_results = {'time': np.zeros((len(container.universe.trajectory),), dtype=float),
                      }

        self._run_container_pre(container, con_results, **kwargs)

        # iterate through trajectory; collect raw data
        container._logger.info("Collecting timeseries...")
        pm = ProgressMeter(container.universe.trajectory.numframes, interval=100)
        container.universe.trajectory[0]
        for ts in container.universe.trajectory:
            pm.echo(ts.frame)
            con_results['time'][container.universe.trajectory.frame - 1] = container.universe.trajectory.time
            self._run_container_loop(container, con_results, **kwargs)

        self._run_container_post(container, con_results, **kwargs)
        self._save(container, con_results)
 
    def _run_container_pre(container, con_results, **kwargs):
        """Operations to be performed before run loop.

        :Arguments:
            *container*
                Sim (or Sim-derived) object
            *con_results*
                dict storing results for *container*

        :Keywords:

        """
        return

    def _run_container_loop(container, con_results, **kwargs):
        """Operations to be performed inside of run loop.

        :Arguments:
            *container*
                Sim (or Sim-derived) object
            *con_results*
                dict storing results for *container*

        :Keywords:
    
        """
        return

    def _run_container_post(container, con_results, **kwargs):
        """Operations to be performed after run loop.

        :Arguments:
            *container*
                Sim (or Sim-derived) object
            *con_results*
                dict storing results for *container*

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

        :method:`_run_container_pre()`
        :method:`_run_container_loop()`
        :method:`_run_container_post()`
        :method:`analyze()`

        at minimum to implement specific results.

        :Arguments:
            *args
                Group (or Group-derived) objects to analyze
        """
        super(MetaAnalysis, self).__init__(*args, **kwargs)
