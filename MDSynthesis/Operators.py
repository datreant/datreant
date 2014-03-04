"""
Operator objects that perform some kind of operation using Container objects
as input.

"""
import numpy as np

import MDAnalysis
from MDAnalysis.core.log import ProgressMeter

from Core import ObjectCore

class Decorators(object):
    """Decorators for use by Operators.

    """
    def parallel(self, f):


    def series(self, f):
        """Apply the method, f, to all Containers loaded into Operator in series.

        """


class OperatorCore(ObjectCore):
    """Core class for all Operators.

    The OperatorCore object is not intended to be useful on its own, but
    instead contains methods and attributes common to all Operator objects.

    """
    _datafile = datafile

    def __init__(self, *args, **kwargs):
        """
        
        """
        super(OperatorCore, self).__init__()
        self.containers = list(args)

    def run(self, **kwargs):
        """Obtain compute-intensive data, usually timeseries.

        This operation is performed on each Container in series by default. 
        Use the keyword *parallel* to operate on each Container as a separate
        process.

        kwargs passed to `:meth:self._run_container()`

        :Keywords:
            *force*
                If True, force recollection of data [``False``]
            *parallel*
                If True, operate on each Container in parallel as 
                separate processes [``False``]
        """
        joblist = []
        force = kwargs.pop('force', False)
        parallel = kwargs.pop('parallel', False)

        # run analysis on each container as a separate process
        if parallel:
            for container in self.containers:
                if (not self._datacheck(container)) or force:
                    p = (Process(target=self._run_container, args=(container,), kwargs=kwargs))
                    p.start()
                    joblist.append(p)
                else:
                    container._logger.info('{} data already present; skipping data collection.'.format(self.__class__.__name__))

            for p in joblist:
                p.join()
        else:
            for container in self.containers:
                self._run_container(container, **kwargs)
    
        # finish up
        for container in self.containers:
            # update analysis list in each object
            if not self.__class__.__name__ in container.metadata['data']:
                container.metadata['data'].append(self.__class__.__name__)
                container.save()

    def analyze(self, **kwargs):
        """Perform analysis of compute-intensive data.

        Does not require stepping through any trajectories.

        """
        # make sure data loaded into each container; should use try/catch here
        self._load()

    def _class(self):
        """Return class name.

        """
        return self.__class__.__name__

    def _save(self, container, cont_results):
        """Save results to main data file.

        :Arguments:
            *container*
                container to save data for
            *cont_results*
                results for container
        """
        outputdir = self._make_savedir(container)
        datafile = self._datapath(container)

        with self.util.open(datafile, 'wb') as f:
            cPickle.dump(cont_results, f)

    def _make_savedir(self, container):
        """Make directory where all output files are placed.

        :Arguments:
            *container*
                Container for which to save data 

        :Returns:
            *outputdir*
                full path to output directory

        """
        outputdir = self._outputdir(container)
        self.util.makedirs(outputdir)
        return outputdir
    
    def _load(self, **kwargs):
        """Load data for each container if not already loaded.

        :Keywords:
            *force*
                If True, force reload of data; default False
        """
        force = kwargs.pop('force', False)

        # make sure data loaded into each system; should use try/catch here
        for container in self.containers:
            if (not self.__class__.__name__ in container.data.keys()) or force:
                container.load(self.__class__.__name__)
    
    def _datacheck(self, container):
        """Check if data file already present.

        :Arguments:
            *container*
                Container object to check
                
        :Returns:
            *present*
                True if data is already present; False otherwise
        """
        return os.path.isfile(self._datapath(container))

    def _outputdir(self, container):
        """Return path to output directory for a particular Container.

        :Arguments:
            *container*
                Container object

        :Returns:
            *analysis_path*
                path to output directory
        """
        return os.path.join(container.metadata['basedir'], self.__class__.__name__)
    
    def _datapath(self, container):
        """Return path to main output datafile for a particular Container.

        :Arguments:
            *container*
                Container object

        :Returns:
            *datafile_path*
                path to datafile
        """
        return os.path.join(self._outputdir(container), self._datafile)

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
        ukey = kwargs.pop('universe', 'main')
        container._logger.info("Running {} analysis on '{}'...".format(self.__class__.__name__, container.metadata['name']))

        container.attach(ukey)

        # set up data storage structure
        con_results = {'time': np.zeros((len(container.universes[ukey].trajectory),), dtype=float),
                      }

        self._run_container_pre(container, con_results, **kwargs)

        # iterate through trajectory; collect raw data
        container._logger.info("Collecting timeseries...")
        pm = ProgressMeter(container.universes[ukey].trajectory.numframes, interval=100)
        container.universes[ukey].trajectory[0]
        for ts in container.universes[ukey].trajectory:
            pm.echo(ts.frame)
            con_results['time'][container.universes[ukey].trajectory.frame - 1] = container.universes[ukey].trajectory.time
            self._run_container_loop(container, con_results, **kwargs)

        self._run_container_post(container, con_results, **kwargs)
        self._save(container, con_results)
 
    def _run_container_pre(self, container, con_results, **kwargs):
        """Operations to be performed before run loop.

        :Arguments:
            *container*
                Sim (or Sim-derived) object
            *con_results*
                dict storing results for *container*

        :Keywords:

        """
        return

    def _run_container_loop(self, container, con_results, **kwargs):
        """Operations to be performed inside of run loop.

        :Arguments:
            *container*
                Sim (or Sim-derived) object
            *con_results*
                dict storing results for *container*

        :Keywords:
    
        """
        return

    def _run_container_post(self, container, con_results, **kwargs):
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
