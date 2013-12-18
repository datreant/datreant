"""
Lower level mixins.

"""
import os
import yaml
import sys
import cPickle
import logging

class ContainerCore(object):
    """Mixin class for all containers.

    The ContainerCore object is not intended to be useful on its own, but
    instead contains methods and attributes common to all Container objects.

    """
    def __init__(self):
        """
        
        """
        self.metadata = dict()              # information about object; defines base object
        self.analysis = dict()              # analysis data 'modular dock'
    
    def save(self):
        """Save base object metadata.

        """
        basedir = self._rel2abspath(self.metadata["basedir"])
        self._makedirs(basedir)

        with open(os.path.join(basedir, self.metadata['metafile']), 'w') as f:
            yaml.dump(self.metadata, f)

    def refresh(self):
        """Reloads metadata from file.

        """
        basedir = self._rel2abspath(self.metadata['basedir'])
        metafile = os.path.join(basedir, self.metadata['metafile'])
        with open(metafile, 'r') as f:
            self.metadata = yaml.load(f)

    def load(self, *args):
        """Load data instances into object.

        If 'all' is in argument list, every available dataset is loaded.

        :Arguments:
            *args*
                datasets to load as attributes
        """
        if 'all' in args:
            self._logger.info("Loading all known data into object '{}'...".format(self.metadata['name']))
            for i in self.metadata['analysis_list']:
                self._logger.info("Loading {}...".format(i))
                with open(os.path.join(self._rel2abspath(self.metadata['basedir']), '{}/{}.pkl'.format(i, i)), 'rb') as f:
                    self.analysis[i] = cPickle.load(f)
            self._logger.info("Object '{}' loaded with all known data.".format(self.metadata['name']))
        else:
            self._logger.info("Loading selected data into object '{}'...".format(self.metadata['name']))
            for i in args:
                self._logger.info("Loading {}...".format(i))
                with open(os.path.join(self._rel2abspath(self.metadata['basedir']), '{}/{}.pkl'.format(i, i)), 'rb') as f:
                    self.analysis[i] = cPickle.load(f)
            self._logger.info("Object '{}' loaded with selected data.".format(self.metadata['name']))

    def unload(self, *args):
        """Unload data instances from object.

        If 'all' is in argument list, every loaded dataset is unloaded.

        :Arguments:
            *args*
                datasets to unload
        """
        if 'all' in args:
            self.analysis.clear()
            self._logger.info("Object '{}' unloaded of all data.".format(self.metadata['name']))
        else:
            self._logger.info("Unloading selected data from object {}...".format(self.metadata['name']))
            for i in args:
                self._logger.info("Unloading {}...".format(i))
                self.analysis.pop(i, None)
            self._logger.info("Object '{}' unloaded of selected data.".format(self.metadata['name']))

    def _abs2relpath(self, abspath):
        """Return path to file relative to project path.
        
        """
        return abspath.replace(self.metadata['projectdir'], '$PROJECT')

    def _rel2abspath(self, relpath):
        """Return realpath given a path relative to project directory. The
            opposite of _project_relpath.
        """
        return relpath.replace('$PROJECT', self.metadata['projectdir'])
    
    def _update_projectdir(self, basedir_abs):
        """Update projectdir based on given (and current) absolute path of basedir.

        """
        self.metadata['projectdir'] = basedir_abs.partition('/MDSynthesis')[0]
        
    def _makedirs(self, p):
        if not os.path.exists(p):
            os.makedirs(p)
    
    def _build_metadata(self, **kwargs):
        """Build metadata. Runs each time object is generated.
        
        Only adds keys; never modifies existing ones.

        :Keywords:
            *name*
                desired name of object, used for logging and referring to
                object in some analyses; default is trajectory file directory
                basename
        """
        # building core items
        attributes = {'name': kwargs.pop('name', self.__class__.__name__),
                      'logfile': '{}.log'.format(self.__class__.__name__),
                      'analysis_list': [],
                      'type': self.__class__.__name__,
                      }

        for key in attributes.keys():
            if not key in self.metadata:
                self.metadata[key] = attributes[key]
    
    def _build_attributes(self):
        """Build attributes. Needed each time object is generated.

        """
        self._start_logger()
    
    def _start_logger(self):
        """Start up the logger.

        """
        # set up logging
        self._logger = logging.getLogger('{}.{}'.format(self.__class__.__name__, self.metadata['name']))
        self._logger.setLevel(logging.INFO)

        # file handler
        logfile = self._rel2abspath(os.path.join(self.metadata['basedir'], self.metadata['logfile']))
        fh = logging.FileHandler(logfile)
        ff = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
        fh.setFormatter(ff)
        self._logger.addHandler(fh)

        # output handler
        ch = logging.StreamHandler(sys.stdout)
        cf = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
        ch.setFormatter(cf)
        self._logger.addHandler(ch)
