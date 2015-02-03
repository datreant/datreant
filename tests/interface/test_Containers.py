"""Interface tests for Sims and Groups.

"""

import pytest
import os, shutil

class TestSim:
    """Test the Sim object.

    """
    @pytest.fixture(scope='module')
    def sim(request):
        import MDSynthesis as mds
        return mds.Sim('testsim')

    @pytest.mark.parametrize(
    def test_sim_init():

        


