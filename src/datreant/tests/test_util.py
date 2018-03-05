import mock
import os
import pytest

import datreant as dtr


def test_makedirs(tmpdir):
    with tmpdir.as_cwd():
        dtr.util.makedirs('this/should/exist', exist_ok=True)
        assert os.path.exists('this/should/exist')
        assert os.path.isdir('this/should/exist')


def test_makedirs_raise(tmpdir):
    with tmpdir.as_cwd():
        dtr.util.makedirs('this/should/exist')
        with pytest.raises(OSError):
            dtr.util.makedirs('this/should/exist')


def test_makedirs_exists(tmpdir):
    # try and make a dir twice
    with tmpdir.as_cwd():
        os.mkdir('this/')
        dtr.util.makedirs('this/', exist_ok=True)
        assert os.path.exists('this/')
        assert os.path.isdir('this/')


def test_makedirs_error_catch(tmpdir):
    # mock a disk full error
    # and make sure it gets propagated through properly
    with tmpdir.as_cwd():
        with mock.patch('os.makedirs') as mp:
            mp.side_effect = OSError(os.errno.ENOSPC, 'Mock - disk full')
            # check the specific error code
            # ie check we don't mangle it enroute
            with pytest.raises(OSError) as error:
                dtr.util.makedirs('this/should/fail', exist_ok=True)
                assert error.errno == os.errno.ENOSPC
