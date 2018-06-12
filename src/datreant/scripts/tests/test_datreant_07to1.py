import pytest
from glob import glob
from os import path
import json

import datreant as dtr
from datreant.scripts.datreant_07to1 import convert


@pytest.fixture
def old_treant():
    return {
        "categories": {
            "totaltime": 2336232.0,
            "selection": " name CA ",
            "reference": "../../ave_struc_all_ca.pdb",
            "time_steps": 8333,
            "sim-path": "box1/production/traj-2/",
            "protein": "myoglobin",
            "dt": 6.0
        },
        "tags": ["foo", "BAR"]
    }


def test_convert(old_treant, tmpdir):
    treant_folder = str(tmpdir)
    with open(path.join(treant_folder, 'Treant-uuid.json'), 'w') as fh:
        json.dump(old_treant, fh)

    convert(treant_folder)
    treant = dtr.Treant(treant_folder)

    assert treant.tags == old_treant['tags']
    assert treant.categories == old_treant['categories']


def test_convert_no_treant(tmpdir):
    treant_folder = str(tmpdir)
    with pytest.warns(UserWarning, match="No treant found.*"):
        convert(treant_folder)


def test_convert_to_many_treants(tmpdir, old_treant):
    treant_folder = str(tmpdir)
    with open(path.join(treant_folder, 'Treant-uuid.json'), 'w') as fh:
        json.dump(old_treant, fh)

    with open(path.join(treant_folder, 'Treant-uuid2.json'), 'w') as fh:
        json.dump(old_treant, fh)

    with pytest.warns(UserWarning, match="Multiple treants found.*"):
        convert(treant_folder)
