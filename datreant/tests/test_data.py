"""Data structures for storage and retrieval.

"""
import datreant as dtr
import pandas as pd
import numpy as np
import pytest
import os
import shutil
import py


class Series():
    @pytest.fixture
    def datastruct(self):
        data = np.random.rand(10000)
        return pd.Series(data)


class DataFrame():
    @pytest.fixture
    def datastruct(self):
        data = np.random.rand(10000, 3)
        return pd.DataFrame(data, columns=('A', 'B', 'C'))


class Blank_DataFrame():
    @pytest.fixture
    def datastruct(self):
        return pd.DataFrame(np.zeros((10, 10)))


class Wide_Blank_DataFrame():
    @pytest.fixture
    def datastruct(self):
        return pd.DataFrame(np.zeros((1, 10)))


class Thin_Blank_DataFrame():
    @pytest.fixture
    def datastruct(self):
        return pd.DataFrame(np.zeros((10, 1)))


class Panel():
    @pytest.fixture
    def datastruct(self):
        data = np.random.rand(4, 10000, 3)
        return pd.Panel(data, items=('I', 'II', 'III', 'IV'),
                        minor_axis=('A', 'B', 'C'))


class Panel4D():
    @pytest.fixture
    def datastruct(self):
        data = np.random.rand(2, 4, 10000, 3)
        return pd.Panel4D(data, labels=('gallahad', 'lancelot'),
                          items=('I', 'II', 'III', 'IV'),
                          minor_axis=('A', 'B', 'C'))


class NumpyScalar():
    @pytest.fixture
    def datastruct(self):
        return np.array(20)


class Numpy1D():
    @pytest.fixture
    def datastruct(self):
        return np.random.rand(10000)


class Numpy2D():
    @pytest.fixture
    def datastruct(self):
        return np.random.rand(10000, 500)


class Wide_Numpy2D():
    @pytest.fixture
    def datastruct(self):
        return np.zeros((1, 10))


class Thin_Numpy2D():
    @pytest.fixture
    def datastruct(self):
        return np.zeros((10, 1))


class Numpy3D():
    @pytest.fixture
    def datastruct(self):
        return np.random.rand(4, 10000, 45)


class Numpy4D():
    @pytest.fixture
    def datastruct(self):
        return np.random.rand(2, 4, 10000, 45)


class List():
    @pytest.fixture
    def datastruct(self):
        return ['galahad', 'lancelot', 42, 'arthur', 3.14159]


class Dict():
    @pytest.fixture
    def datastruct(self):
        return {'pure': 'galahad', 'brave': 'lancelot', 'answer': 42,
                'king': 'arthur', 'pi-ish': 3.14159}


class Tuple():
    @pytest.fixture
    def datastruct(self):
        return ('arthur', 3.14159)


class Set():
    @pytest.fixture
    def datastruct(self):
        return {'arthur', 3.14159, 'seahorses'}


class Dict_Mix():
    @pytest.fixture
    def datastruct(self):
        return {'an array': np.random.rand(100, 46),
                'another': np.random.rand(3, 45, 2)}
