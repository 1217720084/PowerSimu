"""

"""
import system
from devices.base_device import base_device

class syn2():
    def __init__(self):
        return


class syn3():
    def __init__(self):
        return


class syn4():
    def __init__(self):
        return


class syn5a():
    def __init__(self):
        return


class syn5b():
    def __init__(self):
        return


class syn5c():
    def __init__(self):
        return


class syn5d():
    def __init__(self):
        return


class syn6a(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({'bus': None, 'fn': 50, 'm_model': 6, 'xl': 0, 'ra': 0, 'xd': 0, 'xd1': 0, 'xd2': 0,
                          'Td01': 0, 'Td02': 0, 'xq': 0, 'xq1': 0, 'xq2': 0, 'Tq01': 0, 'Tq02': 0, 'M': 0, 'D': 0})
        self._type = 'Syn6'
        self._name = 'Syn6'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['Va', 'V0']
        self._states = ['delta', 'omega', 'e1q', 'e1d', 'e2q', 'e2d']
        self._params.extend(['fn', 'xl', 'ra', 'xd', 'xd1', 'xd2', 'Td01', 'Td02',
                             'xq', 'xq1', 'xq2', 'Tq01', 'Tq02', 'M', 'D'])
        self._voltages = ['V0']  # 后续得修改
        self._powers = ['Pg']    # 后续得修改




class syn6b():
    def __init__(self):
        return

class syn6(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update(
            {'bus': None, 'fn': 50, 'm_model': 6, 'xl': 0, 'ra': 0, 'xd': 0, 'xd1': 0, 'xd2': 0,
             'Td01': 0, 'Td02': 0, 'xq': 0, 'xq1': 0, 'xq2': 0, 'Tq01': 0, 'Tq02': 0, 'M': 0,
             'D': 0, 'pm0': 0, 'vf0': 0, 'J11': 0, 'J12': 0, 'J21': 0, 'J22': 0})
        self._type = 'Syn6'
        self._name = 'Syn6'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['Va', 'V0']
        self._states = ['delta', 'omega', 'e1q', 'e1d', 'e2q', 'e2d']
        self._params.extend(['fn', 'xl', 'ra', 'xd', 'xd1', 'xd2', 'Td01', 'Td02',
                             'xq', 'xq1', 'xq2', 'Tq01', 'Tq02', 'M', 'D'])
        self._voltages = ['V0']  # 后续得修改
        self._powers = ['Pg']    # 后续得修改