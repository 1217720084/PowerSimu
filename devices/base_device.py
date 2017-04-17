"""


"""
# import

from cvxopt.base import mul, div, matrix
import system

class base_device:

    def __init__(self):
        self.n = 0      # 设备的数量
        self.u = []     # 设备的状态
        self.name = []  # 设备的名字
        self.int = {}   # 索引字典（idx:index）

        # 元属性

        self._type = None  # 设备类型
        self._name = None  # 设备名字
        self._bus = {}      # 母线索引
        self._states = []   # 状态变量列表
        self._algebs = []   # 代数变量列表
        self._powers = []   # 功率
        self._voltages = [] # 电压
        self._currents = [] # 电流
        self._z = []        # 阻抗
        self._y = []        # 导纳

        #非零参数
        self._zeros = ['Sn', 'Vn']

        #参数
        self._params = ['u', 'Sn', 'Vn']

        # 数据字典

        self._data = {'u': 1, 'Sn': 100.0, 'Vn': 220.0}

        # 参数单位

        self._units = {'Sn': 'MVA', 'Vn': 'kV', 'u': 'boolean'}

        # 参数描述

        self.descr = {'Sn': '基准功率', 'Vn': '基准电压', 'u': '连接状态'}

#初始化数据

    def _init_data(self):

        for arg in self._data:
            self.__dict__[arg] = []

        for key in self._bus:
            for item in self._bus[key]:
                self.__dict__[item] = []

        for arg in self._states:
            self.__dict__[arg] = []

        for arg in self._algebs:
            self.__dict__[arg] = []

        if self._name is None: self._name = self._type

#把参数化为系统标幺值

    def base(self, Sb = 100.0, Vb = None):

        for var in self._voltages:
            self.__dict__[var] = mul(self.__dict__[var], self.Vn)
            self.__dict__[var] = div(self.__dict__[var], Vb)

        for var in self._powers:
            self.__dict__[var] = mul(self.__dict__[var], self.Sn)
            self.__dict__[var] = div(self.__dict__[var],Sb)

        for var in self._currents:
            self.__dict__[var] = mul(self.__dict__[var], self.Sn)
            self.__dict__[var] = div(self.__dict__[var], self.Vn)
            self.__dict__[var] = mul(self.__dict__[var], Vb)
            self.__dict__[var] /= Sb

        if len(self._z) or len(self._y):
            Zn = div(self.Vn**2, self.Sn)
            #print(Zn)
            Zb = div(Vb**2, Sb)
            #print(Zb)

            for var in self._z:
                self.__dict__[var] = mul(self.__dict__[var], Zn)
                self.__dict__[var] = div(self.__dict__[var], Zb)

            for var in self._y:
                if self.__dict__[var].typecode == 'd':
                    self.__dict__[var] = div(self.__dict__[var], Zn)
                    self.__dict__[var] = mul(self.__dict__[var], Zb)
                    #print(self.__dict__[var])
                if self.__dict__[var].typecode == 'z':
                    self.__dict__[var] = div(self.__dict__[var], Zn + 0j)
                    self.__dict__[var] = mul(self.__dict__[var], Zb + 0j)

    def add(self, idx=None, name=None, **kwargs):

        if idx is None:
            idx = self._type + '_' + str(self.n + 1)

        self.int[idx] = self.n
        self.n += 1

        if name is None:
            self.name.append(self._type + ' ' + str(self.n))
        else:
            self.name.append(name)

        #设置默认值
        for key, value in self._data.items():
            self.__dict__[key].append(value)

        # 根据输入参数kwargs覆盖对应值
        for key, value in kwargs.items():

            if not key in self._data:#self._data.has_key(key):
                print('这个设备没有参数<%s>.' % key)
                continue

            self.__dict__[key][-1] = value

            # 检查非零数据
            if not value and key in self._zeros:
                if key == 'Sn':
                    default = system.Settings.mva
                elif key == 'Vn':
                    default = self._data[key]
                elif key == 'fn':
                    default = system.Settings.freq
                else:
                    default = self._data[key]
                self.__dict__[key][-1] = default

    #
    def _xy_index(self):

        zeros = [0] * self.n
        for item in self._states:
            self.__dict__[item] = zeros[:]
        for item in self._algebs:
            self.__dict__[item] = zeros[:]

        for item in self._states:

            for var in range(self.n):
                self.__dict__[item][var] = system.DAE.nx
                system.DAE.nx += 1

        for item in self._algebs:
            for var in range(self.n):
                self.__dict__[item][var] = system.DAE.ny
                system.DAE.ny += 1

       # for var in range(self.n):

           # for item in self._states:
             #   self.__dict__[item][var] = system.DAE.nx
              #  system.DAE.nx += 1
           # for item in self._algebs:
              #  self.__dict__[item][var] = system.DAE.ny
               # system.DAE.ny += 1

    def _bus_index(self):

        for index in self._bus.keys():
            for item in self.__dict__[index]:
                if not item in system.Bus.int:
                    continue
                    #self.message('Bus index <%s> does not exist', data_tuple = item, level = self.ERROR)
                else:
                    idx = system.Bus.int[item]
                    self.__dict__[self._bus[index][0]].append(system.Bus.a[idx])
                    self.__dict__[self._bus[index][1]].append(system.Bus.v[idx])

    def _list2matrix(self):

        for item in self._params:
            self.__dict__[item] = matrix(self.__dict__[item])

    def _matrix2list(self):

        for item in self._params:
            self.__dict__[item] = list(self.__dict__[item])

    # 删除某一个元件
    def remove(self, idx = None):
        return
