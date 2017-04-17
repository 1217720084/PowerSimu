"""
 测试模型
"""

# import
import system

# Bus
bus_case = {'Vb': 10.5, 'bus': 'Bus_1'}
system.Bus._init_data()
system.Bus.add(**bus_case)
bus_case2 = {'Vb': 220, 'bus': 'Bus_2'}
system.Bus.add(**bus_case2)
bus_case3 = {'Vb': 110, 'bus': 'Bus_3'}
system.Bus.add(**bus_case3)
system.Bus._xy_index()
system.Bus._bus_index()
system.Bus._list2matrix()
system.Bus.yinit(system.DAE)
print('Bus:')
print(system.Bus.__dict__)
print(system.DAE.__dict__)

# PV
PV_case = {'bus': 'Bus_1', 'Sn': 100, 'Vn': 10.5, 'Pg': 3, 'qgmax': 6, 'qgmin': -6, 'V0': 1.05,'Vmax': 1.1, 'Vmin': 0.95}
system.PV._init_data()
system.PV.add(**PV_case)
#PV_case2 = {'bus': 'Bus_2', 'Sn': 100, 'Vn': 220, 'Pg': 4, 'qgmax': 6, 'qgmin': -6, 'V0': 1.05,'Vmax': 1.1, 'Vmin': 0.95}
#system.PV.add(**PV_case2)
#PV_case3 = {'bus': 'Bus_3', 'Sn': 100, 'Vn': 110, 'Pg': 5, 'qgmax': 6, 'qgmin': -6, 'V0': 1.05,'Vmax': 1.1, 'Vmin': 0.95}
#system.PV.add(**PV_case3)
#system.PV._xy_index()
system.PV._bus_index()
system.PV._list2matrix()
system.PV.base(Vb=system.Bus.Vb[system.PV.a])
#system.PV.base(Vb=system.Bus.Vb[system.PV.a[1]])
system.PV.yinit(system.DAE)
print('PV:')
print(system.PV.__dict__)
print(system.DAE.__dict__)

# SW
SW_case = {'bus': 'Bus_3', 'Sn': 100, 'Vn': 110, 'Pg': 3, 'qgmax': 6, 'qgmin': -6, 'V0': 1.05,'Vmax': 1.1, 'Vmin': 0.95, 'Va': 0.0}
system.SW._init_data()
system.SW.add(**SW_case)
system.SW._bus_index()
#system.SW._list2matrix()
#system.SW.base(Vb=system.Bus.Vb[system.SW.a])
#system.PV.base(Vb=system.Bus.Vb[system.PV.a[1]])
system.SW.yinit(system.DAE)
print('SW:')
print(system.SW.__dict__)
print(system.DAE.__dict__)

# PQ
PQ_case = {'bus': 'Bus_2', 'Sn': 100, 'Vn': 220, 'Pl': 3, 'Vmax': 1.1, 'Vmin': 0.95, 'Ql': 1.0}
system.PQ._init_data()
system.PQ.add(**PQ_case)
system.PQ._bus_index()
#system.PQ._list2matrix()
#system.PQ.base(Vb=system.Bus.Vb[system.SW.a])
system.PQ.yinit(system.DAE)
print('PQ:')
print(system.PQ.__dict__)
print(system.DAE.__dict__)

# Shunt
Shunt_case = {'bus': 'Bus_1', 'Sn': 100, 'Vn': 110, 'g': 3.0, 'b': 5.0}
system.Shunt._init_data()
system.Shunt.add(**Shunt_case)
system.Shunt._bus_index()
system.Shunt._list2matrix()
system.Shunt.base(Vb=system.Bus.Vb[system.Shunt.a])
system.Shunt.yinit(system.DAE)
system.Shunt._matrix2list()
print('Shunt:')
print(system.Shunt.__dict__)
print(system.DAE.__dict__)

system.Bus._matrix2list()
print(system.Bus.__dict__)

# Line
# Line_case = {'f_bus': 'Bus_2', 'to_bus': 'Bus_3', 'r': 1, 'x': 2}
# system.Line._init_data()
# system.Line.add(**Line_case)
# print(system.Line.__dict__)
# system.Line._bus_index()
#print('Line:')
# print(system.Line.__dict__)