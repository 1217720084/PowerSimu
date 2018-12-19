import system


class device():

    def __init__(self, device_list):
        self.devices = []  # 当前仿真系统含有的模型
        self.models = device_list  # 软件所含有的所有模型
        self.n = 0
        self.gcall = []
        self.Gycall = []
        self.fcall = []
        self.Fxcall = []
        self.windup = []
        self.pflow = []
        self.xinit = []
        self.shunt = []
        self.series = []
        self.flows = []
        self.connection = []
        self.times = []
        self.stagen = []
        self.dyngen = []
        self.gmcall = []
        self.fmcall = []
        self.dcseries = []
        self.opf = []
        self.obj = []

    def setup(self):
        self.n = 0
        for item in self.models:
            if system.__dict__[item].n:
                self.n += 1
                self.devices.append(item)
                properties = system.__dict__[item].properties
                for key in properties.keys():
                    self.__dict__[key].append(properties[key])

        # 写call函数
        string = '"""\n'
        for gcall, device in zip(self.gcall, self.devices):
            if gcall: string += 'system.' + device + '.gcall()\n'
        string += '\n'
        for gycall, device in zip(self.Gycall, self.devices):
            if gycall: string += 'system.' + device + '.gycall()\n'
        string += '\n'
        for fcall, device in zip(self.fcall, self.devices):
            if fcall: string += 'system.' + device + '.fcall()\n'
        string += '\n'
        for fxcall, device in zip(self.Fxcall, self.devices):
            if fxcall: string += 'system.' + device + '.Fxcall()\n'
        string += '"""'
        self.call_int = compile(eval(string), '', 'exec')

        string = '"""\n'
        for gcall, device in zip(self.gcall, self.devices):
            if gcall: string += 'system.' + device + '.gcall()\n'
        string += '\n'
        for gycall, device in zip(self.Gycall, self.devices):
            if gycall: string += 'system.' + device + '.Gycall()\n'
        string += '\n'

        string += '"""'
        self.call_pf = compile(eval(string), '', 'exec')

