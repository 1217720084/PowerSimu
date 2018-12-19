#

# import subprocess
#
# print('$ nslookup')
# p = subprocess.Popen(['nslookup'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# output, err = p.communicate(b'set q=mx\npython.org\nexit\n')
# print(output.decode('utf-8'))
# print('Exit code:', p.returncode)


# import threading
#
# # 创建全局ThreadLocal对象:
# local_school = threading.local()
#
#
# def process_student():
#     # 获取当前线程关联的student:
#     std = local_school.student
#     print('Hello, %s (in %s)' % (std, threading.current_thread().name))
#
#
# def process_thread(name):
#     # 绑定ThreadLocal的student:
#     local_school.student = name
#     process_student()
#
#
# t1 = threading.Thread(target=process_thread, args=('Alice',), name='Thread-A')
# t2 = threading.Thread(target=process_thread, args=('Bob',), name='Thread-B')
# t1.start()
# t2.start()
# t1.join()
# t2.join()

# 进程共享数据

# import multiprocessing
#
#
# def func(num):
#     num.value = 10.78  # 子进程改变数值的值，主进程跟着改变
#
#
# if __name__ == "__main__":
#     num = multiprocessing.Value("d", 10.0)  # d表示数值,主进程与子进程共享这个value。（主进程与子进程都是用的同一个value）
#     print(num.value)
#
#     p = multiprocessing.Process(target=func, args=(num,))
#     p.start()
#     p.join()
#
#     print(num.value)


# 多进程共享数组
# import multiprocessing
#
#
# def func(num):
#     num[2] = 9999  # 子进程改变数组，主进程跟着改变
#
#
# if __name__ == "__main__":
#     num = multiprocessing.Array("i", [1, 2, 3, 4, 5])  # 主进程与子进程共享这个数组
#     print(num[:])
#
#     p = multiprocessing.Process(target=func, args=(num,))
#     p.start()
#     p.join()
#
#     print(num[:])

# import multiprocessing
# import time
# import os
#
# datalist = ['+++']  # 全局变量，主进程与子进程是并发执行的，他们不能共享全局变量(子进程不能改变主进程中全局变量的值)
#
#
# def adddata():
#     global datalist
#     datalist.append(1)
#     datalist.append(2)
#     datalist.append(3)
#     print("子进程", os.getpid(), datalist)
#
#
# if __name__ == "__main__":
#     p = multiprocessing.Process(target=adddata, args=())
#     p.start()
#     p.join()
#     datalist.append("a")
#     datalist.append("b")
#     datalist.append("c")
#     print("主进程", os.getpid(), datalist)

# # -*- coding:utf-8 -*-
# from multiprocessing import Process, Pipe
#
# def proc1(pipe):
#     for i in range(10):
#        s = ('Hello,This is proc%i' % i)
#        pipe.send(s)
#
# def proc2(pipe):
#     while True:
#         print("proc2 recieve:", pipe.recv())
#
# if __name__ == "__main__":
#     pipe = Pipe()
#     p1 = Process(target=proc1, args=(pipe[0],))
#     p2 = Process(target=proc2, args=(pipe[1],))
#     p3 = Process(target=proc1, args=(pipe[0],))
#     p1.start()
#     p2.start()
#     p3.start()
#     p1.join()
#     p3.join()
#     p2.join(2)   #限制执行时间最多为2秒
#     print('\nend all processes.')

from multiprocessing import Process, Manager
def func(dt, lt):

    for i in range(10):
        key = 'arg' + str(i)
        dt[key] = i * i

    lt += range(11, 16)

def func2(dt, lt):

    for i in range(10):
        key = 'arg' + str(i+10)
        dt[key] = (i+10) ** 2
    lt += range(30, 35)

if __name__ == "__main__":
    manager = Manager()
    dt = manager.dict()
    lt = manager.list()

    p1 = Process(target=func, args=(dt, lt))
    p2 = Process(target=func2, args=(dt, lt))
    p1.start()
    p2.start()
    p1.join()
    p2.join()
    print(dt, '\n', lt)