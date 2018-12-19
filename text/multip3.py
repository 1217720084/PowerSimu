# from multiprocessing import Pool
# import time
# import os
#
#
# class A:
#
#     def __init__(self):
#         self.name = 'syn'
#         self.n = 0
#
#     def addsyn(self):
#         self.n = self.n + 1
#
#
# def task(name):
#
#     print('进程%s(%s)开始时间为：%s s' % (name, os.getpid(), time.time()))
#     syn = A()
#     syn.addsyn()
#     print('syn数值为：%s' % syn.n)
#     # print(time.time())
#     # time.sleep(3)
#
#
# if __name__ == '__main__':
#     print('Parent process %s' % os.getpid())
#     print(time.time())
#
#     p = Pool()
#     for i in range(9):
#         p.apply_async(task, args=(i,))
#     print('Waiting for all subprocess done ...')
#     p.close()
#     ts = time.time()
#     p.join()
#     te = time.time()
#     print('多进程结束时间为：%s' % te)
#     print('多进程开断需要耗时： %s' % (te-ts))
#     print('All subprocess done')
#     ts = time.time()
#     task(1)
#     te = time.time()
#     print('进程1结束时间为 %s s' % te)
#     print('单进程运行需要耗时：%s s' % (te-ts))

import multiprocessing

# def foo(q):
#     q.put([11,'hello',True])
#     print(q.qsize())
#
# q=multiprocessing.Queue() #全局定义一个q进程队列，在产生子进程时候会在子进程里生成，可以指定最大数，限制队列长度
# if __name__ == '__main__':
#     p=multiprocessing.Process(target=foo,args=(q, )) #因为名称空间不同，子进程的主线程创建的q队列，主进程get不到，所以会阻塞住
#     p.start()
#     # foo() #主进程执行一下函数就可以访问到了
#     print(q.get())
#     p.join()

# from multiprocessing import Manager,Process
# def foo(l,i):
#     l.append(i**i)
# if __name__ == '__main__':
#     man=Manager()
#     ml=man.list([11,22,33])
#     l=[]
#     for i in range(5):
#         p=Process(target=foo,args=(ml,i))
#         p.start()
#         l.append(p)
#     for i in l: #必须要join，不然会执行报错，处理一个数据必须要一个个来，不能同时处理一个数据
#         i.join()
#     print(ml)

from multiprocessing import Pool
import time

def foo(n):
    print(n)
    time.sleep(1)

if __name__ == '__main__':
    pool_obj=Pool(8)    #
    for i in range(47):
        pool_obj.apply_async(func=foo,args=(i,))
        # pool_obj.apply(func=foo,args=(i,))    #子进程的生成是靠进程池对象维护的
        # apply同步，子进程一个个执行
        # apply_async异步，多个子进程一起执行
    pool_obj.close()
    pool_obj.join()
    print('ending')

