from multiprocessing import Pool
import time
import os


def task(name):

    print('进程%s(%s)开始时间为：%s s' % (name, os.getpid(), time.time()))
    # print(time.time())
    # time.sleep(3)


if __name__ == '__main__':
    print('Parent process %s' % os.getpid())
    print(time.time())
    p = Pool()
    for i in range(9):
        p.apply_async(task, args=(i,))
    print('Waiting for all subprocess done ...')
    p.close()
    ts = time.time()
    p.join()
    te = time.time()
    print('多进程结束时间为：%s' % te)
    print('多进程开断需要耗时： %s' % (te-ts))
    print('All subprocess done')
    ts = time.time()
    task(1)
    te = time.time()
    print('进程1结束时间为 %s s' % te)
    print('单进程运行需要耗时：%s s' % (te-ts))
