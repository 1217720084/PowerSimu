import numpy as np
import matplotlib.pyplot as plt
import time
class OrnsteinUhlenbeckActionNoise:
    def __init__(self, mu, sigma=0.3, theta=1, dt=1e-2, x0=None):
        self.theta = theta
        self.mu = mu
        self.sigma = sigma
        self.dt = dt
        self.x0 = x0
        self.reset()

    def __call__(self):
        ts = time.time()
        x = self.x_prev + self.theta * (self.mu - self.x_prev) * self.dt + \
            self.sigma * np.sqrt(self.dt) * np.random.normal(size=self.mu.shape)
        self.x_prev = x
        te = time.time()
        print('一次耗时：%s' % (te-ts))
        return x

    def reset(self):
        self.x_prev = self.x0 if self.x0 is not None else np.zeros_like(self.mu)

    def __repr__(self):
        return 'OrnsteinUhlenbeckActionNoise(mu={}, sigma={})'.format(self.mu, self.sigma)


def main():
    ts = time.time()
    ou_noise = OrnsteinUhlenbeckActionNoise(mu=np.zeros(2))
    plt.figure('data')
    y = []
    t = np.linspace(0, 2, 200)
    # np.random.seed(1)
    for _ in t:
        y.append(ou_noise())
    te = time.time()
    print('耗时：%s s' % (te - ts))
    plt.plot(t, y)
    plt.show()


if __name__ == "__main__":

    main()
