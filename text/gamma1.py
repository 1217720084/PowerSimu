from scipy.stats import gamma
import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots(1, 1)

a = 0.01
mean, var, skew, kurt = gamma.stats(a, moments='mvsk')
x = np.linspace(gamma.ppf(0.01, a),
                gamma.ppf(0.99, a), 100)
ax.plot(x, gamma.pdf(x, a, loc=0.0, scale=2.0),
       'r-', lw=5, alpha=0.6, label='gamma pdf')
# ax.plot(x, gamma.pdf(x, a),
#        'r-', lw=5, alpha=0.6, label='gamma pdf')


rv = gamma(a, loc=0.0, scale=2.0)
# rv = gamma(a)
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')

vals = gamma.ppf([0.001, 0.5, 0.999], a)
np.allclose([0.001, 0.5, 0.999], gamma.cdf(vals, a))

r = gamma.rvs(a, size=1000, loc=0.0, scale=2.0)
# r = gamma.rvs(a, size=1000)

# ax.hist(r, density='True', histtype='stepfilled', alpha=0.2)
# ax.legend(loc='best', frameon=False)
plt.show()