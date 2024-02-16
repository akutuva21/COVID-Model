import numpy as np
import matplotlib.pyplot as plt

# irf3 = 1, irf9 = 0, stat1 = 0, stat3 = 1, nfkb (rela) = 1, foxo = 1, creb1 = 1
a = np.array([1, 0, 0, 1, 1, 1, 1])
# considers foxo is 1
b = np.random.randint(0, 2, size=(1000,7))
s = np.sum((b - a)**2, axis = 1)
h = plt.hist(s)
plt.show()