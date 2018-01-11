import numpy as np

a = np.zeros(10000000)

def modify_a(indexes, value=1):
    a[indexes] += value


modify_a(range(9000000))

#for index in range(9000000):
#    modify_a(index)