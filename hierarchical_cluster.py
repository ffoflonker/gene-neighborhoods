#!/usr/bin/python
import sys
import numpy as np
from scipy.cluster import hierarchy
from matplotlib import pyplot as plt


file1 = sys.argv[-2]
file2 = sys.argv[-1]

X = np.loadtxt(file1)

Z = hierarchy.linkage(X, 'average')
leaf=hierarchy.leaves_list(Z)

#print (leaforder)
#fig = plt.figure(figsize=(10, 10))
#dn = hierarchy.dendrogram(Z, orientation='left')
#plt.show()

transX= np.transpose(X)
Z2 = hierarchy.linkage(transX, 'average')
leaf2=hierarchy.leaves_list(Z2)


#fig2 = plt.figure(figsize=(10, 10))
#dn2 = hierarchy.dendrogram(Z2)
#plt.show()

Y = np.loadtxt(file2)
print (Y)
ordercol= Y[:,leaf2]

orderrow= ordercol[leaf,:]

np.savetxt("ordered.tmp.txt",orderrow, delimiter ='\t', fmt= '%i')
np.savetxt("order_spec.tmp.txt",leaf, delimiter ='\t', fmt= '%i')
np.savetxt("order_col.tmp.txt",leaf2, delimiter ='\t', fmt= '%i')