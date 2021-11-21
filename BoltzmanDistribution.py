import numpy as np
import matplotlib.pyplot as plt

N=int(200)
Censemble=np.full(N,fill_value=1)

Maxepoch=70
rng=np.random.default_rng() 

for i in range(Maxepoch):
    index=rng.integers(low=0,high=N,size=2)
    Censemble[index[1]]+=Censemble[index[0]]
    Censemble[index[0]]-=Censemble[index[0]]

plt.hist(Censemble,bins=max(Censemble))
plt.show()
