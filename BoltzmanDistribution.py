import numpy as np
import matplotlib.pyplot as plt
"""
exsercise 4.9 ,page 46 
"Blundell, Stephen J., and Katherine M. Blundell. Concepts in Thermal Physics. 2nd ed. Oxford: Oxford University Press, 2009. Oxford Scholarship Online, 2010. doi: 10.1093/acprof:oso/9780199562091.001.0001."

"""
N=int(400)
Censemble=np.full(N,fill_value=1)

Maxepoch=1000
rng=np.random.default_rng() 

for i in range(Maxepoch):
    index=rng.integers(low=0,high=N,size=2)
    Censemble[index[1]]+=Censemble[index[0]]
    Censemble[index[0]]-=Censemble[index[0]]

plt.hist(Censemble,bins=max(Censemble))
plt.show()
