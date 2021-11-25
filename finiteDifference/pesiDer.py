
import numpy as np
import math

def pesiDerUniform(h:float,point : list, derOrder: int)-> np.ndarray:

    """
    Calc finite difference coefficient over a uniform stencil
    
    Parameters:
    - h: stencil spacing
    - point: list of stencil point
    - derOrder: order of the derivative to approxiamate

    Return:
    - weights: finite difference coefficient ordered as the point list

    Example:

    - First derivative central scheme
    >>> w= pesiDerU(1,[-1,+1],1)
    
    """

    #sanity check
    if derOrder >= len(point): 
        raise RuntimeError("derivative order greate than number of stencil points")
    # if len(point) >= derOrder + 2:
    #     raise RuntimeError("reduce stencil: H matrix is singular")

    H=np.zeros((len(point),len(point)))
    

    H_elRows = np.array(point)*h

    for k in range(0,len(point)):
        H[k,:]= (H_elRows**k)/ math.factorial(k)

    Delta= np.zeros((len(point),1))
    Delta[derOrder]= 1

    weight= np.linalg.solve(H,Delta)

    return weight.T

if __name__ == "__main__":
    w= pesiDerUniform(1,[-1,1],1)
    print(w)