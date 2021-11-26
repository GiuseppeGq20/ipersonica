
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
        raise RuntimeError("derivative order greater than number of stencil points")
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



def pesiDer(x: np.ndarray,xc: float,derOrder=1)-> np.ndarray: # x stencil, xc punto di derivazione, p ordine di derivata
    
    """
    Calc finite difference coefficient over a arbitrary stencil
    
    Parameters:
    - x: stencil 
    - xc: point of derivative collocation
    - derOrder: order of the derivative to approxiamate

    Return:
    - weights: finite difference coefficient ordered as the point list

    Example:

    - First derivative central scheme
    >>> x=np.linspace(0,1,3)
    >>> xc=x[1] 
    >>> w= pesiDerU(1,[-1,+1],1)
    
    """
    n=np.size(x)
    H=np.zeros([n,n])
    
    Delta=np.zeros((n,1))
    Delta[derOrder]=1.0

    for i in range(n):
        H[i,:]=((x-xc)**(i))/math.factorial(i)

    weight= np.linalg.solve(H,Delta)

    return weight.T



if __name__ == "__main__":
    w= pesiDerUniform(1/2,[-1,0,1],1)
    print(w)
    
    x=np.linspace(0,1,3)
    xc=x[1] 
    w= pesiDer(x,xc,1)
    print("\n\n",w)
