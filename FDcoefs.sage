###############################################################################
#
#	Finite Difference Coefficients derivation using SAGE http://sagemath.org
#	Copyright 2015 Nicolas Fressengeas 
#	http://fressengeas.net
#
###############################################################################
#
#
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
###############################################################################

import numpy as np

def FDcoeff(i):
    """
    Returns a valid variable name for a coefficient numbered i (a_i)
    """
    return(var(('a'+'_'+str(i)).replace('-','m') ))
    
def FDcoefflist((imin,imax)):
    """
    Returns a coefficient list for i between imin and imax (imax included).
    i=0 means the point where the discrete dervative is evaluated
    """
    return( list( (FDcoeff(i)) for i in range(imin,imax+1)) )
    
def FDequation(f,x,h,(imin,imax)):
    """
    f is a function
    x is a symbolic variable
    h is a symbolic space step
    imin and imax are the leftmost and rightmost points included in the Finite Difference
    FDequation(f,x,h,(imin,imax)) sums up the coeficients times the function f evaluated on each point
    """
    equation = 0;
    for i in range(imin,imax+1):
        equation+=f(x+i*h)*FDcoeff(i)
    return(equation)

def FDequationList(DiffOrder,x,h,(imin,imax)):
    """
    Returns a system of equation that is fullfilled by the coefficients assuming the discrete derivatives is exact for monomials from order 0 to imax-imin included
    """
    nbcoef=imax-imin+1
    PolyList=[(lambda i:(lambda x:x**i))(i) for i in range(nbcoef)]
    EqList=[FDequation(PolyList[i],x,h,(imin,imax))==diff(PolyList[i](x),x,DiffOrder) for i in range(nbcoef)]
    return(EqList)

def HighOrderFiniteDiff(DiffOrder,(imin,imax)):
    """Computes the coeeficients of a Finite Difference Approximation of the DiffOrder order derivative using the points in between imin and imax. The order of the approximation is given by imax-imin.
    Returns a list of couples (i,a) where i is the discretization point number refered to the point where the derivative is actually evaluated and a the corresponding coeficient in the Finite Difference Approximation. To get the actual coefficient, remeber that it should be divided by h^DiffOrder, where h is the discretization step.
    The method used is called Indeterminate Coefficient
    DiffOrder : differential order
    imin : with respect to the point where the derivative is computed, the leftmost neighbor (e.g. -1 for a centered 3 points first or second order derivative)
    imax : with respect to the point where the derivative is computed, the rightmost neighbor (e.g. +1 for a centered 3 points first or second order derivative)
    """
    if imax<imin :
        print("Rightmost point should be at the right of the leftmost one. And vice versa.")
        return([])
    if imax-imin<DiffOrder :
        print("Finite Difference of Nth order derivative require at least N+1 points")
        return([])
    var('x h')
    EqList=FDequationList(DiffOrder,x,h,(imin,imax))
    ResEq=solve(EqList,FDcoefflist((imin,imax)),solution_dict=True)
    ResCoefList=[coef.substitute(ResEq[0]) for coef in FDcoefflist((imin,imax))]
    ResICoef = [(i,ResCoefList[i-imin]*h^DiffOrder) for i in range(imin,imax+1)]
    return(ResICoef)

def DDiff(tab,h,DiffOrder,Accuracy):
    """
    DDiff(tab,h,DiffOrder,Accuracy) returns the discrete derivative of Order DiffOrder of the 1D array tab at the given Accuracy assuming a Finite Difference homogeneous h step. The size on the returned array is the same as that of the input one.
    It is en enhancement of the numpy.diff function.
    """
    tabN=tab.shape[0]
    npoints=Accuracy+DiffOrder+1  # because Accuracy==(imax-imin)-DiffOrder==npoints-1-DiffOrder
    #Far away from the borders
    # if npoints is odd, take a centered npoints scheme
    # if npoints is even, take a centered npoints+1 scheme
    maxshift=int(npoints/2)
    coefs=HighOrderFiniteDiff(DiffOrder,(-maxshift,+maxshift))
    Dtab=np.zeros(tabN-2*maxshift)
    for i in range(2*maxshift+1):
        Dtab+=coefs[i][1]*tab[i:tabN-2*maxshift+i]
    #let us now car about the borders
    #derivatives will be off centered while keeping the same number of points
    left=[]
    right=[]
    for i in range(maxshift):
        #left
        imin=-i
        imax=imin+npoints-1
        coefs=HighOrderFiniteDiff(DiffOrder,(imin,imax))
        deriv=0
        for j in range(npoints) :
            deriv+=coefs[j][1]*tab[j]
        left.append(deriv)
        Dleft=np.array(left)
        #right
        imax=+i
        imin=imax-npoints+1
        coefs=HighOrderFiniteDiff(DiffOrder,(imin,imax))
        deriv=0
        for j in range(npoints) :
            deriv+=coefs[j][1]*tab[tabN-npoints+j]
        right.insert(0,deriv)
        Dright=np.array(right)
    return(np.concatenate((Dleft,Dtab,Dright))/h^DiffOrder)
