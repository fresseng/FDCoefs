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
