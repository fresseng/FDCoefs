###############################################################################
#
#	Finite Difference Coefficients derivation using SAGE http://sagemath.org
#	Usage Example Code
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
#
#	These examples show how to calculate the tables which are given
#	by Wikipedia : 
#	https://en.wikipedia.org/wiki/Finite_difference_coefficient
#
###############################################################################
#
#	Centered first order derivative up to 9 points
#
#
for j in range(1,5):
    show(HighOrderFiniteDiff(1,(-j,j)))
#
###############################################################################
#
#	Centered second order derivative up to 9 points
#
for j in range(1,5):
    show(HighOrderFiniteDiff(2,(-j,j)))
#
###############################################################################
#
#	Right first order derivative up to 6 points
#
for j in range(1,7):
    show(HighOrderFiniteDiff(1,(0,j)))
#
###############################################################################
#
#	Right second order derivative up to 8 points
#
for j in range(2,7):
    show(HighOrderFiniteDiff(2,(0,j)))
#
###############################################################################
#
#	And so on :
#		Change the first argument to get the derivative order wanted
#		Give your leftmost and rightmost points at a couple 
#			in the second argument
#		But remember :
#			N+1 points are needed to approximate the Nth order
#				derivative		
#			Leftmost point should be at the left of rightmost one
#			The accuracy order of the approximation is the number of
#				points minus the derivative order
#			In HighOrderFiniteDiff(N,(imin,imax)) :
#				Accuracy==(imax-imin)-N
#
###############################################################################
#
#
#	And now one example on how to use DDiff, which is an extension
#	of the numpy.diff discrete differential tool
#
# Create a 1D array of size 100
t=np.array([i**3 for i in range(100)])
# Assuming a space step of 0.1, compute the first order derivative
#  	at accuracy 3
dt=DDiff(t,0.1,1,3)
# Check that the size is unchanged
print(dt.shape)
# Check the derivative value
list_plot(dt)
#
###############################################################################
