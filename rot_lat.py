from sympy import *


n0,n1,a=symbols("n0,n1,a",cls=Symbol,real=True)
m0,m1=symbols("m0,m1",cls=Symbol,real=True)
x=Rational(1,2)
y=sqrt(3)/2

R3=Matrix(
    [[-x,-y],[y,-x]]
)

v0=Matrix([a,0])
v1=Matrix([-x,y])*a

O=Matrix([0,0])
D=Matrix([0,sqrt(3)/3*a])
C=Matrix([a,0])
B=Matrix([x,y])*a
A=Matrix([-x,y])*a


tD=D+n0*v0+n1*v1

pprint(expand(tD[1]))
