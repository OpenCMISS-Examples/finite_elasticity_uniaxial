from sympy import *

numberOfDimensions = 2

c = MatrixSymbol('c',3,3)

c1111 = Symbol('c1111')
c1122 = Symbol('c1122')
c1112 = Symbol('c1112')
c2211 = Symbol('c2211')
c2222 = Symbol('c2222')
c2212 = Symbol('c2212')
c1211 = Symbol('c1211')
c1222 = Symbol('c1222')
c1212 = Symbol('c1212')

c = Matrix([[ c1111, c1122, c1112 ],
            [ c2211, c2222, c2212 ],
            [ c1211, c1222, c1212 ]])

t11 = Symbol('t11')
t22 = Symbol('t22')
t12 = Symbol('t12')

t = Matrix([[ t11, t22, t12 ]])

i = Matrix([[ 1, 1, 0 ]])

ti = t.T*i
print(" ti = ",ti)
it = i.T*t
print(" it = ",it)

tiit = ti+it
print(" tiit = ",tiit)

P = Matrix([[  1/2, -1/2, 0 ],
            [ -1/2,  1/2, 0 ],
            [    0,    0, 1 ]])

Pc = P*c
PcP = Pc*P

print(" PcP = ", PcP)

print(" PcP1111 = ",simplify(PcP[0,0]))
print(" PcP1122 = ",simplify(PcP[0,1]))
print(" PcP1112 = ",simplify(PcP[0,2]))
print(" PcP2211 = ",simplify(PcP[1,0]))
print(" PcP2222 = ",simplify(PcP[1,1]))
print(" PcP2212 = ",simplify(PcP[1,2]))
print(" PcP1211 = ",simplify(PcP[2,0]))
print(" PcP1222 = ",simplify(PcP[2,1]))
print(" PcP1212 = ",simplify(PcP[2,2]))

trt = t11 + t22

trtP = trt*P

ce2 = 1/2*trtP - 1/2*tiit

print(" ce2 = ",ce2)

c = PcP + ce2

print(" c = ", c)


print(" c1111 = ",simplify(c[0,0]))
print(" c1122 = ",simplify(c[0,1]))
print(" c1112 = ",simplify(c[0,2]))
print(" c2211 = ",simplify(c[1,0]))
print(" c2222 = ",simplify(c[1,1]))
print(" c2212 = ",simplify(c[1,2]))
print(" c1211 = ",simplify(c[2,0]))
print(" c1222 = ",simplify(c[2,1]))
print(" c1212 = ",simplify(c[2,2]))
