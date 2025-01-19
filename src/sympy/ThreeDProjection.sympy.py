from sympy import *

numberOfDimensions = 3

c = MatrixSymbol('c',6,6)

c1111 = Symbol('c1111')
c1122 = Symbol('c1122')
c1133 = Symbol('c1133')
c1123 = Symbol('c1123')
c1113 = Symbol('c1113')
c1112 = Symbol('c1112')
c2211 = Symbol('c2211')
c2222 = Symbol('c2222')
c2233 = Symbol('c2233')
c2223 = Symbol('c2223')
c2213 = Symbol('c2213')
c2212 = Symbol('c2212')
c3311 = Symbol('c3311')
c3322 = Symbol('c3322')
c3333 = Symbol('c3333')
c3323 = Symbol('c3323')
c3313 = Symbol('c3313')
c3312 = Symbol('c3312')
c2311 = Symbol('c2311')
c2322 = Symbol('c2322')
c2333 = Symbol('c2333')
c2323 = Symbol('c2323')
c2313 = Symbol('c2313')
c2312 = Symbol('c2312')
c1311 = Symbol('c1311')
c1322 = Symbol('c1322')
c1333 = Symbol('c1333')
c1323 = Symbol('c1323')
c1313 = Symbol('c1313')
c1312 = Symbol('c1312')
c1211 = Symbol('c1211')
c1222 = Symbol('c1222')
c1233 = Symbol('c1233')
c1223 = Symbol('c1223')
c1213 = Symbol('c1213')
c1212 = Symbol('c1212')

c = Matrix([[ c1111, c1122, c1133, c1123, c1113, c1112 ],
            [ c2211, c2222, c2233, c2223, c2213, c2212 ],
            [ c3311, c3322, c3333, c3323, c3313, c3312 ],
            [ c2311, c2322, c2333, c2323, c2313, c2312 ],
            [ c1311, c1322, c1333, c1323, c1313, c1312 ],
            [ c1211, c1222, c1233, c1223, c1213, c1212 ]])


t11 = Symbol('t11')
t22 = Symbol('t22')
t33 = Symbol('t33')
t23 = Symbol('t23')
t13 = Symbol('t13')
t12 = Symbol('t12')

t = Matrix([[ t11, t22, t33, t23, t13, t12 ]])

i = Matrix([[ 1, 1, 1, 0, 0, 0]])

ti = t.T*i
print(" ti = ",ti)
it = i.T*t
print(" it = ",it)

tiit = ti+it
print(" tiit = ",tiit)

P = Matrix([[  2/3, -1/3, -1/3, 0, 0, 0],
            [ -1/3,  2/3, -1/3, 0, 0, 0],
            [ -1/3, -1/3,  2/3, 0, 0, 0],
            [    0,    0,    0, 1, 0, 0],
            [    0,    0,    0, 0, 1, 0],
            [    0,    0,    0, 0, 0, 1]])

Pc = P*c
PcP = Pc*P

print(" PcP = ", PcP)

print(" PcP1111 = ",simplify(PcP[0,0]))
print(" PcP1122 = ",simplify(PcP[0,1]))
print(" PcP1133 = ",simplify(PcP[0,2]))
print(" PcP1123 = ",simplify(PcP[0,3]))
print(" PcP1113 = ",simplify(PcP[0,4]))
print(" PcP1112 = ",simplify(PcP[0,5]))
print(" PcP2211 = ",simplify(PcP[1,0]))
print(" PcP2222 = ",simplify(PcP[1,1]))
print(" PcP2233 = ",simplify(PcP[1,2]))
print(" PcP2223 = ",simplify(PcP[1,3]))
print(" PcP2213 = ",simplify(PcP[1,4]))
print(" PcP2212 = ",simplify(PcP[1,5]))
print(" PcP3311 = ",simplify(PcP[2,0]))
print(" PcP3322 = ",simplify(PcP[2,1]))
print(" PcP3333 = ",simplify(PcP[2,2]))
print(" PcP3323 = ",simplify(PcP[2,3]))
print(" PcP3313 = ",simplify(PcP[2,4]))
print(" PcP3312 = ",simplify(PcP[2,5]))
print(" PcP2311 = ",simplify(PcP[3,0]))
print(" PcP2322 = ",simplify(PcP[3,1]))
print(" PcP2333 = ",simplify(PcP[3,2]))
print(" PcP2323 = ",simplify(PcP[3,3]))
print(" PcP2313 = ",simplify(PcP[3,4]))
print(" PcP2312 = ",simplify(PcP[3,5]))
print(" PcP1311 = ",simplify(PcP[4,0]))
print(" PcP1322 = ",simplify(PcP[4,1]))
print(" PcP1333 = ",simplify(PcP[4,2]))
print(" PcP1323 = ",simplify(PcP[4,3]))
print(" PcP1313 = ",simplify(PcP[4,4]))
print(" PcP1312 = ",simplify(PcP[4,5]))
print(" PcP1211 = ",simplify(PcP[5,0]))
print(" PcP1222 = ",simplify(PcP[5,1]))
print(" PcP1233 = ",simplify(PcP[5,2]))
print(" PcP1223 = ",simplify(PcP[5,3]))
print(" PcP1213 = ",simplify(PcP[5,4]))
print(" PcP1212 = ",simplify(PcP[5,5]))



trt = t11 + t22 + t33

trtP = trt*P

ce2 = 2/3*trtP - 2/3*tiit

print(" ce2 = ",ce2)

c = PcP + ce2

print(" c = ", c)

print(" c1111 = ",simplify(c[0,0]))
print(" c1122 = ",simplify(c[0,1]))
print(" c1133 = ",simplify(c[0,2]))
print(" c1123 = ",simplify(c[0,3]))
print(" c1113 = ",simplify(c[0,4]))
print(" c1112 = ",simplify(c[0,5]))

print(" c2211 = ",simplify(c[1,0]))
print(" c2222 = ",simplify(c[1,1]))
print(" c2233 = ",simplify(c[1,2]))
print(" c2223 = ",simplify(c[1,3]))
print(" c2213 = ",simplify(c[1,4]))
print(" c2212 = ",simplify(c[1,5]))

print(" c3311 = ",simplify(c[2,0]))
print(" c3322 = ",simplify(c[2,1]))
print(" c3333 = ",simplify(c[2,2]))
print(" c3323 = ",simplify(c[2,3]))
print(" c3313 = ",simplify(c[2,4]))
print(" c3312 = ",simplify(c[2,5]))

print(" c2311 = ",simplify(c[3,0]))
print(" c2322 = ",simplify(c[3,1]))
print(" c2333 = ",simplify(c[3,2]))
print(" c2323 = ",simplify(c[3,3]))
print(" c2313 = ",simplify(c[3,4]))
print(" c2312 = ",simplify(c[3,5]))

print(" c1311 = ",simplify(c[4,0]))
print(" c1322 = ",simplify(c[4,1]))
print(" c1333 = ",simplify(c[4,2]))
print(" c1323 = ",simplify(c[4,3]))
print(" c1313 = ",simplify(c[4,4]))
print(" c1312 = ",simplify(c[4,5]))

print(" c1211 = ",simplify(c[5,0]))
print(" c1222 = ",simplify(c[5,1]))
print(" c1233 = ",simplify(c[5,2]))
print(" c1223 = ",simplify(c[5,3]))
print(" c1213 = ",simplify(c[5,4]))
print(" c1212 = ",simplify(c[5,5]))


