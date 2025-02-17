from sympy import *

numberOfDimensions = 3

cbar = MatrixSymbol('cbar',6,6)
sigmbar = MatrixSymbol('sigmbar',1,6)
J = Symbol('J')

def VoigtToTensorTwoS(Voigt):
    return Matrix([[simplify(Voigt[0]),simplify(Voigt[5]),simplify(Voigt[4])],
                   [simplify(Voigt[5]),simplify(Voigt[1]),simplify(Voigt[3])],
                   [simplify(Voigt[4]),simplify(Voigt[3]),simplify(Voigt[2])]])

def TensorTwoSToVoigt(TensorTwo):
    return Matrix([[simplify(TensorTwo[0,0])],
                   [simplify(TensorTwo[1,1])],
                   [simplify(TensorTwo[2,2])],
                   [simplify(TensorTwo[1,2])],
                   [simplify(TensorTwo[0,2])],
                   [simplify(TensorTwo[0,1])]])

def TensorTwoCrossSToVoigt(TensorTwo1,TensorTwo2):
    return Matrix([[simplify(TensorTwo1[0,0]*TensorTwo2[0,0]),simplify(TensorTwo1[0,0]*TensorTwo2[1,1]),simplify(TensorTwo1[0,0]*TensorTwo2[2,2]),simplify(TensorTwo1[0,0]*TensorTwo2[1,2]),simplify(TensorTwo1[0,0]*TensorTwo2[0,2]),simplify(TensorTwo1[0,0]*TensorTwo2[0,1])],
                   [simplify(TensorTwo1[1,1]*TensorTwo2[0,0]),simplify(TensorTwo1[1,1]*TensorTwo2[1,1]),simplify(TensorTwo1[1,1]*TensorTwo2[2,2]),simplify(TensorTwo1[1,1]*TensorTwo2[1,2]),simplify(TensorTwo1[1,1]*TensorTwo2[0,2]),simplify(TensorTwo1[1,1]*TensorTwo2[0,1])],
                   [simplify(TensorTwo1[2,2]*TensorTwo2[0,0]),simplify(TensorTwo1[2,2]*TensorTwo2[1,1]),simplify(TensorTwo1[2,2]*TensorTwo2[2,2]),simplify(TensorTwo1[2,2]*TensorTwo2[1,2]),simplify(TensorTwo1[2,2]*TensorTwo2[0,2]),simplify(TensorTwo1[2,2]*TensorTwo2[0,1])],
                   [simplify(TensorTwo1[1,2]*TensorTwo2[0,0]),simplify(TensorTwo1[1,2]*TensorTwo2[1,1]),simplify(TensorTwo1[1,2]*TensorTwo2[2,2]),simplify(TensorTwo1[1,2]*TensorTwo2[1,2]),simplify(TensorTwo1[1,2]*TensorTwo2[0,2]),simplify(TensorTwo1[1,2]*TensorTwo2[0,1])],
                   [simplify(TensorTwo1[0,2]*TensorTwo2[0,0]),simplify(TensorTwo1[0,2]*TensorTwo2[1,1]),simplify(TensorTwo1[0,2]*TensorTwo2[2,2]),simplify(TensorTwo1[0,2]*TensorTwo2[1,2]),simplify(TensorTwo1[0,2]*TensorTwo2[0,2]),simplify(TensorTwo1[0,2]*TensorTwo2[0,1])],
                   [simplify(TensorTwo1[0,1]*TensorTwo2[0,0]),simplify(TensorTwo1[0,1]*TensorTwo2[1,1]),simplify(TensorTwo1[0,1]*TensorTwo2[2,2]),simplify(TensorTwo1[0,1]*TensorTwo2[1,2]),simplify(TensorTwo1[0,1]*TensorTwo2[0,2]),simplify(TensorTwo1[0,1]*TensorTwo2[0,1])]])
 
cbar1111 = Symbol('cbar1111')
cbar1122 = Symbol('cbar1122')
cbar1133 = Symbol('cbar1133')
cbar1123 = Symbol('cbar1123')
cbar1113 = Symbol('cbar1113')
cbar1112 = Symbol('cbar1112')
cbar2211 = Symbol('cbar2211')
cbar2222 = Symbol('cbar2222')
cbar2233 = Symbol('cbar2233')
cbar2223 = Symbol('cbar2223')
cbar2213 = Symbol('cbar2213')
cbar2212 = Symbol('cbar2212')
cbar3311 = Symbol('cbar3311')
cbar3322 = Symbol('cbar3322')
cbar3333 = Symbol('cbar3333')
cbar3323 = Symbol('cbar3323')
cbar3313 = Symbol('cbar3313')
cbar3312 = Symbol('cbar3312')
cbar2311 = Symbol('cbar2311')
cbar2322 = Symbol('cbar2322')
cbar2333 = Symbol('cbar2333')
cbar2323 = Symbol('cbar2323')
cbar2313 = Symbol('cbar2313')
cbar2312 = Symbol('cbar2312')
cbar1311 = Symbol('cbar1311')
cbar1322 = Symbol('cbar1322')
cbar1333 = Symbol('cbar1333')
cbar1323 = Symbol('cbar1323')
cbar1313 = Symbol('cbar1313')
cbar1312 = Symbol('cbar1312')
cbar1211 = Symbol('cbar1211')
cbar1222 = Symbol('cbar1222')
cbar1233 = Symbol('cbar1233')
cbar1223 = Symbol('cbar1223')
cbar1213 = Symbol('cbar1213')
cbar1212 = Symbol('cbar1212')

sigmabar11 = Symbol('sigmabar11')
sigmabar22 = Symbol('sigmabar22')
sigmabar33 = Symbol('sigmabar33')
sigmabar23 = Symbol('sigmabar23')
sigmabar13 = Symbol('sigmabar13')
sigmabar12 = Symbol('sigmabar12')

cbarV = Matrix([[ cbar1111, cbar1122, cbar1133, cbar1123, cbar1113, cbar1112 ],
                [ cbar2211, cbar2222, cbar2233, cbar2223, cbar2213, cbar2212 ],
                [ cbar3311, cbar3322, cbar3333, cbar3323, cbar3313, cbar3312 ],
                [ cbar2311, cbar2322, cbar2333, cbar2323, cbar2313, cbar2312 ],
                [ cbar1311, cbar1322, cbar1333, cbar1323, cbar1313, cbar1312 ],
                [ cbar1211, cbar1222, cbar1233, cbar1223, cbar1213, cbar1212 ]])

sigmabarV = Matrix([[sigmabar11],
                   [sigmabar22],
                   [sigmabar33],
                   [sigmabar23],
                   [sigmabar13],
                   [sigmabar12]])

mV = Matrix([[  2/3, -1/3, -1/3, 0, 0, 0],
            [ -1/3,  2/3, -1/3, 0, 0, 0],
            [ -1/3, -1/3,  2/3, 0, 0, 0],
            [    0,    0,    0, 1, 0, 0],
            [    0,    0,    0, 0, 1, 0],
            [    0,    0,    0, 0, 0, 1]])

isharpV = Matrix([[  1, 0, 0,   0,   0,   0],
                  [  0, 1, 0,   0,   0,   0],
                  [  0, 0, 1,   0,   0,   0],
                  [  0, 0, 0, 1/2,   0,   0],
                  [  0, 0, 0,   0, 1/2,   0],
                  [  0, 0, 0,   0,   0, 1/2]])

msharpV = Matrix([[  2/3, -1/3, -1/3, 0, 0, 0],
                  [ -1/3,  2/3, -1/3, 0, 0, 0],
                  [ -1/3, -1/3,  2/3, 0, 0, 0],
                  [    0,    0,    0, 1/2, 0, 0],
                  [    0,    0,    0, 0, 1/2, 0],
                  [    0,    0,    0, 0, 0, 1/2]])

g = Matrix([[ 1, 0, 0 ],
            [ 0, 1, 0 ],
            [ 0, 0, 1 ]])


mcbarV = mV*cbarV
mcbarmTV = mcbarV*mV.T

trsigmabar = sigmabarV[1] + sigmabarV[2] + sigmabarV[3]

msharpisharpV = msharpV*isharpV

sigmadevV = J**(-2/3)*mV*sigmabarV

sigmadev = VoigtToTensorTwoS( sigmadevV )

sigmadevxgV = TensorTwoCrossSToVoigt( sigmadev, g )
gxsigmadevV = TensorTwoCrossSToVoigt( g, sigmadev )

cdevV = mcbarmTV -2/3*J**(-2/3)*trsigmabar*msharpisharpV -2/3*(sigmadevxgV+gxsigmadevV)

print(" cdevV1111 = ",simplify(cdevV[0,0]))
print(" cdevV1122 = ",simplify(cdevV[0,1]))
print(" cdevV1133 = ",simplify(cdevV[0,2]))
print(" cdevV1123 = ",simplify(cdevV[0,3]))
print(" cdevV1113 = ",simplify(cdevV[0,4]))
print(" cdevV1112 = ",simplify(cdevV[0,5]))
print(" cdevV2211 = ",simplify(cdevV[1,0]))
print(" cdevV2222 = ",simplify(cdevV[1,1]))
print(" cdevV2233 = ",simplify(cdevV[1,2]))
print(" cdevV2223 = ",simplify(cdevV[1,3]))
print(" cdevV2213 = ",simplify(cdevV[1,4]))
print(" cdevV2212 = ",simplify(cdevV[1,5]))
print(" cdevV3311 = ",simplify(cdevV[2,0]))
print(" cdevV3322 = ",simplify(cdevV[2,1]))
print(" cdevV3333 = ",simplify(cdevV[2,2]))
print(" cdevV3323 = ",simplify(cdevV[2,3]))
print(" cdevV3313 = ",simplify(cdevV[2,4]))
print(" cdevV3312 = ",simplify(cdevV[2,5]))
print(" cdevV2311 = ",simplify(cdevV[3,0]))
print(" cdevV2322 = ",simplify(cdevV[3,1]))
print(" cdevV2333 = ",simplify(cdevV[3,2]))
print(" cdevV2323 = ",simplify(cdevV[3,3]))
print(" cdevV2313 = ",simplify(cdevV[3,4]))
print(" cdevV2312 = ",simplify(cdevV[3,5]))
print(" cdevV1311 = ",simplify(cdevV[4,0]))
print(" cdevV1322 = ",simplify(cdevV[4,1]))
print(" cdevV1333 = ",simplify(cdevV[4,2]))
print(" cdevV1323 = ",simplify(cdevV[4,3]))
print(" cdevV1313 = ",simplify(cdevV[4,4]))
print(" cdevV1312 = ",simplify(cdevV[4,5]))
print(" cdevV1211 = ",simplify(cdevV[5,0]))
print(" cdevV1222 = ",simplify(cdevV[5,1]))
print(" cdevV1233 = ",simplify(cdevV[5,2]))
print(" cdevV1223 = ",simplify(cdevV[5,3]))
print(" cdevV1213 = ",simplify(cdevV[5,4]))
print(" cdevV1212 = ",simplify(cdevV[5,5]))



