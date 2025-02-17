from sympy import *

numberOfDimensions = 3
numberOfVoigt = Matrix([[1,3,6]])

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
 

def PrintTensorTwo(TensorTwo,Name):
    print("")
    for iIdx in range(0,numberOfDimensions):
        for jIdx in range(0,numberOfDimensions):
            print(" %s[%1d,%1d] = %s" % (Name,iIdx+1,jIdx+1,TensorTwo[iIdx,jIdx]))
            
def PrintTensorTwoVoigt(TensorTwoV,Name):
    print("")
    for iIdx in range(0,numberOfVoigt[numberOfDimensions-1]):
        print(" %s[%1d] = %s" % (Name,iIdx+1,TensorTwoV[iIdx]))

def PrintTensorFourVoigt(TensorFourV,Name):
    print("")
    for iIdx in range(0,numberOfVoigt[numberOfDimensions-1]):
        for jIdx in range(0,numberOfVoigt[numberOfDimensions-1]):
            print(" %s[%1d,%1d] = %s" % (Name,iIdx+1,jIdx+1,TensorFourV[iIdx,jIdx]))

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

iV = Matrix([[  1, 0, 0, 0, 0, 0],
             [  0, 1, 0, 0, 0, 0],
             [  0, 0, 1, 0, 0, 0],
             [  0, 0, 0, 1, 0, 0],
             [  0, 0, 0, 0, 1, 0],
             [  0, 0, 0, 0, 0, 1]])

isharpV = Matrix([[  1, 0, 0,   0,   0, 0  ],
                  [  0, 1, 0,   0,   0,   0],
                  [  0, 0, 1,   0,   0,   0],
                  [  0, 0, 0, 1/2,   0,   0],
                  [  0, 0, 0,   0, 1/2,   0],
                  [  0, 0, 0,   0,   0, 1/2]])

kV = Matrix([[ 1/3, 1/3, 1/3, 0, 0, 0],
             [ 1/3, 1/3, 1/3, 0, 0, 0],
             [ 1/3, 1/3, 1/3, 0, 0, 0],
             [   0,   0,   0, 0, 0, 0],
             [   0,   0,   0, 0, 0, 0],
             [   0,   0,   0, 0, 0, 0]])

ksharpV = kV

mV = iV - kV

msharpV = isharpV - ksharpV

g = Matrix([[ 1, 0, 0 ],
            [ 0, 1, 0 ],
            [ 0, 0, 1 ]])

PrintTensorFourVoigt(mV,"mV")

mcbarV = mV*cbarV
ce1V = mcbarV*mV.T
PrintTensorFourVoigt( ce1V, "ce1V" )

trsigmabar = sigmabarV[0] + sigmabarV[1] + sigmabarV[2]

isharpmTV = isharpV*mV.T

twothirds = 2/3
J23 = J**(-twothirds)
ce2V = -twothirds*J23*trsigmabar*isharpmTV
PrintTensorFourVoigt( ce2V, "ce2V" )

sigmadevV = J23*mV*sigmabarV

sigmadev = VoigtToTensorTwoS( sigmadevV )

sigmadevxgV = TensorTwoCrossSToVoigt( sigmadev, g )
gxsigmadevV = TensorTwoCrossSToVoigt( g, sigmadev )

ce3V = -twothirds*(sigmadevxgV+gxsigmadevV)
PrintTensorFourVoigt( ce3V, "ce3V" )

cdevV = simplify(simplify(ce1V) + simplify(ce2V) + simplify(ce3V))
PrintTensorFourVoigt( cdevV, "cdevV" )
