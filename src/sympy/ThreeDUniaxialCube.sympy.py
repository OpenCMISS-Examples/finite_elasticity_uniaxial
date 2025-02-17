from sympy import *

numberOfDimensions = 3
numberOfVoigt = Matrix([[1,3,6]])
numberOfXi = 3
numberOfNodes = 8

dx, dy, dz = symbols('dx,dy,dz')
xi1, xi2, xi3  = symbols('xi1,xi2,xi3')
alpha, beta, gamma = symbols('alpha,beta,gamma')
F = MatrixSymbol('F',numberOfDimensions,numberOfDimensions)
Fbar = MatrixSymbol('Fbar',numberOfDimensions,numberOfDimensions)
C = MatrixSymbol('C',numberOfDimensions,numberOfDimensions)
Cbar = MatrixSymbol('Cbar',numberOfDimensions,numberOfDimensions)
B = MatrixSymbol('B',numberOfDimensions,numberOfDimensions)
Bbar = MatrixSymbol('Bbar',numberOfDimensions,numberOfDimensions)
b = MatrixSymbol('b',numberOfDimensions,numberOfDimensions)
J, I1, I2, I3 = symbols('J,I1,I2,I3')
I1bar, I2bar, I3bar = symbols('I1bar,I2bar,I3bar')
delI1delC = MatrixSymbol('delI1delC',numberOfDimensions,numberOfDimensions)
delI2delC = MatrixSymbol('delI2delC',numberOfDimensions,numberOfDimensions)
delI3delC = MatrixSymbol('delI3delC',numberOfDimensions,numberOfDimensions)
W, c1, c2, p = symbols('W,c1,c2,p')
delWdelI1, delWdelI2, delWdelI3 = symbols('delWdelI1,delWdelI2,delWdelI3')
I = MatrixSymbol('I',numberOfDimensions,numberOfDimensions)
S = MatrixSymbol('S',numberOfDimensions,numberOfDimensions)
Sbar = MatrixSymbol('Sbar',numberOfDimensions,numberOfDimensions)
Sdev = MatrixSymbol('Sdev',numberOfDimensions,numberOfDimensions)
Ssph = MatrixSymbol('Ssph',numberOfDimensions,numberOfDimensions)
sigma = MatrixSymbol('sigma',numberOfDimensions,numberOfDimensions)
cauchy = zeros(numberOfDimensions,numberOfDimensions)
dPhimdx = zeros(numberOfDimensions,1)
resid = zeros(numberOfDimensions*numberOfNodes+1,1)

I = Identity(numberOfDimensions)

symi=Matrix([[ 1, 0, 0, 0, 0, 0],
             [ 0, 1, 0, 0, 0, 0],
             [ 0, 0, 1, 0, 0, 0],
             [ 0, 0, 0, 1, 0, 0],
             [ 0, 0, 0, 0, 1, 0],
             [ 0, 0, 0, 0, 0, 1]])
print(" symi = ",symi)

symisharp=Matrix([[ 1, 0, 0, 0, 0, 0],
                  [ 0, 1, 0, 0, 0, 0],
                  [ 0, 0, 1, 0, 0, 0],
                  [ 0, 0, 0, 1/2, 0, 0],
                  [ 0, 0, 0, 0, 1/2, 0],
                  [ 0, 0, 0, 0, 0, 1/2]])
print(" symisharp = ",symisharp)

k=Matrix([[ 1/3, 1/3, 1/3, 0, 0, 0],
          [ 1/3, 1/3, 1/3, 0, 0, 0],
          [ 1/3, 1/3, 1/3, 0, 0, 0],
          [   0,   0,   0, 0, 0, 0],
          [   0,   0,   0, 0, 0, 0],
          [   0,   0,   0, 0, 0, 0]])
print(" k = ",k)

ksharp=Matrix([[ 1/3, 1/3, 1/3, 0, 0, 0],
               [ 1/3, 1/3, 1/3, 0, 0, 0],
               [ 1/3, 1/3, 1/3, 0, 0, 0],
               [   0,   0,   0, 0, 0, 0],
               [   0,   0,   0, 0, 0, 0],
               [   0,   0,   0, 0, 0, 0]])
print(" ksharp = ",ksharp)


m = symi - k
print(" m = ",m)

def PushForwardTensorTwo(TensorTwo):
    return (1/J)*(F*TensorTwo*F.T)

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

def TensorFourSToVoigt(TensorFour):
    return Matrix([[simplify(TensorFour[0,0,0,0]),simplify(TensorFour[0,0,1,1]),simplify(TensorFour[0,0,2,2]),simplify(TensorFour[0,0,1,2]),simplify(TensorFour[0,0,0,2]),simplify(TensorFour[0,0,0,1])],
                   [simplify(TensorFour[1,1,0,0]),simplify(TensorFour[1,1,1,1]),simplify(TensorFour[1,1,2,2]),simplify(TensorFour[1,1,1,2]),simplify(TensorFour[1,1,0,2]),simplify(TensorFour[1,1,0,1])],
                   [simplify(TensorFour[2,2,0,0]),simplify(TensorFour[2,2,1,1]),simplify(TensorFour[2,2,2,2]),simplify(TensorFour[2,2,1,2]),simplify(TensorFour[2,2,0,2]),simplify(TensorFour[2,2,0,1])],
                   [simplify(TensorFour[1,2,0,0]),simplify(TensorFour[1,2,1,1]),simplify(TensorFour[1,2,2,2]),simplify(TensorFour[1,2,1,2]),simplify(TensorFour[1,2,0,2]),simplify(TensorFour[1,2,0,1])],
                   [simplify(TensorFour[0,2,0,0]),simplify(TensorFour[0,2,1,1]),simplify(TensorFour[0,2,2,2]),simplify(TensorFour[0,2,1,2]),simplify(TensorFour[0,2,0,2]),simplify(TensorFour[0,2,0,1])],
                   [simplify(TensorFour[0,1,0,0]),simplify(TensorFour[0,1,1,1]),simplify(TensorFour[0,1,2,2]),simplify(TensorFour[0,1,1,2]),simplify(TensorFour[0,1,0,2]),simplify(TensorFour[1,1,0,1])]])

def TensorTwoCrossSToVoigt(TensorTwo1,TensorTwo2):
    return Matrix([[simplify(TensorTwo1[0,0]*TensorTwo2[0,0]),simplify(TensorTwo1[0,0]*TensorTwo2[1,1]),simplify(TensorTwo1[0,0]*TensorTwo2[2,2]),simplify(TensorTwo1[0,0]*TensorTwo2[1,2]),simplify(TensorTwo1[0,0]*TensorTwo2[0,2]),simplify(TensorTwo1[0,0]*TensorTwo2[0,1])],
                   [simplify(TensorTwo1[1,1]*TensorTwo2[0,0]),simplify(TensorTwo1[1,1]*TensorTwo2[1,1]),simplify(TensorTwo1[1,1]*TensorTwo2[2,2]),simplify(TensorTwo1[1,1]*TensorTwo2[1,2]),simplify(TensorTwo1[1,1]*TensorTwo2[0,2]),simplify(TensorTwo1[1,1]*TensorTwo2[0,1])],
                   [simplify(TensorTwo1[2,2]*TensorTwo2[0,0]),simplify(TensorTwo1[2,2]*TensorTwo2[1,1]),simplify(TensorTwo1[2,2]*TensorTwo2[2,2]),simplify(TensorTwo1[2,2]*TensorTwo2[1,2]),simplify(TensorTwo1[2,2]*TensorTwo2[0,2]),simplify(TensorTwo1[2,2]*TensorTwo2[0,1])],
                   [simplify(TensorTwo1[1,2]*TensorTwo2[0,0]),simplify(TensorTwo1[1,2]*TensorTwo2[1,1]),simplify(TensorTwo1[1,2]*TensorTwo2[2,2]),simplify(TensorTwo1[1,2]*TensorTwo2[1,2]),simplify(TensorTwo1[1,2]*TensorTwo2[0,2]),simplify(TensorTwo1[1,2]*TensorTwo2[0,1])],
                   [simplify(TensorTwo1[0,2]*TensorTwo2[0,0]),simplify(TensorTwo1[0,2]*TensorTwo2[1,1]),simplify(TensorTwo1[0,2]*TensorTwo2[2,2]),simplify(TensorTwo1[0,2]*TensorTwo2[1,2]),simplify(TensorTwo1[0,2]*TensorTwo2[0,2]),simplify(TensorTwo1[0,2]*TensorTwo2[0,1])],
                   [simplify(TensorTwo1[0,1]*TensorTwo2[0,0]),simplify(TensorTwo1[0,1]*TensorTwo2[1,1]),simplify(TensorTwo1[0,1]*TensorTwo2[2,2]),simplify(TensorTwo1[0,1]*TensorTwo2[1,2]),simplify(TensorTwo1[0,1]*TensorTwo2[0,2]),simplify(TensorTwo1[0,1]*TensorTwo2[0,1])]])
 
def BarToDeviatoricVoigt(TensorBarV):
    return J**(-2/3)*m*TensorBarV

def SphericalVoigt(TensorV):
    return J**(-2/3)*m*TensorV

def Evaluate(Object):
    return Object.evalf(subs={dx:1,dy:1,dz:1,alpha:0.1,beta:0.1,gamma:0.1,c1:0.5,c2:0.1,p:0.1})

def PrintTensorTwo(TensorTwo,Name):
    print("")
    for iIdx in range(0,numberOfDimensions):
        for jIdx in range(0,numberOfDimensions):
            print(" %s[%1d,%1d] = %s" % (Name,iIdx+1,jIdx+1,simplify(TensorTwo[iIdx,jIdx])))
            
def PrintTensorTwoVoigt(TensorTwoV,Name):
    print("")
    for iIdx in range(0,numberOfVoigt[numberOfDimensions-1]):
        print(" %s[%1d] = %s" % (Name,iIdx+1,simplify(TensorTwoV[iIdx])))

def PrintTensorFourVoigt(TensorFourV,Name):
    print("")
    for iIdx in range(0,numberOfVoigt[numberOfDimensions-1]):
        for jIdx in range(0,numberOfVoigt[numberOfDimensions-1]):
            print(" %s[%1d,%1d] = %s" % (Name,iIdx+1,jIdx+1,simplify(TensorFourV[iIdx,jIdx])))

W = c1*(I1-numberOfDimensions) + c2*(I2-numberOfDimensions)
print(" W = ",W)
Wbar = c1*(I1bar-numberOfDimensions) + c2*(I2bar-numberOfDimensions)
print(" Wbar = ",Wbar)

delWdelI1 = diff(W,I1)
delWdelI2 = diff(W,I2)
delWdelI3 = diff(W,I3)
print(" delWdelI1 = ",delWdelI1)
print(" delWdelI2 = ",delWdelI2)
print(" delWdelI3 = ",delWdelI3)

delWbardelI1bar = diff(Wbar,I1bar)
delWbardelI2bar = diff(Wbar,I2bar)
delWbardelI3bar = diff(Wbar,I3bar)
print(" delWbardelI1bar = ",delWbardelI1bar)
print(" delWbardelI2bar = ",delWbardelI2bar)
print(" delWbardelI3bar = ",delWbardelI3bar)

F = Matrix([[ 1.0+alpha,      0.0,       0.0 ],
            [       0.0, 1.0-beta,       0.0 ],
            [       0.0,      0.0, 1.0-gamma ]])
PrintTensorTwo(F,"F")

J = F.det()
print(" J = ",J)

Fbar = J**(-1/3)*F
PrintTensorTwo(Fbar,"Fbar")

C = F.T*F
PrintTensorTwo(C,"C")

b = F*F.T

B = C.inv()

Cbar = J**(-2/3)*C
eCbar = Evaluate(Cbar)
PrintTensorTwo(eCbar,"eCbar")

Bbar = Cbar.inv()

I1 = simplify(Trace(C))
I2 = simplify(Trace(B)*C.det())
I3 = simplify(C.det())
print(" I1 = ",I1)
print(" I2 = ",I2)
print(" I3 = ",I3)

I1bar = simplify(Trace(Cbar))
I2bar = simplify(Trace(Bbar)*Cbar.det())
I3bar = simplify(Cbar.det())
print(" I1bar = ",I1bar)
print(" I2bar = ",I2bar)
print(" I3bar = ",I3bar)

delI1delC = I
delI2delC = I1*I - C
delI3delC = I3*C

delI1bardelCbar = I
delI2bardelCbar = I1bar*I - Cbar
delI3bardelCbar = I3bar*Cbar

Sbar = simplify(2.0*(delWbardelI1bar*delI1bardelCbar+delWbardelI2bar*delI2bardelCbar))
PrintTensorTwo(Sbar,"Sbar")

SbarV=TensorTwoSToVoigt(Sbar)
PrintTensorTwoVoigt(SbarV,"SbarV")

eSbarV = Evaluate(SbarV)
PrintTensorTwoVoigt(eSbarV,"eSbarV")

SdevV=BarToDeviatoricVoigt(SbarV)
PrintTensorTwoVoigt(SdevV,"SdevV")

eSdevV = Evaluate(SdevV)
PrintTensorTwoVoigt(eSdevV,"eSdevV")

Ssph = simplify(p*J*B)
PrintTensorTwo(Ssph,"Ssph")

SsphV = TensorTwoSToVoigt(Ssph)
eSsphV = Evaluate(SsphV)
PrintTensorTwoVoigt(eSsphV,"eSsphV")

SV = SdevV + SsphV
eSV = Evaluate(SV)
PrintTensorTwoVoigt(eSV,"eSV")

S = VoigtToTensorTwoS(SV)
PrintTensorTwo(S,"S")

sigmadevbar = PushForwardTensorTwo(Sbar)
PrintTensorTwo(sigmadevbar,"sigmadevbar")

sigmadevbarV = TensorTwoSToVoigt(sigmadevbar)
PrintTensorTwoVoigt(sigmadevbarV,"sigmadevbarV")

esigmadevbarV = Evaluate(sigmadevbarV)
PrintTensorTwoVoigt(esigmadevbarV,"esigmadevbarV")

sigmadevV=BarToDeviatoricVoigt(sigmadevbarV)

sigmadev=VoigtToTensorTwoS(sigmadevV)

PrintTensorTwo(sigmadev,"sigmadev")

PrintTensorTwoVoigt(sigmadevV,"sigmadevV")

esigmadev = Evaluate(sigmadev)
PrintTensorTwo(esigmadev,"esigmadev")

sigmasph = -p*I
PrintTensorTwo(sigmasph,"sigmasph")

sigmasphV = TensorTwoSToVoigt(sigmasph)
PrintTensorTwoVoigt(sigmasphV,"sigmasphV")

esigmasphV = Evaluate(sigmasphV)
PrintTensorTwoVoigt(esigmasphV,"esigmasphV")

sigmaV = sigmadevV+sigmasphV
PrintTensorTwoVoigt(sigmaV,"sigmaV")

cauchy = VoigtToTensorTwoS(sigmaV)
PrintTensorTwo(cauchy,"cauchy")

ecauchy = Evaluate(cauchy)
PrintTensorTwo(ecauchy,"ecauchy")

S = simplify(2.0*(delWdelI1*delI1delC+delWdelI2*delI2delC)-p*J*B)
sigmaS = PushForwardTensorTwo(S)
PrintTensorTwo(sigmaS,"sigmaS")
 
dPhidxi = Matrix([[ -(1-xi2)*(1-xi3), -(1-xi1)*(1-xi3), -(1-xi1)*(1-xi2) ],
                  [  (1-xi2)*(1-xi3),     -xi1*(1-xi3),     -xi1*(1-xi2) ],
                  [     -xi2*(1-xi3),  (1-xi1)*(1-xi3),     -(1-xi1)*xi2 ],
                  [      xi2*(1-xi3),       xi1*(1-xi3),        -xi1*xi2 ],
                  [     -(1-xi2)*xi3,     -(1-xi1)*xi3,  (1-xi1)*(1-xi2) ],
                  [      (1-xi2)*xi3,         -xi1*xi3,      xi1*(1-xi2) ],
                  [         -xi2*xi3,      (1-xi1)*xi3,      (1-xi1)*xi2 ],
                  [          xi2*xi3,          xi1*xi3,          xi1*xi2 ]])
gu = zeros(numberOfDimensions,numberOfDimensions)

JB = dx*dy*dz

dxidx = Matrix([[ 1/dx,    0,    0 ],
                [    0, 1/dy,    0 ],
                [    0,    0, 1/dz ]])       
 
for xiIdx1 in range(0,numberOfXi):
    for xiIdx2 in range(0,numberOfXi):
        gu[xiIdx1,xiIdx2]=0
        for coordIdx in range(0,numberOfDimensions):
            gu[xiIdx1,xiIdx2] = gu[xiIdx1,xiIdx2] + dxidx[xiIdx1,coordIdx]*dxidx[xiIdx2,coordIdx]
#print("gu = ",gu)

index = 0
for jIdx in range(0,numberOfDimensions):
    for mIdx in range(0,numberOfNodes):
        integrand = 0
        for iIdx in range(0,numberOfDimensions):
            dPhimdx[iIdx]=0
            for xiIdx in range(0,numberOfXi):
                dPhimdx[iIdx] = dPhimdx[iIdx] + dPhidxi[mIdx,xiIdx]*dxidx[xiIdx,iIdx]
            print(" dPhimdx[%1d,%2d,%1d] = %s" % (iIdx+1,mIdx+1,jIdx+1,dPhimdx[iIdx]))
            integrand = integrand + cauchy[iIdx,jIdx]*dPhimdx[iIdx]
        integrand = integrand*JB
        integral1 = integrate(integrand,(xi3,0,1))
        integral2 = integrate(integral1,(xi2,0,1))
        integral3 = simplify(integrate(integral2,(xi1,0,1)))
        resid[index] = simplify(integral3)
        #print("index = %2d, j=%2d, m=%2d, K = %s" % (index,jIdx+1,mIdx+1,integral3))
        index = index + 1
#Find pressure residual        
integrand = -(J-1)*JB
integral1 = integrate(integrand,(xi3,0,1))
integral2 = integrate(integral1,(xi2,0,1))
integral3 = simplify(integrate(integral2,(xi1,0,1)))
resid[index] = simplify(integral3)

for idx in range(0,numberOfDimensions*numberOfNodes+1):
    sresid = simplify(resid[idx])
    print(" resid[%2d] = %s" % (idx+1,sresid))
    
eresid = Evaluate(resid)
for idx in range(0,numberOfDimensions*numberOfNodes+1):
    print(" eresid[%2d] = %s" % (idx+1,eresid[idx]))

CbarV = Matrix([[    0, 4*c2, 4*c2,     0,     0,     0 ],
                [ 4*c2,    0, 4*c2,     0,     0,     0 ],
                [ 4*c2, 4*c2,    0,     0,     0,     0 ],
                [    0,    0,    0, -4*c2,     0,     0 ],
                [    0,    0,    0,     0, -4*c2,     0 ],
                [    0,    0,    0,     0,     0, -4*c2 ]])
    
FxF = TensorTwoCrossSToVoigt(F,F)

cbarV = (1/J)*(FxF*CbarV)*FxF.T               
PrintTensorFourVoigt(cbarV,"cbarV")

cdev1V = (m*cbarV)*m.T
PrintTensorFourVoigt(cdev1V,"cdev1V")

trsigmabar = sigmadevbarV[0]+sigmadevbarV[1]+sigmadevbarV[2]
cdev2V = 2/3*J**(-2/3)*trsigmabar*symisharp*m.T
PrintTensorFourVoigt(cdev2V,"cdev2V")
    
IxsigmadevV = TensorTwoCrossSToVoigt(I,sigmadev)
PrintTensorFourVoigt(IxsigmadevV,"IxsigmadevV")
sigmadevxIV = TensorTwoCrossSToVoigt(sigmadev,I)
PrintTensorFourVoigt(sigmadevxIV,"sigmadevxIV")
cdev3V = 2/3*(IxsigmadevV + sigmadevxIV)
PrintTensorFourVoigt(cdev3V,"cdev3V")



cdevV = cdev1V + cdev2V + cdev3V
PrintTensorFourVoigt(cdevV,"cdevV")

ecdevV = Evaluate(cdevV)
PrintTensorFourVoigt(ecdevV,"ecdevV")

csph1V = -3*J*p*ksharp
csph2V = 2*J*p*symisharp

csphV = csph1V + csph2V
PrintTensorFourVoigt(csphV,"csphV")

ecsphV = Evaluate(csphV)
PrintTensorFourVoigt(ecsphV,"ecsphV")

cV = cdevV + csphV

ecV = Evaluate(cV)
PrintTensorFourVoigt(ecV,"ecV")
