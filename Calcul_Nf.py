# Créé par pjphi, le 22/03/2025 en Python 3.7
import scipy.optimize as resol

sigmaf = 3950
epsilonf = 1.5
b = -0.151
c = -0.761
E = 170000
SWT = 1

def f(x):
    return SWT - (((sigmaf**2)/E)*((2*x)**(2*b)) + sigmaf*epsilonf*((2*x)**(b+c)))

print(resol.brentq(f, 1, 100000000000))
