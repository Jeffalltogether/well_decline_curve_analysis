# compute partial derivitaves for gradient descent
from sympy import Symbol, Derivative

qi,di,b,y,t = Symbol('qi'), Symbol('di'), Symbol('b'), Symbol('t'), Symbol('y')

SSE = 0.5*(y - qi*(1.0 - b*(di)*t)**(-1.0 / b))**2

dSSE_dqi = Derivative(SSE, qi).doit()
dSSE_ddi = Derivative(SSE, di).doit()
dSSE_db = Derivative(SSE, b).doit()

print dSSE_dqi
print dSSE_ddi
print dSSE_db