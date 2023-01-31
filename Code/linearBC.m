function g = linearBC(z,a,b,fa,fb)
g = fa*abs((b-z)/(b-a)) + fb*abs((z-a)/(b-a));