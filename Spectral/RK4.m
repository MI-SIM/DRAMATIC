function w = RK4(f,x,h)

k1 = f(x);
k2 = f(x+0.5*h*k1);
k3 = f(x+0.5*h*k2);
k4 = f(x+h*k3);

w = x +h/6*(k1+2*k2+2*k3+k4);