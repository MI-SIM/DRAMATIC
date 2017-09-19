function [t,x,ts,xs,ximp]=ImpulseA(f,g,imp,x0,ti)

xtmp =@(x) imp(x(1))*x; % Rturns 0 if imp(x) if false
dx =@(t,x) imp(x(1))*f(t,xtmp(x));

t0 = ti(1);
Tend = ti(end);
Ttol = .01;
tolr = 1e-9;
tola = 1e-6;
options = rdpset('RelTol',tolr,'AbsTol',tola);
x0 = reshape(x0,1,[]);
ts = t0;
xs = x0;
ximp = x0;
t = t0;
x = x0;
while(t0<Tend-Ttol)
  [tnew,xnew]=RADAUsolver(dx,[t0,Tend],x0,options);
  x1new = xnew(:,1);
  ind1 = find([~imp(x1new);1],1)-1;
  t = [t; NaN; tnew(1:ind1)];
  x = [x; NaN*x0; xnew(1:ind1,:)];
  ts = [ts;t(end)];
  xs = [xs;g(x(end,:))];
  ximp = [ximp;x(end,:)];
  t0 = t(end);
  x0 = g(x(end,:));
end

end