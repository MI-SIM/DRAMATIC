function [tt,SS,CNH4C,CNO2C,CO2C,fafa,fnfn,fifi,x] = PartialNitritation_NO3(T,N)
%% Anammox and Partial Nitritation in a biofilm
% Using the model of Masic and Eberl (2015)
% Note that this is a function, not a script so that we can use 
% nested functions to simplify the use of parameters. 


%% Parameters
Ar = 0.17; d =5 ; rho = 10000;
V = 0.006; alpha = 1; E = 1000;
Ac = 0.0068; nc=50; A = Ar+Ac*nc;
kaO2 = 0.5; knO2 = 0.5;
muamax= 0.3082; kaNH4 = 0.1690;
munmax = 0.4015; knNO2 = 0.302;
ya = 0.15; yn = 0.041;
dNO2 = 7.93E-5; dNH4 = 9.13E-5; dO2 = 9.93E-5; 
SNO2in = 0; 
SNH4in = 30; SO2in = 5;
fxi = 0.1;  eta = 0.5;
ba = 0.04; bn = 0.08;


%% Discretization of domain; 
[D,x] = cheb(N);                                        % Chebyshev Differentiation Matrix
Dx= D; Dx(N+1,:) = zeros(N+1,1);                          % For the Neumann B.C. x(end) -> z=0
D2 = D^2; D2(N+1,:) = D(N+1,:); D2 = D2(2:N+1,2:N+1);   % Also for the Neumann B.C.
eps = 0.01; %dt = min([10e-6,N^(-4)/eps]);                 % timestep dependent on N
dt=0.001;
t=0;
% Chebyshev integration matrix:
% Sl is the indefinite integral from x to 1,
% int is the indefinite integral from -1 to x
[int,~] = intmat(N); 


%% Numerical Settings
tmax= T; points=9000;           % Maxtime and number of save points.
tplot = T/points;               % plot times
nplots = round(tmax/dt);        % # of iterations for timestepping
plotgap = round(tplot/dt);      % Gap between saving points for plotting
dt = tplot/plotgap;             % Timestep
tol = 1E-9;                     % Tolerance for the BVP

%% Exponential time differencing setup
M = 32; r = 15*exp(1i*pi*((1:M)-.5)/M); h=dt;
A1 = h*zeros(N+1); E1 = eye(N+1); E2 = eye(N+1);
I = eye(N+1); Z = zeros(N+1);
f1 = Z; f2 = Z; f3= Z; Q = Z;
for j = 1:M
    z = r(j);
    zIA = inv(z*I);
    Q = Q+h*zIA*(exp(z/2)-1);
    f1 = f1+h*zIA*(-4-z+exp(z)*(4-3*z+z^2))/z^2;
    f2 = f2+h*zIA*(2+z+exp(z)*(z-2))/z^2;
    f3 = f3+h*zIA*(-4-3*z-z*2+exp(z)*(4-z))/z^2;
end
f1 = real(f1/M); f2 = real(f2/M); f3 = real(f3/M); Q = real(Q/M);


%% Response Functions & Respiration Rates
muaO = @(y) y./(kaO2+y); % Oxygen consumption by AOB
munO = @(y) y./(knO2+y); % Oxygen consumption by NOB
mua = @(y,z) muamax*y./(kaNH4+y).*muaO(z); % Ammonium consumption by AOB
mun = @(y,z) munmax*y./(knNO2+y).*munO(z); % Nitrite consumption by NOB
Rxa = @(y,z) mua(y,z)-muamax*(1+eta)*ba*muaO(z); % 
Rxn = @(y,z) mun(y,z)-munmax*(1+eta)*bn*munO(z); % 
mui = @(xa,xn,y) (fxi+eta)*(ba*xa.*muaO(y)+bn*xn.*munO(y));
Rf = @(C,F) F(:,1).*Rxa(C(:,1),C(:,3))+F(:,2).*Rxn(C(:,2),C(:,3))+mui(F(:,1),F(:,2),C(:,3));
%% Biofilm Velocities and transport flux
function w = v(C,~,F) % biofilm velocity 
    w = int*(Rf(C,F));
end
function w = Lp(C,S,F) % L prime
        vel = v(C,S,F);
        w = 0.5*vel(1)*S(6)+alpha/(A*rho)*(S(3)+S(4)+S(5))-E*S(6)^2;
end
        

%% Biofilm Substrate Update Function
  function w = biofilmbvp(C,S,F)
    %CNH4 = C(:,1); CNO2 = C(:,2); CNO3 = C(:,3); CO2 = C(:,4);
    %SNH4 = S(1); SNO2=S(2); SNO3 = S(3); L = S(7);
    %fa = F(:,1); fn = F(:,2); fi = F(:,3); 
    change = 1;
      function w = RDE(C,S,F)
        factor= S(6)^2/4;
        RDENH4 = factor*mua(C(2:N+1,1),C(2:N+1,3)).*F(2:N+1,1)*rho/ya/dNH4;
        RDENO2 = factor*(mun(C(2:N+1,2),C(2:N+1,3)).*F(2:N+1,2)*rho/yn-mua(C(2:N+1,1),C(2:N+1,3))*rho/ya.*F(2:N+1,1))/dNO2;
        RDEO2 = factor/dO2*(((3.43-ya)/ya*mua(C(2:N+1,1),C(2:N+1,3))+(1+fxi)*ba*muaO(C(2:N+1,3))).*F(2:N+1,1)*rho+((1.14-yn)/yn*mun(C(2:N+1,2),C(2:N+1,3))+(1-fxi)*bn*munO(C(2:N+1,3))).*F(2:N+1,2)*rho);
        w = [RDENH4,RDENO2,RDEO2];
      end
    while change > tol
        CNew = RDE(C,S,F);
        CNH4new = [0;D2\CNew(:,1)]+S(1);
        CNO2new = [0;D2\CNew(:,2)]+S(2);
        CO2new = [0;D2\CNew(:,3)]+SO2;
        change = norm(C-[CNH4new,CNO2new,CO2new],inf);
        C = [CNH4new,CNO2new,CO2new];
    end
    w = C;
  end
%% Transport Equations for Biofilm
    function w = transport(C,S,F)
    % Transport equations
            vel = v(C,S,F);
            vp = Rf(C,F);
            fat = F(:,1).*Rxa(C(:,1),C(:,3))-F(:,1).*vp+(x-1)*(Lp(C,S,F)/S(6)).*Dx*F(:,1)-vel.*(Dx*F(:,1));
            fnt = F(:,2).*Rxn(C(:,2),C(:,3))-F(:,2).*vp+(x-1)*(Lp(C,S,F)/S(6)).*Dx*F(:,2)-vel.*(Dx*F(:,2));
            fit = mui(F(:,1),F(:,2),C(:,3))-F(:,3).*vp+(x-1)*(Lp(C,S,F)/S(6)).*Dx*F(:,3)-vel.*(Dx*F(:,3));
            w = [fat,fnt,fit];
        end
%% Planktonic Equations
    function w = planktonic(C,S,F)
    %CNH4 = C(:,1); CNO2 = C(:,2); CO2 = C(:,3); CNO3 = C(:,4);
    %SNH4 = S(1); SNO2=S(2); SNO3 = S(3);
    %Xa = S(4); Xn = S(5); Xi = S(6); L = S(7);
    %fa = F(:,1); fn = F(:,2); fi = F(:,3);
    % Biofilm Velocity and Flux Equations
    %vel = v(C,S,F);
    JNH4 = S(6)/2*rho*int*(mua(C(:,1),C(:,3)).*F(:,1)/ya)/dNH4; 
    JNO2 = S(6)/2*rho*int*(mun(C(:,2),C(:,3)).*F(:,2)/yn -mua(C(:,1),C(:,3)).*F(:,1)/ya)/dNO2;
    %JNO3 = S(6)/2*rho*int*(mun(C(:,2),C(:,3)).*F(:,2)/yn)/dNO3;
    % Main function;
    w = [d*(SNH4in-S(1))-1/V*(1/ya*mua(S(1),SO2)*S(3))-1/V*A*dNH4*JNH4(1);
         d*(SNO2in-S(2))-1/V*(1/yn*mun(S(2),SO2)*S(4)-1/ya*mua(S(1),SO2)*S(3))-1/V*A*dNO2*JNO2(1);
         S(3)*(mua(S(1),SO2)-d-alpha)+A*rho*F(1,1)*E*S(6)^2;
         S(4)*(mun(S(2),SO2)-d-alpha)+A*rho*F(1,2)*E*S(6)^2;
         (fxi+eta)*(S(3)*ba*muaO(SO2)+S(4)*bn*munO(SO2))-S(5)*(d+alpha)+A*rho*F(1,3)*E*S(6)^2;
         Lp(C,S,F)];
    end
%% Initial Conditions
% Planktonic
  SNO2 = SNO2in; SNH4 = SNH4in ;   % Substrates
  SO2 = SO2in ; %SNO3 = SNO3in;
  Xa = 20E-6;  Xn = 20E-6; Xi = 0;  % Biomass
  L = 1E-5;
  S = [SNH4;SNO2;Xa;Xn;Xi;L];
% Biofilm
  CO2 = SO2*ones(size(x));   CNO2 = SNO2*ones(size(x));
  CNH4 = SNH4*ones(size(x));%  CNO3 = SNO3*ones(size(x));
  C = [CNH4,CNO2,CO2];
  fa = 0.5*ones(size(x)); fn = 0.5*ones(size(x)); fi = zeros(size(x));
  F = [fa,fn,fi];
  C = biofilmbvp(C,S,F);
%% Set storage, preallocate size for speed
SS = zeros(6,points); tt = zeros(1,points);
fafa = zeros(N+1,points); fafa(:,1)=fa;
fnfn = zeros(N+1,points); fnfn(:,1)=fn;
fifi = zeros(N+1,points); fifi(:,1)=fi;
CNH4C = zeros(N+1,points); CNH4C(:,1) = C(:,1);
CNO2C = zeros(N+1,points); CNO2C(:,1) = C(:,2);
%CNO3C = zeros(N+1,points); CNO3C(:,1) = CNO3';
CO2C = zeros(N+1,points); CO2C(:,1) = C(:,3);
reverseStr = '';

%% time stepping function
%% Using the exponential timestepping method outlined in KassamTrefethen (2005);
    function [S,F] = timestep(C,S,F)
        Sk1 = planktonic(C,S,F);
        Nu = transport(C,S,F);
        a = F+Q*Nu;
        Sk2 = planktonic(C,S+0.5*dt*Sk1,a);
        Na = transport(C,S+0.5*dt*Sk1,a);
        b = F+Q*Na;
        Sk3= planktonic(C,S+0.5*dt*Sk2,b);
        Nb = transport(C,S+0.5*dt*Sk2,b);
        c = a+Q*(2*Nb-Nu);
        Sk4 = planktonic(C,S+dt*Sk3,c);
        Nc = transport(C,S+dt*Sk3,c);
        S = S+dt/6*(Sk1+2*Sk2+2*Sk3+Sk4);
        F = E1*F+f1*Nu+2*f2*(Na+Nb)+f3*Nc;
    end
        
        
%% Main Timestepping Loop
  for i = 1:nplots
    % Save into temporary Storage
    t= t+dt;
    Ct = biofilmbvp(C,S,F); 
    [S,F] = timestep(Ct,S,F);
    C = Ct;
    %% Save the data sometimes 
    if mod(i,plotgap)==0
      j = i/plotgap+1;
      tt(j) = t ; SS(:,j) = S;
      fafa(:,j) = F(:,1); fnfn(:,j) = F(:,2);
      fifi(:,j) = F(:,3);
      CNH4C(:,j) = C(:,1); CNO2C(:,j) = C(:,2);
      CO2C(:,j) = C(:,3);   %CNO3C = [CNO3C;C(:,4)']; 
      percentDone = 100 * i/nplots;
      msg = sprintf('Percent done: %3.3f', percentDone); %Don't forget this semicolon
      fprintf([reverseStr, msg]);
      reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
  end
fprintf([reverseStr, 'Percent done: 100', '\n']);
x = (1+x)/2;
end
