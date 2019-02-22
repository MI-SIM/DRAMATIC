%% competition.m - Microbial competition in a biofilm
function [tt,x,SS,CC,f1f1,f2f2] = competition(T,N)

%% Parameters
Ar = 0.17; d =5 ; rho = 10000;
V = 0.006; alpha = 1; E = 1000;
Ac = 0.0068; nc=50; A = Ar+Ac*nc;
mu1max= 0.3082; k1 = 0.1690;
mu2max = 0.4015; k2 = 0.302;
y1 = 0.15; y2 = 0.041;
diffuse_C = 7.93E-5;
Sin = 30;
%% Initial conditions & Storage
S0 = Sin; x1_0= 2E-5; x2_0 = 2E-5;
L = 1E-5; 
S = [S0;x1_0;x2_0;L];
C = S0*ones(N+1,1);
f1_0 = 0.5*ones(N+1,1);
f2_0=0.5*ones(N+1,1);
F = [f1_0,f2_0];
%% Spatial Grid & Differentiation Matrix
% Rather than computing everything on the interval [0,L(t)], we change
% variables to compute solutions on the grid [-1,1] and make use of
% spectral methods to do most of our heavy lifting. See Trefethen (2001).
[D,x] = cheb(N);                  % Set up chebyshev differentiation Matrix
[int,~] = intmat(N);              % Integration matrix
Dx = D; Dx(N+1,:) = zeros(N+1,1);   % Neumann B.C. at x=-1 (x(1))
D2 = D^2; D2(N+1,:) = D(N+1,:);       % Neumann B.C. at x=-1 (x(1))
D2 = D2(2:N+1,2:N+1);                 % Dirichlet B.C at x=1 (x(end)) 

%% Time stepping Settings
dt = 0.005; t = 0;                 % Timestep and Initial Time
points = 900;                     % Number of points to save
tplot = T/points;                 % time at which to plot
nplots = round(T/dt);             % Number of total time points 
plotgap = round(tplot/dt);        % Gap between saved points
dt = tplot/plotgap;               % Adjusted timestep

%% Exponential time differencing setup
% Exponential time differencing is a stiff method for PDEs
% and can be used effectively for ODES coupled to ODEs. This version is
% similar to the standard RK4 Solver. Here we precompute various matrix and
% scalar quantities to save time. See Kassam Trefethen (2004)
% For PDEs:
    M = 32; r = 15*exp(1i*pi*((1:M)-.5)/M);
    I = eye(N+1); Z = zeros(N+1);
    fU1 = Z; fU2 = Z; fU3= Z; QU = Z;
    for jk = 1:M
        z = r(jk);
        zIA = inv(z*I);
        QU = QU+dt*zIA*(exp(z/2)-1);
        fU1 = fU1+dt*zIA*(-4-z+exp(z)*(4-3*z+z^2))/z^2;
        fU2 = fU2+dt*zIA*(2+z+exp(z)*(z-2))/z^2;
        fU3 = fU3+dt*zIA*(-4-3*z-z*2+exp(z)*(4-z))/z^2;
    end
    fU1 = real(fU1/M); fU2 = real(fU2/M); fU3 = real(fU3/M); QU = real(QU/M);
% For ODEs:
    FS1 = 0; FS2=0; FS3=0; QS=0;
    for jk = 1:M
        z = r(jk);
        zIA = 1/z;
        QS = QS +dt*zIA*(exp(z/2)-1);
        FS1 = FS1+dt*zIA*(-4-z+exp(z)*(4-3*z+z^2))/z^2;
        FS2 = FS2+dt*zIA*(2+z+exp(z)*(z-2))/z^2;
        FS3 = FS3+dt*zIA*(-4-3*z-z*2+exp(z)*(4-z))/z^2;
    end
    FS1 = real(FS1/M); FS2 = real(FS2/M); FS3 = real(FS3/M); QS=real(QS/M);
    
%% Uptake and growth functions
mu1 = @(y) mu1max*y./(k1+y.^2);
mu2 = @(y) mu2max*y./(k2+y.^2);
function w = vel(C,S,F)
    w = S(4)/2*int*(F(:,1).*mu1(C)+F(:,2).*mu2(C));
end
function w = Lp(C,S,F)
    v = vel(C,S,F);
    w = v(1)+alpha*S(2)/rho-E*S(4)^2;
end


%% Biofilm Community Equations
function w = transport(C,S,F)
    dx = (x+1)/S(4)*Lp(C,S,F);
    v = 2/S(4)*vel(C,S,F);
    DF1 = mu1(C).*F(:,1)-F(:,1).*(mu1(C).*F(:,1)+mu2(C).*F(:,2))-dx.*Dx*F(:,1)-v.*Dx*F(:,1);
    DF2 = mu2(C).*F(:,2)-F(:,2).*(mu1(C).*F(:,1)+mu2(C).*F(:,2))-dx.*Dx*F(:,2)-v.*Dx*F(:,2);
    w = [DF1,DF2];
end
%% Planktonic Equations
function DS = planktonic(C,S,F)
    % Flux
    J =rho*S(4)/2*int*(1/y1*mu1(C).*F(:,1)+1/y2*mu2(C).*F(:,2));
    DS = [d*(Sin-S(1))-1/y1*mu1(S(1))*S(2)-1/y2*mu2(S(1))*S(3)-J(1);
        -d*S(2)+mu1(S(1))*S(2)-alpha*S(2)+rho*E*S(4)^2*F(1,1);
        -d*S(3)+mu2(S(1))*S(3)-alpha*S(3)+rho*E*S(4)^2*F(1,2);
        Lp(C,S,F)];
end

%% Biofilm Concentration Equations
function w = BVP(C,S,F)
    DC = S(4)^2/4*(mu1(C(1:N)).*F(1:N,1)+mu2(C(1:N)).*F(1:N,2))/diffuse_C;
    w = [0;D2\DC]+S(1);
end

%% Solve Boundary value problem
% Boundary value numerical parameters
tol = 1E-5;
% Simple fixed point solver for the BVP
function w = BVPSolve(C,S,F)
    change = 1; 
    while change>tol
        Cnew = BVP(C,S,F);
        change = norm(C-Cnew,inf);
        C = Cnew;
    end
    w = C;
end

%% timestepping loop (ETDRK4)
function [w,ww] = ETDRK4(C,S,F)
        NSu = planktonic(C,S,F);
        Nu = transport(C,S,F);
        au = F+QU*Nu;
        as = S+QS*NSu;
        NSa = planktonic(C,as,au);
        Na = transport(C,as,au);
        bu = F+QU*Na;
        bs = S+QS*NSa;
        NSb= planktonic(C,bs,bu);
        Nb = transport(C,bs,bu);
        cu = au+QU*(2*Nb-Nu);
        cs = as+QS*(2*NSb-NSu);
        NSc = planktonic(C,cs,cu);
        Nc = transport(C,cs,cu);
        w = S+FS1*NSu+2*FS2*(NSa+NSb)+FS3*NSc;
        ww = F+fU1*Nu+2*fU2*(Na+Nb)+fU3*Nc;
end

%% Set Storage
SS = zeros(4,points); SS(:,1)= S(:);
tt = zeros(1,points);
f1f1 = zeros(N+1,points); f1f1(:,1) = F(:,1);
f2f2 = zeros(N+1,points); f2f2(:,1) = F(:,2);
CC = zeros(N+1,points); CC(:,1) =  C;
reverseStr = '';

%% Main Solver
for i = 1:nplots
    % Save into temporary Storage
    t= t+dt; 
    C = BVPSolve(C,S,F);
    [S,F] = ETDRK4(C,S,F);
    %% Save the data sometimes 
    if mod(i,plotgap)==0
      j = i/plotgap+1;
      tt(j) = t ; SS(:,j) = S;
      f1f1(:,j) = F(:,1); f2f2(:,j) = F(:,2);
      CC(:,j) = C; 
      percentDone = 100 * i/nplots;
      msg = sprintf('Percent done: %3.3f', percentDone); %Don't forget this semicolon
      fprintf([reverseStr, msg]);
      reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
  end
fprintf([reverseStr, 'Percent done: 100', '\n']);
x = (1+x)/2;
end


