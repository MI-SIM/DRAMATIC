%% Single species biofilm without considering the stirred community
% Using the model of Pritchett, Dockery (2001)
% Equations:
% dS/dt = d2S/dx2 - 1/phi^2 S/(1+S)u
% du/dt = -d/dx(vu) + (S/(1+S)-alpha-beta)*u
% dL/dt = u(L(t),t)
% v(z,t)= int((S(y,t)/(1+S(y,t))-beta)u(y,t)dy
%
% Boundary Conditions
% S'(0,t) =  0
% S(L(t),t) = S_b
% f'(0,t) = 0

%% Parameters
A = ; au = ; aw =;
d = ; dc =; E = ;
F = ; 


% % Scaled Parameters

Lhat = sqrt(d/mu);  phi = sqrt(Y*K/rho);
alpha = Ki/mu; beta = Kd/mu;
Sb = Cb/K;

%% Discretization of domain; 
N = 64;                                                 % # of Chebyshev Points
[D,x] = cheb(N);                                        % Chebyshev Differentiation Matrix
Dx= D; Dx(N+1,:) = zeros(N+1,1);                        %
mu = 1.65E-5; d = 1E-5;
rho = 2.5E-3; Y = 0.071;
K = 3.9E-3; Ki = 1E-6;                                  % For the Neumann B.C.
D2 = D^2; D2(N+1,:) = D(N+1,:); D2 = D2(2:N+1,2:N+1);   % Also for the Neumann B.C.
Drde =  Dx(1:N,1:N);
eps = 0.01; dt = min([.01,N^(-4)/eps]);                 % timestep dependent on N
t=0;
int = clencurt(N);                                      % Clenshaw Curtis Integration

%% Initial Conditions
global u s L v                   % global parameters so the functions update

 u = 0.1*ones(size(x));
 s = Sb/4*(x+1).^2;              %% naive initial guess for s
 L = Lhat;


%% Numerical Settings
tmax= 10; tplot = 1;            % Maxtime and plot times.
nplots = round(tmax/dt);        % # of iterations for timestepping
plotgap = round(tplot/dt);      % Gap between saving points for plotting
dt = tplot/plotgap;             % Timestep
tol = 1E-6;                     % Tolerance for the BVP

%% Differential Equations
% Rescaled to remove the moving boundary
% Dx becomes: 2/L*Dz
% D2x becomes: 4/L^2*D2z
% Dt becomes: Dt-(z+1)*L'/L*Dz

v = int.*((s./(1+s)-beta).*u); v = v-v(end);                  % Velocity Calculation
transport = @(y) -2*v.*Dx*y+v(1)*(x+1).*Dx*y+(s./(1+s)-beta).*y.*(1-y)-alpha*y;  
RDE = @(y) L^2/phi^2*y./(1+y).*u(2:N+1);
thick =  @(y) y*v(1);

%% Correct the initial guess for s:
change = 1;
while change > tol
        Snew = [0;D2\RDE(s(2:N+1))]+Sb; %% +Sb for the right hand boundary value
        change = norm(Snew-s,inf);
        s = Snew;
end


%% Solving
uu = u'; ss = s';
LL = L;  tt = t;
for i = 1:nplots
    change = 1;
    Linop = D2 + 1/2*L*v(1)*(x(1:N)+1).*Drde;
    v = int.*((s./(1+s)-beta).*u); v = v-v(end);    % re-calculate the velocity
    t = t+dt;
    ut = RK4(transport,u,dt);                   % RK4 is a 4th order Runge-Kutte solver
    Lt = RK4(thick,L,dt);
    while change > tol                          %% Solve BVP for S
        Snew = [0;Linop\RDE(s(2:N+1))]+Sb;         %% +Sb for the RHS B.C.
        change = norm(Snew-s,inf);
        s = Snew;
    end
    u = ut;       L = Lt;
    if mod(i,plotgap)==0                        % Only save the data sometimes
       LL = [LL;L];  tt = [tt;t];
       uu = [uu;u'];  ss= [ss;s'];
    end
end
x = 1/2*(1+x);                  %% Rescale the plots so they are correct.
clf 
figure(1)
surf(x,tt,uu), lighting phong, axis tight
view([-45 60]), colormap(cool), light('col',[1 1 0],'pos',[-10 0 10])
xlabel  'x', ylabel 't', zlabel 'u', axis([ 0 1 0 tmax 0 1]);
figure(2)
plot(x,s);
xlabel 'x', ylabel 's'
figure(3)
plot(tt,LL)
xlabel 't', ylabel 'L'
