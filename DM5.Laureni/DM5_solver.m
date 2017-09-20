%ODE solver for DM5.Laureni

%% Set parameters and impulse
%options = odeset('RelTol',10e-6,'AbsTol',10e-6,'OutputFcn',@odeplot); %Solver options
init = [20,0,0,0,300,300]; %Initial conditions

%Decay rates
bA = 0.014;
bN = 0.02;

%Nitrogen incorporation
iNa = 0.083;
iNn = 0.083;
iNx = 0.058;

%AMX concentration dependent on rAMX
rAmx = 86; %mgCOD/L/d
muAmx = 0.017; %1/d
XAmx = rAmx/muAmx; %mgCOD/L

%DO concentration
DO = 1.15; %mgO2/L

%Biomass retention
WAS = 4e-3;
Bret = 1-WAS;

params = [0.297,0.337,0.017,2.4,0.5,0.6,0.4,0.03,0.005,0.18,0.08,0.17,iNa,iNn,iNx,bA,bN,DO,XAmx]; %Growth parameters and fixed state values (O2 and AMX)
tspan = [0,150];

s1in = init(1);
s2in = init(2);
s3in = init(3);
s4in = init(4);
s1bar = 2;
imp =@(s1) s1 > s1bar;
xin = [s1in s2in s3in s4in];
g =@(x) [(x(1:4)+xin)/2 x(5:6)*Bret];

%% Run solver
[t,x,t0,x0,ximp] = ImpulseA(@(t,x)DM5_growth(t,x,params),g,imp,init,tspan);

%% Plot dynamics
%Post-impulse
figure(1)
plot(t0,x0(:,1:4),'.','markersize',4)
xlabel('Time (d)')
ylabel('Substrate concentrations (mgN L^{-1})')
legend({'NH_4^+','NO_2^-','NO_3^-','N_2'})
figure(2)
plot(t0,x0(:,5:6),'.','markersize',4)
xlabel('Time (d)')
ylabel('Biomass concentrations (mgCOD L^{-1})')
legend({'AOB','NOB'})

%Pre-impulse
figure(3)
plot(t0,ximp(:,1:4),'-','markersize',4)
xlabel('Time (d)')
ylabel('Substrate concentrations (mgN L^{-1})')
legend({'NH_4^+','NO_2^-','NO_3^-','N_2'})
figure(4)
plot(t0,ximp(:,5:6),'-','markersize',4)
xlabel('Time (d)')
ylabel('Biomass concentrations (mgCOD L^{-1})')
legend({'AOB','NOB'})
