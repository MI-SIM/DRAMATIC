%ODE solver for DM5.Laureni

%% Set parameters and impulse
%options = odeset('RelTol',10e-6,'AbsTol',10e-6,'OutputFcn',@odeplot); %Solver options
init = [20,0,0,0,300,300]; %Initial conditions

%Nitrogen incorporation
iNa = 0;
iNn = 0;

%AMX concentration dependent on rAMX
rAmx = 86; %mg/L/d
muAmx = 0.017; %1/d
XAmx = rAmx/muAmx; %mg/L

%DO concentration
DO = 1.5;

%Biomass retention
WAS = 5e-3;
Bret = 1-WAS;

params = [0.297,0.337,0.017,2.4,0.5,0.6,0.4,0.03,0.005,0.18,0.08,0.17,0.083,0.083,0.058,iNa,iNn,DO,XAmx]; %Growth parameters and fixed state values (O2 and AMX)
tspan = [0,400];

s1in = 20;
s2in = 0;
s3in = 0;
s4in = 0;
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
plot(t0,ximp(:,1:4),'.','markersize',4)
xlabel('Time (d)')
ylabel('Substrate concentrations (mgN L^{-1})')
legend({'NH_4^+','NO_2^-','NO_3^-','N_2'})
figure(4)
plot(t0,ximp(:,5:6),'.','markersize',4)
xlabel('Time (d)')
ylabel('Biomass concentrations (mgCOD L^{-1})')
legend({'AOB','NOB'})
