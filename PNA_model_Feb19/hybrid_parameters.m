%Parameters file for DM6.Wade

%Inputs
SOin = 15;
Sins = [20,SOin]; %Ammonia and oxygen input into reactor (mg/L)

%Kinetic growth parameters
    %AOB
muA = 0.297; %Max. specific growth rate of AOB (1/d)
YA = 0.18; %Growth yield of AOB (gCOD/gN)
KamA = 2.4; %Ammonium half-saturation constant of AOB (gNH4-N/m3)
Ko2A = 0.6; %Oxygen half-saturation constant of AOB (gCOD/m3)
bA = 0.014; %Decay rate of AOB (1/d)
iNA = 0.083; %Nitrogen incorporation in AOB (gN/gCOD)

    %NOB
muN = 0.337; %Max. specific growth rate of NOB (1/d)
YN = 0.08; %Growth yield of NOB (gCOD/gN)
KniN = 0.5; %Ammonium half-saturation constant of NOB (gNO2-N/m3)
Ko2N = 0.4; %Oxygen half-saturation constant of NOB (gCOD/m3)
bN = 0.02; %Decay rate of NOB (1/d)
iNN = 0.083; %Nitrogen incorporation in NOB (gN/gCOD)

    %AMX
muX = 0.017; %Max. specific growth rate of AMX (1/d)
YX = 0.17; %Growth yield of AMX (gCOD/gN)
KamX = 0.03; %Ammonium half-saturation constant of AMX (gNH4-N/m3)
KniX = 0.005; %Nitrite half-saturation constant of AMX (gNO2/m3)
bX = 0.0003; %Decay rate of AMX (1/d)
iNX = 0.058; %Nitrogen incorporation in AMX (gN/gCOD)
KIamx = 0.01; %Oxygen affinity coefficient on AMX (gO2/m3) [Van Hulle, 2005, Diss.; Trojanowicz et al., 2017, Env. Tech.)

%Spatial Discretisation
T = 200; %Number of spatial discretisation points
L0 = 1e-5; %Initial biofilm height (m)
xpos = linspace(0,L0,T); %Initial grid positions

%Reactor
HRT=5.9;
SRT=24.5;
D = 1/(HRT*24); %Dilution rate (1/d) (HRT = 1/D*24 h)
Ds = 1/SRT; %Solids washout rate (1/d) (SRT = 1/Ds d)
%nc = 50; %Number of carriers
%Ac = 0.0068; %Carrier area (m2)
%A = nc*Ac + 0.17; %Total surface area with nc carriers (m2) 
V = 0.012; %reactor volume - 12L bulk liquid volume (m3)
cV = V/3; %Volumetric loading of K5 carriers (m3)
A = cV*800; %Total surface area of carriers (m2)

%Physical parameters
alpha = 1; %Attachment rate (1/d)
E = 1000; %Erosion parameter (1/m.d)
Do2 = 9.93e-5; %O2 diffusion coefficient (m2/d) 
Dnh4 = 9.13e-5; %NH4 diffusion coefficient (m2/d) 
Dno2 = 7.93e-5; %NO2 diffusion coefficient (m2/d)
Dno3 = 7.93e-5; %NO3 diffusion coefficient (m2/d) 
Dn2 = 1.73e-4; %N2 diffusion coefficient (m2/d) - Cadogan et al. (2014) J. Chem. Eng. Data

%Initial conditions
    %Chemical species
so2 = SOin; %Oxygen
snh4 = 20; %Ammonia
sno2 = 0; %Nitrite
sno3 = 1; %Nitrate
sn2 = 1; %Nitrogen gas

    %Planktonic species
uA = 0.003; %AOB (mgCOD/L)
uN = 0.003; %NOB (mgCOD/L)
uX = 0;    %AMX (mgCOD/L)
uI = 0;    %Inerts (mgCOD/L)

    %Biofilm fractions
fa = 0.1; %Initial fraction of AOB in biofilm
fn = 0.1; %Initial fraction of NOB in biofilm 
fx = 0.8; %Initial fraction of AMX in biofilm 
fxi = 0.1; %Fraction of inert biomass produced by endogenous respiration
rho = 10000; %Biomass density (g/m3)
eta = 0.5; %Anoxic reduction factor
etaX = 0.1; %CHECK - Nitrite reduction factor for AMX

%Initial fraction of species in biofilm
%Initial AMX fraction = 1 upto 20% below surface
%Weibull distribution dictates relative fractions in remaining 20% 

nf = 1; %Fractional depth of AOB/NOB
AMXf = ones(T*nf,1)*0.95;
AOBf = ones(T*nf,1)*0.025;
NOBf = ones(T*nf,1)*0.025;
% Weib_X = 1:0.2*T;
% pdfW = wblpdf(Weib_X,50,0.3); %Shallow slope
% scale_pdfW = pdfW./max(pdfW);

% AMXf100 = [AMXf80;scale_pdfW'];
% AOBf100 = [AOBf80;(1-scale_pdfW')/2];
% NOBf100 = [NOBf80;(1-scale_pdfW')/2];

fi_init = [AOBf,NOBf,AMXf,zeros(T,1)];
fi_init_rs = reshape(fi_init,1,T*4);

