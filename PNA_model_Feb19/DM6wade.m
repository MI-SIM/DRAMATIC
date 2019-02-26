%DM6.Wade - suspended flocs & wall attachment, multiple species (June, 2018)
function [dP] = DM6wade(t,y,Sins,D,Ds,Nx,Bfguess1)

persistent Bfguess

%% Solver iterations
%Inputs
So2 = y(1); %Oxygen
Snh4 = y(2); %Ammonia
Sno2 = y(3); %Nitrite
Sno3 = y(4); %Nitrate
Sn2 = y(5); %N2
uA = y(6); %Planktonic AOB
uN = y(7); %Planktonic NOB
uX = y(8); %Planktonic AMX
uI = y(9); %Planktonic inert
L = y(10); %Biofilm thickness
fA = y(11:10+Nx);
fN = y(11+Nx:10+(2*Nx));
fX = y(11+(2*Nx):10+(3*Nx));
fI = y(11+(3*Nx):10+(4*Nx));
xpos = y(11+(4*Nx):end);

fi_old=[fA,fN,fX,fI];

%% Parameters

%Kinetic growth parameters
    %AOB
muA = 0.297; %Max. specific growth rate of AOB (1/d)
YA = 0.18; %Growth yield of AOB (gCOD/gN)
KamA = 2.4; %Ammonium half-saturation constant of AOB (gNH4-N/m3)
Ko2A = 0.6; %Oxygen half-saturation constant of AOB (gCOD/m3)
bA = 0.014; %Decay rate of AOB (1/d)
iNA = 0.083; %Nitrogen content in AOB (gN/gCOD)

    %NOB
muN = 0.337; %Max. specific growth rate of NOB (1/d)
YN = 0.08; %Growth yield of NOB (gCOD/gN)
KniN = 0.5; %Nitrite half-saturation constant of NOB (gNO2-N/m3)
Ko2N = 0.4; %Oxygen half-saturation constant of NOB (gCOD/m3)
bN = 0.02; %Decay rate of NOB (1/d)
iNN = 0.083; %Nitrogen content in NOB (gN/gCOD)

    %AMX
muX = 0.017; %Max. specific growth rate of AMX (1/d)
YX = 0.17; %Growth yield of AMX (gCOD/gN)
KamX = 0.03; %Ammonium half-saturation constant of AMX (gNH4-N/m3)
KniX = 0.005; %Oxygen half-saturation constant of AMX (gNO2/m3)
bX = 0.0003; %Decay rate of AMX (1/d)
iNX = 0.058; %Nitrogen content in AMX (gN/gCOD)
KIamx = 0.01; %Oxygen affinity coefficient on AMX (gO2/m3) [Van Hulle, 2005, Diss.; Trojanowicz et al., 2017, Env. Tech.)

%%%CONSIDER Nitrite inhibition of AMX too?

    %Inerts
iNI = 0.02; %Nitrogen content in inerts (gN/gCOD)

%Reactor
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

fxi = 0.1; %Fraction of inert biomass produced by endogenous respiration
rho = 10000; %Biomass density (g/m3)
eta = 0.5; %Anoxic reduction factor
etaX = 0; %CHECK - Nitrite reduction factor for AMX

disp(['Time = ',num2str(t),' (d)'])

%Influent conditions
Snh4_in = Sins(1);
So2_in = Sins(2);
Sno2_in = 0;
Sno3_in = 0;
Sn2_in = 0;
yL=[];

%Parameter vectors
Dk = [Do2,Dnh4,Dno2,Dno3,Dn2]; %Diffusion coefficients
S = [So2,Snh4,Sno2,Sno3,Sn2];  %Substrate concentrations in bulk liquid

%Boundary value problem
x_init = xpos; %Discretisation
if t==0
    yi = Bfguess1;
    solinit = bvpinit(x_init, yi); %Initial guess
else
    solinit.solver='bvpinit';
    solinit.x = x_init;
    solinit.y = Bfguess;
    solinit.yinit = Bfguess(:,end)';
end

%BVP settings
nmax = 5000;
Xx = solinit.x; %Discretization calculated based on ICs

options = bvpset('RelTol',1.0E-4,'Nmax',nmax);
try
    sol = bvp5c( @twoode, @twobc, solinit, options);
catch
    sol = bvp4c( @twoode, @twobc, solinit, options);
end

%Fraction of biomass types in biofilm
fA = fi_old(:,1);
fN = fi_old(:,2);
fX = fi_old(:,3);
fI = fi_old(:,4);

%Check that fractions sum to 1
if std([sum(fA+fN+fX+fI)-200,0])>1e-12
    warning('fractions don''t sum to one')
    dbstop if warning
end

%Evaluate solutions at each discretisation point
yL = deval(sol,x_init);

%Substrates in biofilm (Cbar)
Co2 = yL(1,:)';
Cnh4 = yL(3,:)';
Cno2 = yL(5,:)';
Cno3 = yL(7,:)';
Cn2 = yL(9,:)';

BF_conc = yL([1 3 5 7 9],:)';

Bfguess = yL;

%Growth functions in biofilm
muOaBF = Co2./(Ko2A+Co2); %Growth of AOB on oxygen
muOnBF = Co2./(Ko2N+Co2); %Growth of NOB on oxygen
muNHBF = Cnh4./(KamA+Cnh4); %Growth of AOB on NH4
muNHBFX = Cnh4./(KamX+Cnh4); %Growth of AMX on NH4
muNO2BF = Cno2./(KniN+Cno2); %Growth of NOB on NO2
muNO2BFX = Cno2./(KniX+Cno2); %Growth of AMX on NO2

%Anammox inhibition by oxygen
AMX_IBF = KIamx./(KIamx + Co2);

dsAOBBF = muA.*muOaBF.*muNHBF; %Dual substrate limitation (oxygen and ammonia)
dsNOBBF = muN.*muOnBF.*muNO2BF; %Dual substrate limitation (oxygen and nitrite)
dsAMXBF = muX.*muNHBFX.*muNO2BFX.*AMX_IBF; %Dual substrate limitation (ammonia and nitrite) and oxygen inhibition

muA_NBF = dsAOBBF - muA.*muOaBF.*(1+eta)*bA;
muN_NBF = dsNOBBF - muN.*muOnBF.*(1+eta)*bN;
muX_NBF = dsAMXBF - muX.*AMX_IBF.*(1+etaX)*bX; %CHECK THIS

muIBF(fI==0,1)=1e20;
muIBF(fI~=0,1)=(1./fI(fI~=0)).*(fxi + eta).*(fA(fI~=0).*bA.*muOaBF(fI~=0) + fN(fI~=0).*bN.*muOnBF(fI~=0))+(1./fI(fI~=0)).*(fxi+etaX).*(fX(fI~=0).*bX.*AMX_IBF(fI~=0));

mubar0=[muA_NBF,muN_NBF,muX_NBF,muIBF];

mu1a = muOaBF;
mu2a = muOnBF;
mu3a = muNHBF;
mu4a = muNHBFX;
mu5a = muNO2BF;
mu6a = muNO2BFX;

dsONHa = muA.*mu1a.*mu3a; %Dual substrate limitation (oxygen and ammonia)
dsONO2a = muN.*mu2a.*mu5a; %Dual substrate limitation (oxygen and nitrite)
dsNHNO2a = muX.*mu4a.*mu6a.*AMX_IBF; %Dual substrate limitation (ammonia and nitrite) and oxygen inhibition

RT1a = ((3.43 - YA)/YA).*dsONHa.*fA.*rho + (1 - fxi).*bA.*mu1a.*fA.*rho + ((1.14-YN)/YN).*dsONO2a.*fN.*rho + (1 - fxi).*bN.*mu2a.*fN.*rho;
RT2a = (1/YA + iNA).*dsONHa.*fA.*rho - (iNA - iNI*fxi).*bA.*mu1a.*fA.*rho + iNN.*dsONO2a.*fN.*rho - (iNN - iNI*fxi).*bN.*mu2a.*fN.*rho + (1/YX + iNX).*dsNHNO2a.*fX.*rho - (iNX - iNI*fxi).*bX.*AMX_IBF.*fX.*rho;
RT3a = (-1/YA).*dsONHa.*fA.*rho + (1/YN).*dsONO2a.*fN.*rho + (1/YX + 1/1.14).*dsNHNO2a.*fX.*rho;
RT4a = (-1/YN).*dsONO2a.*fN.*rho - (1/1.14).*dsNHNO2a.*fX.*rho;
RT5a = (-2/YX).*dsNHNO2a.*fX.*rho;

%Substrate flux
if  L==0
    j=zeros(4,1);
else
    j = (1./Dk).*trapz(xpos,[RT1a,RT2a,RT3a,RT4a,RT5a]);
end

%Suspended biomass growth functions
muOA = So2/(Ko2A+So2); %Growth of AOB on oxygen
muON = So2/(Ko2N+So2); %Growth of NOB on oxygen
muNH4A = Snh4/(KamA+Snh4); %Growth of AOB on NH4
muNH4X = Snh4/(KamX+Snh4); %Growth of AMX on NH4
muNO2N = Sno2/(KniN+Sno2); %Growth of NOB on NO2
muNO2X = Sno2/(KniX+Sno2); %Growth of AMX on NO2

%Anammox inhibition by oxygen
AMX_I = KIamx/(KIamx + So2);

dsAOBPL = muA*muOA*muNH4A; %Dual substrate limitation (oxygen and ammonia)
dsNOBPL = muN*muON*muNO2N; %Dual substrate limitation (oxygen and nitrite)
dsAMXPL = muX*muNH4X*muNO2X*AMX_I; %Dual substrate limitation (ammonia and nitrite) and oxygen inhibition

muA_N1A = dsAOBPL - muA*muOA*(1+eta)*bA;
muN_N1A = dsNOBPL - muN*muON*(1+eta)*bN;
muX_N1A = dsAMXPL - muX*AMX_I*(1+etaX)*bX; %CHECK THIS

%Bulk liquid reaction rates (CHECK DECAY RATES - just bX?)
ro2 = ((3.43 - YA)/YA)*dsAOBPL*uA + (1 - fxi)*bA*muOA*uA + ((1.14 - YN)/YN)*dsNOBPL*uN + (1 - fxi)*bN*muON*uN; 
rnh4 = (1/YA + iNA)*dsAOBPL*uA - (iNA - iNI*fxi)*bA*muOA*uA + iNN*dsNOBPL*uN - (iNN - iNI*fxi)*bN*muON*uN + (1/YX + iNX)*dsAMXPL*uX - (iNX - iNI*fxi)*bX*AMX_I*uX;
rno2 = (-1/YA)*dsAOBPL*uA + (1/YN)*dsNOBPL*uN + (1/YX + 1/1.14)*dsAMXPL*uX;
rno3 = (-1/YN)*dsNOBPL*uN - (1/1.14)*dsAMXPL*uX;
rn2 = (-2/YX)*dsAMXPL*uX;

%Substrate ODEs with diffusion
dSo2dt = D*(So2_in - So2) - V^-1*(ro2 + A*Do2*j(1)); %Oxygen
dSnh4dt = D*(Snh4_in - Snh4) - V^-1*(rnh4 + A*Dnh4*j(2)); %Ammonium
dSno2dt = D*(Sno2_in - Sno2) - V^-1*(rno2 + A*Dno2*j(3)); %Nitrite
dSno3dt = D*(Sno3_in - Sno3) - V^-1*(rno3 + A*Dno3*j(4)); %Nitrate
dSn2dt = D*(Sn2_in - Sn2) - V^-1*(rn2 + A*Dn2*j(5)); %N2

%Detachment from biofilm
d = E*L;

duAdt = uA*(muA_N1A - Ds - alpha) + A*rho*d*L*fi_old(end,1);
duNdt = uN*(muN_N1A - Ds - alpha) + A*rho*d*L*fi_old(end,2);
duXdt = uX*(muX_N1A - Ds - alpha) + A*rho*d*L*fi_old(end,3);
duIdt = (fxi + eta)*(uA*bA*muOA + uN*bN*muON) + (fxi + etaX)*(uX*bX*AMX_I) - uI*Ds - alpha*uI + A*rho*d*L*fi_old(end,4); 
exf = trapz(xpos',sum(fi_old.*mubar0,2));
dLdt = exf + (alpha/(A*rho))*(uA + uN + uX + uI) - d*L*sum(fi_old(end,:));

%Biomass fractions in biofilm differentials dfi/dt
%Calculate biofilm fractions 

% invA = (alpha/(A*rho))*uA;
% invN = (alpha/(A*rho))*uN;
% invX = (alpha/(A*rho))*uX;
% invI = (alpha/(A*rho))*uI;
% invasion = [invA,invN,invX,invI];
% 
% detA = d*L*sum(fi_old(end,1));
% detN = d*L*sum(fi_old(end,2));
% detX = d*L*sum(fi_old(end,3));
% detI = d*L*sum(fi_old(end,4));
% detachment = [detA,detN,detX,detI];

v(1) = 0;
lxp = length(xpos);

for kk = 2:lxp  
    v(kk) = trapz(xpos(1:kk)',sum([fi_old(1:kk,:).*mubar0(1:kk,:)],2));
end

newxpos = v';

% plot(t,L,'r.','markersize',12)
% hold on
% plot(t,xpos(end),'k.')
% drawnow
plot(t,uX,'k.','markersize',20)
hold on
drawnow
dvdx = sum(fi_old.*mubar0,2);
dfi_new = fi_old.*mubar0 - fi_old.*dvdx;

%Biomass fractions in biofilm differentials dfi/dt
bfA = dfi_new(:,1);
bfN = dfi_new(:,2);
bfX = dfi_new(:,3);
bfI = dfi_new(:,4);

%Write final fractions to file for plotting
if t==90
     dlmwrite('fractions.txt',fi_old)
     dlmwrite('conc_bf.txt',BF_conc)
end

dP = [dSo2dt;dSnh4dt;dSno2dt;dSno3dt;dSn2dt;duAdt;duNdt;duXdt;duIdt;dLdt;bfA;bfN;bfX;bfI;newxpos];

    function dydx = twoode(x,ny)

        %Growth functions in planktonic phase
        %Find closest node in discretised space equal to adaptive BVP grid
        %discretisation and use this index for fractional biomass in solver
        [minm,indz] = min(abs(Xx-x));

        if isempty(indz)
            dydx=[0;0;0;0;0;0;0;0;0;0];
            dydx=0;
        else
            S1 = ny(1); %O2
            S2 = ny(3); %NH4
            S3 = ny(5); %NO2
            
            mu1 = S1/(Ko2A+S1); %Growth of AOB on oxygen
            mu2 = S1/(Ko2N+S1); %Growth of NOB on oxygen
            mu3 = S2/(KamA+S2); %Growth of AOB on NH4
            mu4 = S2/(KamX+S2); %Growth of AMX on NH4
            mu5 = S3/(KniN+S3); %Growth of NOB on NO2
            mu6 = S3/(KniX+S3); %Growth of AMX on NO2
            
            %Anammox inhibition by oxygen
            AMXI = KIamx/(KIamx+S1);
            
            dsONH = muA*mu1*mu3; %Dual substrate limitation (oxygen and ammonia)
            dsONO2 = muN*mu2*mu5; %Dual substrate limitation (oxygen and nitrite)
            dsNHNO2 = muX*mu4*mu6*AMXI; %Dual substrate limitation (ammonia and nitrite) and oxygen inhibition

            RT1 = ((3.43 - YA)/YA)*dsONH*fA(indz)*rho + (1 - fxi)*bA*mu1*fA(indz)*rho + ((1.14-YN)/YN)*dsONO2*fN(indz)*rho + (1 - fxi)*bN*mu2*fN(indz)*rho;
            RT2 = (1/YA + iNA)*dsONH*fA(indz)*rho - (iNA - iNI*fxi)*bA*mu1*fA(indz)*rho + iNN*dsONO2*fN(indz)*rho - (iNN - iNI*fxi)*bN*mu2*fN(indz)*rho + (1/YX + iNX)*dsNHNO2*fX(indz)*rho - (iNX - iNI*fxi)*bX*AMXI*fX(indz)*rho;
            RT3 = (-1/YA)*dsONH*fA(indz)*rho + (1/YN)*dsONO2*fN(indz)*rho + (1/YX + 1/1.14)*dsNHNO2*fX(indz)*rho;
            RT4 = (-1/YN)*dsONO2*fN(indz)*rho - (1/1.14)*dsNHNO2*fX(indz)*rho;
            RT5 = (-2/YX)*dsNHNO2*fX(indz)*rho;

            ny2 = [RT1,RT2,RT3,RT4,RT5]./Dk;

            dydx = [ny(2);ny2(1);ny(4);ny2(2);ny(6);ny2(3);ny(8);ny2(4);ny(10);ny2(5)];
        end
    end

    function res = twobc(ya,yb)
        %Residual function for BVP solver
        yb2=[yb(1);yb(3);yb(5);yb(7);yb(9)] - S';
        res=[ya(2);yb2(1);ya(4);yb2(2);ya(6);yb2(3);ya(8);yb2(4);ya(10);yb2(5)];
    end

end




