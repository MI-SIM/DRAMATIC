%1D hybrid MBBR model describing Partial Nitritation/Anammox with
%controllable oxygen and SRT 
%M. Wade (August, 2018), after Masic & Eberl (2014)

%Call parameters
hybrid_parameters

%Set solver options
options=odeset('reltol',1e-6,'abstol',1e-6,'stats','off','outputFcn',[]);

%Run solver
Bfini=[so2,0,snh4,0,sno2,0,sno3,0,sn2,0];
init=[so2,snh4,sno2,sno3,sn2,uA,uN,uX,uI,L0,fi_init_rs,xpos];
[tout,yout] = ode45(@DM6wade,[0:0.1:90],init,options,Sins,D,Ds,T,Bfini);
