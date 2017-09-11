%ODE solver for DM5.Laureni

options = odeset('RelTol',10e-6,'AbsTol',10e-6,'OutputFcn',@odeplot); %Solver options
init = [20,0,0,0,1e-8,1e-8]; %Initial conditions

params = [0.297,0.337,0.017,2.4,0.5,0.6,0.4,0.03,0.005,0.18,0.08,0.17,0.083,0.083,0.058,0.014,0.02,1,700]; %Growth parameters and fixed state values (O2 and AMX)

t = linspace(0,1000,100); %Time vector

%ODE solver
[tout,yout] = ode23s(@DM5_growth,t,init,options,params);