function [tt,SS,CNH4C,CNO2C,CO2C,fafa,fnfn,fifi,x] = PartialNitritation_NO3(T,N)
%% Anammox and Partial Nitritation in a biofilm
% Using the model of Masic and Eberl (2015)
% Note that this is a function, not a script so that we can use 
% nested functions to simplify the use of parameters. 


%% Parameters
Ar = 0.17; d =0.5 ; rho = 10000;
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
ia = 0; ii = 0;


%% Discretization of domain; 
[D,x] = cheb(N);                                        % Chebyshev Differentiation Matrix
Dx= D; Dx(N+1,:) = zeros(N+1,1);                        % For the Neumann B.C. x(end) -> z=0
D2 = D^2; D2(N+1,:) = D(N+1,:); D2 = D2(2:N+1,2:N+1);   % Also for the Neumann B.C.
eps = 0.01; %dt = min([10e-6,N^(-4)/eps]);              % timestep dependent on N
dt=0.001;
t=0;
% Chebyshev integration matrix:
% Sl is the indefinite integral from x to 1,
% int is the indefinite integral from -1 to x
[int,~] = intmat(N); 


%% Numerical Settings
tmax= T; points=900;           % Maxtime and number of save points.
tplot = T/points;               % plot times
nplots = round(tmax/dt);        % # of iterations for timestepping
plotgap = round(tplot/dt);      % Gap between saving points for plotting
dt = tplot/plotgap;             % Timestep
tol = 1E-9;                     % Tolerance for the BVP

%% Exponential time differencing setup
% For Transport PDE
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
% For the ODEs:
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

%% Response Functions & Respiration Rates
muaO = @(y) y./(kaO2+y); % Oxygen consumption by AOB
munO = @(y) y./(knO2+y); % Oxygen consumption by NOB
mua = @(y,z) muamax*y./(kaNH4+y).*muaO(z); % Ammonium consumption by AOB
mun = @(y,z) munmax*y./(knNO2+y).*munO(z); % Nitrite consumption by NOB
Rxa = @(y,z) mua(y,z)-muamax*(1+eta)*ba*muaO(z); % Called mu_a in the paper
Rxn = @(y,z) mun(y,z)-munmax*(1+eta)*bn*munO(z); % Called mu_n in the paper
mui = @(xa,xn,y) (fxi+eta)*(ba*xa.*muaO(y)+bn*xn.*munO(y)); % Reaction term in inert biomass equation
Rf = @(C,F) F(:,1).*Rxa(C(:,1),C(:,3))+F(:,2).*Rxn(C(:,2),C(:,3))+mui(F(:,1),F(:,2),C(:,3)); % Sum of biomass reactions


%% Biofilm Velocities and transport flux
function w = v(C,~,F) % biofilm velocity 
    w = int*(Rf(C,F));
end
function w = Lp(C,S,F) % L prime
        vel = v(C,S,F);
        w = S(6)/2*vel(1)+alpha/(A*rho)*(S(3)+S(4)+S(5))-E*S(6)^2;
end
%% Biofilm Substrate Update Function
% Biofilm Substrate equations
function w = RDE(C,~,F)
             dz = (S(6)/2)^2;
        % Terms in NH4 equation
             RNH4a = (1/ya+ia)*mua(C(2:N+1,1),C(2:N+1,3));
             RNH4b = (ia-ii*fxi)*ba*muaO(C(2:N+1,3));
             RNH4n = ia*mun(C(2:N+1,2),C(2:N+1,3));
             RNH4m = (ia-ii*fxi)*bn*munO(C(2:N+1,3));
         RDENH4 =  (dz/dNH4)*((RNH4a-RNH4b).*F(2:N+1,1)*rho+(RNH4n-RNH4m).*F(2:N+1,2)*rho);
        % NH4 equation, no i_i's
%         RDENH4 = (1/dNH4)*(1/ya)*mua(C(2:N+1,1),C(2:N+1,3)).*F(2:N+1,1)*rho;
        % Terms in NO2 equation
            RNO2a = 1/ya*mua(C(2:N+1,1),C(2:N+1,3));
            RNO2n = 1/yn*mun(C(2:N+1,2),C(2:N+1,3));
        RDENO2 = (dz/dNO2)*(RNO2n.*F(2:N+1,2)*rho-RNO2a.*F(2:N+1,1)*rho);
        % Terms in the oxygen equation
             RO2a = (3.43-ya)/ya*mua(C(2:N+1,1),C(2:N+1,3));
             RO2b = (1-fxi)*ba*muaO(C(2:N+1,3));
             RO2n = (1.14-yn)/yn*mun(C(2:N+1,2),C(2:N+1,3));
             RO2m = (1-fxi)*bn*munO(C(2:N+1,3));
        RDEO2 = (dz/dO2)*((RO2a+RO2b).*F(2:N+1,1)*rho+(RO2n+RO2m).*F(2:N+1,2)*rho);
        w = [RDENH4,RDENO2,RDEO2];
end
%% QR delete;
    function [Q,R] = qrdelete(Q,R)
        m = size(R,1);
        n = size(Q,1);
        for kk = 1:m-1
            temp = sqrt(R(kk,kk+1)^2+R(kk+1,kk+1)^2);
            c = R(kk,kk+1)/temp; s = R(kk+1,kk+1)/temp;
            R(kk,kk+1)=temp; R(kk+1,kk+1)=0;
            if kk <m-1
                for jj = kk+2:m
                    temp = c*R(kk+1,jj)+s*R(kk+1,jj);
                    R(kk+1,jj)=-s*R(kk,jj)+c*R(kk+1,jj);
                    R(kk,jj)=temp;
                end
            end
            for ll=1:n
                temp = c*Q(ll,kk)+s*Q(ll,kk+1);
                Q(ll,kk+1)=-s*Q(ll,kk)+c*Q(ll,kk+1);
                Q(kk,ll) = temp;
            end
        end
        Q=Q(:,1:m-1); R =R(1:m-1,2:m);
    end
%% Fixed point iteration using Andersen Acceleration 
function w = biofilmbvp(C,S,F)
    %CNH4 = C(:,1); CNO2 = C(:,2); CNO3 = C(:,3); CO2 = C(:,4);
    %SNH4 = S(1); SNO2=S(2); SNO3 = S(3); L = S(7);
    %fa = F(:,1); fn = F(:,2); fi = F(:,3); 
    change = 1; maa= 0; l= 0; droptol = 10e-11;
    mmax=4; GNH4 =[]; GNO2=[]; GO2=[];
    while change > tol
        l = l+1;
        Gtemp = RDE(C,S,F);
        GNH4new = [0;D2\Gtemp(:,1)]+S(1);
        GNO2new = [0;D2\Gtemp(:,2)]+S(2);
        GO2new = [0;D2\Gtemp(:,3)]+SO2;
        Gnew = [GNH4new,GNO2new,GO2new];
        Fnew = Gnew-C; FNH4new=GNH4new-C(:,1);
        FNO2new=GNO2new-C(:,2); FO2new=GO2new-C(:,3);
        if norm(Fnew,inf) < tol
            w = Gnew;
            break;
        end
         if l > 1
             deltaFNH4 = FNH4new-FNH4old; deltaGNH4 = GNH4new-GNH4old;
             deltaFNO2 = FNO2new-FNO2old; deltaGNO2 = GNO2new-GNO2old;
             deltaFO2 = FO2new-FO2old; deltaGO2=GO2new-GO2old;
             if maa<mmax
                 GNH4 = [GNH4,deltaGNH4];
                 GNO2 = [GNO2,deltaGNO2];
                 GO2 = [GO2,deltaGO2];
             else
                 GNH4 = [GNH4(:,2:maa),deltaGNH4];
                 GNO2 = [GNO2(:,2:maa),deltaGNO2];
                 GO2 = [GO2(:,2:maa),deltaGO2];
             end
             maa=maa+1;
         end
         FNH4old=FNH4new; FNO2old=FNO2new;
         FO2old = FO2new; Fold = Fnew;
         GNH4old = GNH4new; GNO2old= GNO2new;
         GO2old = GO2new; Gold = Gnew;
         if maa == 0
             Cnew = Gnew;
         else
             if maa==1
                 QNH4(:,1) = deltaFNH4/norm(deltaFNH4);
                 RNH4(1,1) = norm(deltaFNH4);
                 QNO2(:,1) = deltaFNO2/norm(deltaFNO2);
                 RNO2(1,1) = norm(deltaFNO2);
                 QO2(:,1) = deltaFO2/norm(deltaFO2);
                 RO2(1,1) = norm(deltaFO2);
             else
                 if maa>mmax
                     [QNH4,RNH4] = qrdelete(QNH4,RNH4);
                     [QNO2,RNO2] = qrdelete(QNO2,RNO2);
                     [QO2,RO2] = qrdelete(QO2,RO2);
                     maa = maa-1;
                 end
                 for kk = 1:maa-1
                     RNH4(kk,maa) = QNH4(:,kk)'*deltaFNH4;
                     RNO2(kk,maa) = QNO2(:,kk)'*deltaFNO2;
                     RO2(kk,maa) = QO2(:,kk)'*deltaFO2;
                     deltaFNH4 = deltaFNH4-RNH4(kk,maa)*QNH4(:,kk);
                     deltaFNO2 = deltaFNO2-RNO2(kk,maa)*QNO2(:,kk);
                     deltaFO2 = deltaFO2-RO2(kk,maa)*QO2(:,kk);
                 end
                 QNH4(:,maa) = deltaFNH4/norm(deltaFNH4);
                 RNH4(maa,maa) = norm(deltaFNH4);
                 QNO2(:,maa) = deltaFNO2/norm(deltaFNO2);
                 RNO2(maa,maa) = norm(deltaFNO2);
                 QO2(:,maa) = deltaFO2/norm(deltaFO2);
                 RO2(maa,maa) = norm(deltaFO2);                    
             end
           if rcond(RNH4) > droptol
               gammaNH4 = RNH4\(QNH4'*FNH4new);
           else
               [UNH4,S1NH4,VNH4] = svd(GNH4,'econ');
               gammaNH4 = VNH4*inv(S1NH4)*UNH4'*FNH4new;
           end
           if rcond(RNO2) > droptol
           gammaNO2 = RNO2\(QNO2'*FNO2new);
           else
               [UNO2,S1NO2,VNO2] = svd(GNO2,'econ');
               gammaNO2 =VNO2*inv(S1NO2)*UNO2'*FNO2new;
           end
           if rcond(RO2) > droptol
           gammaO2 = RO2\(QO2'*FO2new);
           else
               [UO2,S1O2,VO2] = svd(GO2,'econ');
               gammaO2 = VO2*inv(S1O2)*UO2'*FO2new;
           end
           Cnew = [GNH4new-GNH4*gammaNH4,GNO2new-GNO2*gammaNO2,GO2new-GO2*gammaO2];          
         end
         change = norm(C-Cnew,inf);
         C = Cnew;
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
    %CNH4 = C(:,1); CNO2 = C(:,2); CO2 = C(:,3);
    %SNH4 = S(1); SNO2=S(2);
    %Xa = S(3); Xn = S(4); Xi = S(5); L = S(6);
    %fa = F(:,1); fn = F(:,2); fi = F(:,3);
    % Biofilm Velocity and Flux Equations
         JNH4a = (1/ya+ia)*mua(C(:,1),C(:,3));
         JNH4b = (ia-ii*fxi)*muaO(C(:,3));
         JNH4n = ia*mun(C(:,2),C(:,3));
         JNH4m = (ia-ii*fxi)*munO(C(:,3));
    JNH4 = (2/S(6))*int*((JNH4a-JNH4b).*F(:,1)*rho+(JNH4n-JNH4m).*F(:,2)*rho);
    JNO2 = (2/S(6))*int*(mun(C(:,2),C(:,3)).*F(:,2)*rho/yn-mua(C(:,1),C(:,3)).*F(:,1)*rho/ya);
    % Pieces of SNH4 Equation
         RNH4a =(1/ya+ia)*mua(S(1),SO2); 
         RNH4b =(ia-ii*fxi)*muaO(SO2);
         RNH4n =ia*mun(S(2),SO2);
         RNH4m =(ia-ii*fxi)*munO(SO2);
    w = [d*(SNH4in-S(1))-1/V*((RNH4a-RNH4b)*S(3)+(RNH4n-RNH4m)*S(4))-1/V*A*JNH4(1);
         d*(SNO2in-S(2))-1/V*(1/yn*mun(S(2),SO2)*S(4)-1/ya*mua(S(1),SO2)*S(3))-1/V*A*JNO2(1);
         S(3)*Rxa(S(1),SO2)-(d+alpha)*S(3)+A*rho*F(1,1)*E*S(6)^2;
         S(4)*Rxn(S(2),SO2)-(d+alpha)*S(4)+A*rho*F(1,2)*E*S(6)^2;
         mui(S(3),S(4),SO2)-S(5)*(d+alpha)+A*rho*F(1,3)*E*S(6)^2;
         Lp(C,S,F)];
    end
%% Initial Conditions
% Planktonic
  SNO2 = SNO2in; SNH4 = SNH4in ;   % Substrates
  SO2 = SO2in ; %SNO3 = SNO3in;
  Xa = 2E-6;  Xn = 2E-6; Xi = 0;  % Biomass
  L = 1E-6;
  S = [SNH4;SNO2;Xa;Xn;Xi;L];
% Biofilm
  CO2 = SO2*ones(size(x));   CNO2 = SNO2*ones(size(x));
  CNH4 = SNH4*ones(size(x));%  CNO3 = SNO3*ones(size(x));
  C = [CNH4,CNO2,CO2];
  fa = 0.5*ones(size(x)); fn = 0.5*ones(size(x)); fi = zeros(size(x));
  F = [fa,fn,fi];
  C = biofilmbvp(C,S,F);
%% Set storage, preallocate size for speed
SS = zeros(6,points); SS(:,1) = S(:);
tt = zeros(1,points);
fafa = zeros(N+1,points); fafa(:,1)=fa;
fnfn = zeros(N+1,points); fnfn(:,1)=fn;
fifi = zeros(N+1,points); fifi(:,1)=fi;
CNH4C = zeros(N+1,points); CNH4C(:,1) = C(:,1);
CNO2C = zeros(N+1,points); CNO2C(:,1) = C(:,2);
%CNO3C = zeros(N+1,points); CNO3C(:,1) = CNO3';
CO2C = zeros(N+1,points); CO2C(:,1) = C(:,3);
reverseStr = '';

%% time stepping function
% Using the exponential timestepping method outlined in KassamTrefethen (2005);
    function [S,F] = timestep(C,S,F)
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
        S = S+FS1*NSu+2*FS2*(NSa+NSb)+FS3*NSc;
        F = F+fU1*Nu+2*fU2*(Na+Nb)+fU3*Nc;
    end
        
        
%% Main Timestepping Loop
  for i = 1:nplots
    % Save into temporary Storage
    t= t+dt; 
    [S,F] = timestep(C,S,F);
    C = biofilmbvp(C,S,F);
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
