
%% Impulsive system: MBBR with SBR cycles - Partial Nitritation/Anammox
% Wade, Meadows & Wolkowicz (2019)
% 0D model by Laureni et al. (2019). Water Research.
% MatContM implementation by Ting-Hao Hsu (Uni. Miami)


%% Preparation
run('./matcontm5p4/init');
global opt cds fpmds
opt = contset();

%% Set Numerics
sys = @pna_sysFinal_NS;
ap = 16; % active parameter

%% Set Parameters
m1 = 0.3;
y1 = 0.18;
ka1 = 2.4;
ko1 = 0.6;
in1 = 0.083;
m2 = 0.34;
y2 = 0.08;
kn2 = 0.5;
ko2 = 0.4;
in2 = 0.083;
%p3 = 6.88;    %Varied in paper [0-24 mgCOD/L.d]
%r3 = 150;    %Varied in paper [0-300 mgN/L.d]
y3 = 0.17;
ka3 = 0.03;
kn3 = 0.005;
in3 = 0.058;

aIN = 20;
oIN = 0.3;  %Varied in paper [0.15-1.5 mg/L]
s1bar = 2;
r1 = 0.005;  %Varied in paper [0.4-1.7%]

%AMX = p3/m3;
R3 = 86;

%d1 = 24/5.9; %HRT (1/h)
%d2 = 1/24.5; %SRT (1/d)

r = .5;

par = [m1,ka1,ko1,y1,in1,m2,kn2,ko2,y2,in2,ka3,kn3,y3,in3,aIN,oIN,R3,s1bar,r,r1].';

%% Set Functions
sp = r*aIN + (1-r)*s1bar;

xini = [0.2717899; 215.5964; 40.76058];

i = 1; % i iteration
[u0,v0] = init_FPm_FPm(sys,xini,par,ap,i);

%% Follow Per1
opt = contset;
opt = contset(opt,'MaxNumPoints',200);
opt = contset(opt,'Singularities',1);
opt = contset(opt,'Multipliers',1);
opt = contset(opt,'backward',0);
opt = contset(opt,'Increment',1e-2);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters', 10);
opt = contset(opt, 'MaxTestIters', 10);
opt = contset(opt, 'FunTolerance', 1e-6);
opt = contset(opt, 'VarTolerance', 1e-6);
opt = contset(opt, 'TestTolerance', 1e-6);
opt = contset(opt, 'InitStepSize', 0.01);

[u1,v1,s1,~,~] = cont(@fixedpointmap,u0,v0,opt);
opt = contset(opt,'backward',1);

u0n = u0;
v0n = v0;
niter = 1;
s2 = 0;
while length(s2)<=2 && niter<20
    
    [u2,v2,s2,~,~] = cont(@fixedpointmap,u0n,v0n,opt);
    [valc,dirc] = min([u2(4,1),u2(4,end)]);
    if dirc == 2
        indx = length(u2(4,:));
    else
        indx = 1;
    end
    u0n = u2(:,indx);
    v0n = v2(:,indx);
    niter = niter+1;
    
end

%Print eigenvalues
for k = 2:length(s1)-1
    fprintf('\n')
    fprintf(['Eigenvalues for ',s1(k).label,'\n'])
    fprintf('%f \n',s1(k).data.eval')
end

for k = 2:length(s2)-1
    fprintf('\n')
    fprintf(['Eigenvalues for ',s2(k).label,'\n'])
    fprintf('%f \n',s2(k).data.eval')
end

% Now continue from bifurcation (BP)
% Find BP index on curve

for ii = 1:length(s2)
    if strmatch(s2(ii).label,'BP  ');
        bpindex = ii;
    end
end

ap2 = [16 17]; %Two parameter [o2, Ramx]

u2i = u2(1:3,s2(bpindex).index); %Starting point
p1 = par; p1(fpmds.ActiveParams) = u2(4,s2(bpindex).index);

[ubp,vbp] = init_LPm_LPm(sys,u2i,p1,ap2,2);
opt = contset(opt,'backward',0);
opt = contset(opt,'singularities',0);
opt = contset(opt,'Increment',1e-2);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters', 10);
opt = contset(opt, 'MaxTestIters', 10);
opt = contset(opt, 'FunTolerance', 1e-6);
opt = contset(opt, 'VarTolerance', 1e-6);
opt = contset(opt, 'TestTolerance', 1e-5);
opt = contset(opt, 'InitStepSize', 0.01);
opt = contset(opt, 'Adapt', 3);
opt = contset(opt,'MaxNumPoints',5);
[xbp1,vbp1,sbp1,hbp1,fbp1]=cont(@limitpointmap,ubp,vbp,opt);

opt = contset(opt,'backward',1);
[xbp2,vbp2,sbp2,hbp2,fbp2]=cont(@limitpointmap,ubp,vbp,opt);

figure(1)
cpl(xbp1,vbp1,sbp1,[4 5])
hold on
cpl(xbp2,vbp2,sbp2,[4 5])
xlabel('O_2')
ylabel('r_{AMX}')
savefig('Par_o2_rAMX.fig')


ap3 = [16 20]; %Two parameter [O2, f_was]
[ubpa,vbpa] = init_LPm_LPm(sys,u2i,p1,ap3,2);

opt = contset(opt,'backward',0);
[xbp3,vbp3,sbp3,hbp3,fbp3]=cont(@limitpointmap,ubpa,vbpa,opt);

opt = contset(opt,'backward',1);
[xbp4,vbp4,sbp4,hbp4,fbp4]=cont(@limitpointmap,ubpa,vbpa,opt);

figure(2)
cpl(xbp3,vbp3,sbp3,[4 5])
hold on
cpl(xbp4,vbp4,sbp4,[4 5])
xlabel('O_2')
ylabel('f_{WAS}')
savefig('Par_o2_fwas.fig')

ap4 = [17 20]; %Two parameter [RAMX, f_was]
[ubpb,vbpb] = init_LPm_LPm(sys,u2i,p1,ap4,2);

opt = contset(opt,'backward',0);
[xbp5,vbp5,sbp5,hbp5,fbp5]=cont(@limitpointmap,ubpb,vbpb,opt);

opt = contset(opt,'backward',1);
[xbp6,vbp6,sbp6,hbp6,fbp6]=cont(@limitpointmap,ubpb,vbpb,opt);

figure(3)
cpl(xbp5,vbp5,sbp5,[4 5])
hold on
cpl(xbp6,vbp6,sbp6,[4 5])
xlabel('r_{AMX}')
ylabel('f_{WAS}')
savefig('Par_rAMX_fwas.fig')

save(['Data_Run_',date])

return

%Now continue from bifurcation (BP)
% ap2 = [16 17]; %Two parameter
% u1i = u1(1:3,s1(2).index); %Starting point
% p1 = par; p1(fpmds.ActiveParams) = u1(4,s1(2).index);
% 
% [ubp,vbp] = init_BPm_FPm(sys,u1i,p1,s1(2),0.01,1);



%% Initialize Per2 - THIS IS FOR E(X1,X2) = [1,O]
uini = [3.061916; 9.92211; 0];
ap=17;
[u0,v0] = init_FPm_FPm(sys,uini,par,ap,i);

%% Follow Per2
opt = contset;
opt = contset(opt,'MaxNumPoints',100);
opt = contset(opt,'Singularities',1);
opt = contset(opt,'Multipliers',1);
opt = contset(opt,'backward',0);
opt = contset(opt,'Increment',1e-3);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters', 10);
opt = contset(opt, 'MaxTestIters', 10);
opt = contset(opt, 'FunTolerance', 1e-7);
opt = contset(opt, 'VarTolerance', 1e-7);
opt = contset(opt, 'TestTolerance', 1e-7);
opt = contset(opt, 'InitStepSize', 0.01);
[u3,v3,s3,~,~] = cont(@fixedpointmap,u0,v0,opt);

%% Follow Backward Per2
opt = contset(opt,'backward',1);
[u4,v4,s4,~,~] = cont(@fixedpointmap,u0,v0,opt);

%% Plot Per2
xx3 = u3(2,:);
yy3 = u3(3,:);
pp3 = u3(end,:);
ss3 = s3;
plot3(xx3,yy3,pp3,'k-','LineWidth',2);
numLable = size(ss3,1) - 2;
for js=2:numLable+1
    ind = ss3(js).index;
    plot3(xx3(ind),yy3(ind),pp3(ind),'r*','MarkerSize',6);
    text(xx3(ind),yy3(ind),pp3(ind),ss3(js).label);
end

xx4 = u4(2,:);
yy4 = u4(3,:);
pp4 = u4(end,:);
ss4 = s4;
plot3(xx4,yy4,pp4,'k-','LineWidth',2);
numLable = size(ss4,1) - 2;
for js=2:numLable+1
    ind = ss4(js).index;
    plot3(xx4(ind),yy4(ind),pp4(ind),'r*','MarkerSize',6);
    text(xx4(ind),yy4(ind),pp4(ind),ss4(js).label);
end



% cpl(u1,v1,s1,[4 2]);
% xlabel('m_1');
% ylabel('s');

% %% Follow Branch Point
% opt = contset;
% opt = contset(opt,'Multipliers',1);
% opt = contset(opt,'Singularities',1);
% 
% % bpind = s1(2).index;
% indbps = find(arrayfun(@(c)strcmp(c.label,'BP  '),s1),1);
% indbp = s1(indbps).index;
% xini = u1(1:end-1,indbp);
% 
% par(fpmds.ActiveParams) = u1(end,indbp);
% 
% % Go Forward
% opt = contset(opt,'backward',0);
% opt = contset(opt,'MaxNumPoints',50);
% 
% h = 0.01; % initial amplitute
% i = 1; % iteration number of the branch
% [u2,v2] = init_BPm_FPm(sys,xini,par,s1(indbps),h,i);
% 
% [u21,v21,s21,h21,f21] = cont(@fixedpointmap,u2,v2,opt);



% %% Go Backward
% opt = contset(opt,'backward',1);
% [u22,v22,s22,h22,f22] = cont(@fixedpointmap,u2,v2,opt);

return
%% Save Figure
h = fig1;
filename = 'fig_cfl_competition';
% set(gca,'XColor','none','YColor','none');
% grid off;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,filename,'-dpdf','-r0')
% set(gca,'XColor','k','YColor','k');
% grid on;