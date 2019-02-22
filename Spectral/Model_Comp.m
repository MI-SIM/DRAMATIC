%% Model_Comp.m

N = 16;
T = 5000;

[t,x,S,C,f1,f2] = competition(T,N);

figure(1);  clf;
plot(t,S(1,:),'LineWidth',2);
legend('S');

figure(2); clf
plot(t,S(2:3,:),'LineWidth',2);
legend('u_1','u_2');