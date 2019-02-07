
T =90; %Maximum time
N = 16;%Number of Chebyshev points

[t,S,CNH4,CNO2,CO2,fa,fn,fi,x] = PartialNitritation_NO3(T,N);

figure(1);
clf;
hold on
plot(t,S(1:2,:),'LineWidth',2)
legend('NH_4','NO_2');
figure(2);
clf;
hold on
plot(t,S(3:5,:),'LineWidth',2)
legend('AOB','NOB','IB');
