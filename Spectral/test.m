% p30.m - spectral integration (Clenshaw Curtis)
% Computation: Various values of N, four functions

Nmax=50; E = zeros(4,Nmax); clf'
for N = 1:Nmax; i = 1:N;
    [x,w] = clencurt(N);
    f = abs(x).^3;     E(1,N) = abs(w*f-0.5);
    f = exp(-x.^(-2));  E(2,N) = abs(w*f-2*(exp(-1)+sqrt(pi)*(erf(1)-1)));
    f = 1./(1+x.^2);    E(3,N) = abs(w*f - pi/2);
    f = x.^10;          E(4,N) = abs(w*f - 2/11);
end

%plot results
labels = {'|x|^3','exp(-x^{-2})','1/(1+x^2)','x^{10}'};
for iplot = 1:4,
    subplot(3,2,iplot)
    semilogy(E(iplot,:)+1e-100,'.','markersize',12), hold on
    plot(E(iplot,:)+1e-100,'linewidth',0.8)
    axis([0 Nmax 1e-18 1e3]), grid on
    set(gca,'xtick',0:10:Nmax,'ytick',(10).^(-15:5:0))
    ylabel error, text(32,.004,labels(iplot))
end
