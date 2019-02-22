function [Sl,Sr] = intmat(N,a,b)

% Spectral integration matrices of order N+1;
j=[0:N/2]';

% Form Ut: the upper N/2+1 rows of U
Ut(:,1:2) = [cos(j*pi/N),1/2*cos(j*pi/N).^2];
Ut(:,3:N+1) = (cos(j*[3:N+1]*pi/N)./(ones(N/2+1,1)*[3:N+1])...
    -cos(j*[1:N-1]*pi/N)./(ones(N/2+1,1)*[1:N-1]))/2;
% Discrete Chebyshev Transform of Ut';
B = real(fft([Ut,Ut(:,N:-1:2)]'))';
At = [B(1:N/2+1,1)/2,B(1:N/2+1,2:N),B(1:N/2+1,N+1)/2]/N;

% Construct A from the ACS Property
A = [At;-rot90(At(1:N/2,:),2)];

% Projection and centro-transpose
Sl = A -ones(N+1,1)*A(N+1,:);  
Sr = rot90(Sl,2); Sr = Sr*(a-b)/2; Sl =Sl*(b-a)/2;