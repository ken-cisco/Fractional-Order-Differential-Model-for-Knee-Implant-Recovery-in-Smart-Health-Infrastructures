% This Numerical Scheme is for Caputo-Fabrizio-Caputo (Exponential Law) Diferential Operator

clc;clear;close all;
%Inputs
h=0.002; t(1)=0; tfinal=1000; t=t(1):h:tfinal; N=ceil((tfinal-t(1))/h);

%Initial Condition
R(1)=0.5;I(1)=0.7;M(1)=0.6;

%Constant fractional order
alpha=0.976;
tau=5

%parameters of the model
a1=0.5;a2=0.3;a3=0.0001;a4=0.2;a5=0.05;a6=0.002;a7=0.6;a8=0.04;a10=0.2;

%Given below is set of ODEs (Knee Implant Rcovery delay Model Control)
f1=@(t,R,I,M) a1+a2.*M-a3.*I.*(t-tau).*R-a4.*R;
f2=@(t,R,I,M) a5.*I-a3.*I.*(t-tau).*R-a6.*M;
f3=@(t,R,I,M) a7-a8.*I-a10.*M;

R(2)=R(1)+((h.^alpha)./(gamma(alpha+1))).*f1(t(1),R(1),I(1),M(1));
I(2)=I(1)+((h.^alpha)./(gamma(alpha+1))).*f2(t(1),R(1),I(1),M(1));
M(2)=M(1)+((h.^alpha)./(gamma(alpha+1))).*f3(t(1),R(1),I(1),M(1));

%Constant version of the Caputo-Fabrizio-Caputo Algorithm starts
for n=2:N
R(n+1)=R(n)+(0.5.*(2-alpha).*(1-alpha)+0.25.*3.*h.*alpha.*(2-alpha)).*f1(t(n),R(n),I(n),M(n))-(0.5.*(2-alpha).*(1-alpha)+0.25.*h.*alpha.*(2-alpha)).*f1(t(n-1),R(n-1),I(n-1),M(n-1));    
I(n+1)=I(n)+(0.5.*(2-alpha).*(1-alpha)+0.25.*3.*h.*alpha.*(2-alpha)).*f2(t(n),R(n),I(n),M(n))-(0.5.*(2-alpha).*(1-alpha)+0.25.*h.*alpha.*(2-alpha)).*f2(t(n-1),R(n-1),I(n-1),M(n-1));
M(n+1)=M(n)+(0.5.*(2-alpha).*(1-alpha)+0.25.*3.*h.*alpha.*(2-alpha)).*f3(t(n),R(n),I(n),M(n))-(0.5.*(2-alpha).*(1-alpha)+0.25.*h.*alpha.*(2-alpha)).*f3(t(n-1),R(n-1),I(n-1),M(n-1));

t(n+1)=t(n)+h;
end
plot3(R,I,M, 'b-','LineWidth', 5);
xlabel('Recovery Dynamics');
ylabel('Inflammation Dynamics');
zlabel('Mechanical Function')
grid on;
title('Knee Implant Recovery Dynamics and Mechanical Function in the Absence of Delay');