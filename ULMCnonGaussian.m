%% This code is used to draw ULMCnonGaussian.eps
%% Use RCADULMCFError and RCDULMCFError to draw graph
%% parameter
clear
d=1000;%dimension
N=1000000;%particle number
Dt = 1e-3;% time step
RCDULMCError=zeros(11,10000);%error for RCD with diffrent steps and different time steps
RCADULMCError=zeros(11,10000);%error for RCAD with diffrent steps and different time steps
RCDULMCFError=zeros(11);%equilibrium error for RCD with different time steps
RCADULMCFError=zeros(11);%equilibrium error for RCAD with different time steps
q=0;
samples=randn(N,2)+0.5*ones(N,2);
RCDULMCsamples=samples;
RCADULMCsamples=samples;
abs(sum(samples(:,1).^2)/N-5)
for p=[1 1.2 1.6 2 2.5 3 4 5 6 7 9 10]
%% initial condition
RCDULMCW2=abs(sum((samples(:,1)).^2)/N-5);
RCADULMCW2=abs(sum((samples(:,1)).^2)/N-5);
dt=Dt/p;
sqrdt = dt^(1/2);
Cov=[dt-3/4-1/4*exp(-4*dt)+exp(-2*dt) 0.5*(1+exp(-4*dt)-2*exp(-2*dt));0.5*(1+exp(-4*dt)-2*exp(-2*dt)) 1-exp(-4*dt)]^0.5;
x=RCADULMCsamples(:,1);
RCADULMCgini=(x-2).*exp(-(x-2).^2/2)+(x+2).*exp(-(x+2).^2/2);%%Derivative calculation
RCADULMCgini=RCADULMCgini./(exp(-(x-2).^2/2)+exp(-(x+2).^2/2));%%Derivative calculation   
%% first flux
RCADULMCg=RCADULMCgini;
m=0;
while m<(10^4*p)
K=randn(2,N);choice=(floor(rand(N,1)*d)==zeros(N,1));
%% RCD-ULMC
RCDULMCflux = zeros(N,1);
x=RCDULMCsamples(:,1);
RCDULMCflux=(x-2).*exp(-(x-2).^2/2)+(x+2).*exp(-(x+2).^2/2);%%Derivative calculation
RCDULMCflux=RCDULMCflux./(exp(-(x-2).^2/2)+exp(-(x+2).^2/2));%%Derivative calculation
RCDULMCflux = RCDULMCflux.*choice;
ExRCDULMCsamples = RCDULMCsamples(:,1)+0.5*(1-exp(-2*dt)).*RCDULMCsamples(:,2)-1/2*(dt-0.5*(1-exp(-2*dt)))*d.*RCDULMCflux;
EvRCDULMCsamples = exp(-2*dt).*RCDULMCsamples(:,2)-0.5*(1-exp(-2*dt))*d.*RCDULMCflux;
RCDULMCsamples = (Cov*K+[ExRCDULMCsamples';EvRCDULMCsamples'])';
%% RCAD-ULMC
flux=(1-d)*RCADULMCg;
x=RCADULMCsamples(:,1);
RCADULMCflux = (x-2).*exp(-(x-2).^2/2)+(x+2).*exp(-(x+2).^2/2);%%Derivative calculation
RCADULMCflux = RCADULMCflux./(exp(-(x-2).^2/2)+exp(-(x+2).^2/2));%%Derivative calculation
RCADULMCflux = RCADULMCflux.*choice;
RCADULMCg = RCADULMCg.*(1-choice)+RCADULMCflux.*choice;
RCADULMCflux=flux+d*RCADULMCg;
ExRCADULMCsamples = RCADULMCsamples(:,1)+0.5*(1-exp(-2*dt)).*RCADULMCsamples(:,2)-1/2*(dt-0.5*(1-exp(-2*dt))).*RCADULMCflux;
EvRCADULMCsamples = exp(-2*dt).*RCADULMCsamples(:,2)-0.5*(1-exp(-2*dt)).*RCADULMCflux;
RCADULMCsamples = (Cov*K+[ExRCADULMCsamples';EvRCADULMCsamples'])';
m=m+1;
if (floor(m/10)==m/10)
       RCADULMCW2=[RCADULMCW2;abs(sum((RCADULMCsamples(:,1)).^2)/N-5)];
       RCDULMCW2=[RCDULMCW2;abs(sum((RCDULMCsamples(:,1)).^2)/N-5)];   
       m
       p
end
end
q=q+1;
RCADULMCError(q,1:length(RCADULMCW2))=RCADULMCW2';
RCDULMCError(q,1:length(RCDULMCW2))=RCDULMCW2';
RCADULMCFError(q)=RCADULMCW2((length(RCADULMCW2)));
RCDULMCFError(q)=RCDULMCW2((length(RCDULMCW2)));
save ULMCNonGussianError
end





