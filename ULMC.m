%% This code is used to draw ULMC.eps
%% Use RCADULMCFError and RCDULMCFError to draw graph
%% parameter
clear
d=1000;%dimension
N=100000;%particle number
Dt = 1e-3;% time step
RCDULMCError=zeros(11,10000);%error for RCD with diffrent steps and different time steps
RCADULMCError=zeros(11,10000);%error for RCAD with diffrent steps and different time steps
RCDULMCFError=zeros(11);%equilibrium error for RCD with different time steps
RCADULMCFError=zeros(11);%equilibrium error for RCAD with different time steps
q=0;
for p=[1 1.5 2 2.5 3 3.5 4 5 6 7 8 10 15]
%% initial condition
samples=randn(N,2)+0.5*ones(N,2);
RCDULMCsamples=samples;
RCADULMCsamples=samples;
RCDULMCW2=norm(samples,'fro').^2/N-1;
RCADULMCW2=norm(samples,'fro').^2/N-1;
dt=Dt/p;
sqrdt = dt^(1/2);
Cov=[dt-3/4-1/4*exp(-4*dt)+exp(-2*dt) 0.5*(1+exp(-4*dt)-2*exp(-2*dt));0.5*(1+exp(-4*dt)-2*exp(-2*dt)) 1-exp(-4*dt)]^0.5;
%% first flux
RCADULMCg=samples(:,1);
m=0;
while m<(2*10^4*p)
W=randn(N,1);K=randn(2,N);choice=(floor(rand(N,1)*d)==zeros(N,1));
%% RCD-ULMC
RCDULMCflux = RCDULMCsamples(:,1).*choice;
ExRCDULMCsamples = RCDULMCsamples(:,1)+0.5*(1-exp(-2*dt))*RCDULMCsamples(:,2)-1/2*(dt-0.5*(1-exp(-2*dt)))*d*RCDULMCflux;
EvRCDULMCsamples = exp(-2*dt)*RCDULMCsamples(:,2)-0.5*(1-exp(-2*dt))*d*RCDULMCflux;
RCDULMCsamples = (Cov*K+[ExRCDULMCsamples';EvRCDULMCsamples'])';
%% RCAD-ULMC
flux=(1-d)*RCADULMCg;
RCADULMCg = RCADULMCg.*(1-choice)+RCADULMCsamples(:,1).*choice;
RCADULMCflux=flux+d*RCADULMCg;
ExRCADULMCsamples = RCADULMCsamples(:,1)+0.5*(1-exp(-2*dt))*RCADULMCsamples(:,2)-1/2*(dt-0.5*(1-exp(-2*dt)))*RCADULMCflux;
EvRCADULMCsamples = exp(-2*dt)*RCADULMCsamples(:,2)-0.5*(1-exp(-2*dt))*RCADULMCflux;
RCADULMCsamples = (Cov*K+[ExRCADULMCsamples';EvRCADULMCsamples'])';
m=m+1;
if (floor(m/10)==m/10)
         m
         p
    RCADULMCW2=[RCADULMCW2;abs(norm(RCADULMCsamples(:,1),'fro').^2/N-1)];
    RCDULMCW2=[RCDULMCW2;abs(norm(RCDULMCsamples(:,1),'fro').^2/N-1)];     
end
end
q=q+1;
RCADULMCError(q,1:length(RCADULMCW2))=RCADULMCW2';
RCDULMCError(q,1:length(RCADULMCW2))=RCDULMCW2';
RCADULMCFError(q)=RCADULMCW2(length(RCADULMCW2));
RCDULMCFError(q)=RCDULMCW2(length(RCDULMCW2));
save ULMCError
end





