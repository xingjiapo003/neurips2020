%% This code is used to draw OLMC.eps
%% Use RCADOLMCFError and RCDOLMCFError to draw graph
%% parameter
clear
d=1000;%dimension
N=1000000;%particle number
Dt = 5*1e-4;% time step
RCDOLMCError=zeros(20,20000);%error for RCD with diffrent steps and different time steps
RCADOLMCError=zeros(20,20000);%error for RCAD with diffrent steps and different time steps
RCDOLMCFError=zeros(20);%equilibrium error for RCD with different time steps
RCADOLMCFError=zeros(20);%equilibrium error for RCAD with different time steps
q=0;
for p=[1 1.5 2 2.5 3 4 5 6 8 9 10 15]
%% initial condition
samples=randn(N,1)+0.5*ones(N,1);
RCDOLMCsamples=samples;
RCADOLMCsamples=samples;
RCDOLMCW2=norm(samples,'fro').^2/N-1;
RCADOLMCW2=norm(samples,'fro').^2/N-1;
dt=Dt/p;
sqrdt = dt^(1/2);
%% first flux
RCADOLMCg=samples;
m=0;
while m<(2*10^4*p)
W=randn(N,1);choice=(floor(rand(N,1)*d)==zeros(N,1));
%% RCD-OLMC
RCDOLMCflux = RCDOLMCsamples.*choice;
RCDOLMCsamples = RCDOLMCsamples - d*RCDOLMCflux * dt + sqrt(2)*sqrdt*W;
%% RCAD-OLMC
flux=(1-d)*RCADOLMCg;
RCADOLMCg = RCADOLMCg.*(1-choice)+RCADOLMCsamples.*choice;
flux=flux+d*RCADOLMCg;
RCADOLMCsamples = RCADOLMCsamples - flux * dt + sqrt(2)*sqrdt*W;
m=m+1;
if (floor(m/10)==m/10)
         m
         p
    RCADOLMCW2=[RCADOLMCW2;abs(norm(RCADOLMCsamples,'fro').^2/N-1)];
    RCDOLMCW2=[RCDOLMCW2;abs(norm(RCDOLMCsamples,'fro').^2/N-1)];     
end
end
q=q+1;
RCADOLMCError(q,1:length(RCADOLMCW2))=RCADOLMCW2';
RCDOLMCError(q,1:length(RCADOLMCW2))=RCDOLMCW2';
RCADOLMCFError(q)=RCADOLMCW2(length(RCADOLMCW2));
RCDOLMCFError(q)=RCDOLMCW2(length(RCDOLMCW2));
save OLMCError
end





