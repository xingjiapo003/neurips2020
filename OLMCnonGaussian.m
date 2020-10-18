%% This code is used to draw OLMCnonGaussian.eps
%% Use RCADOLMCFError and RCDOLMCFError to draw graph
%% parameter
clear
d=1000;%dimension
N=1000000;%particle number
Dt = 1e-3;% time step
RCDOLMCError=zeros(20,20000);%error for RCD with diffrent steps and different time steps
RCADOLMCError=zeros(20,20000);%error for RCAD with diffrent steps and different time steps
RCDOLMCFError=zeros(20);%equilibrium error for RCD with different time steps
RCADOLMCFError=zeros(20);%equilibrium error for RCAD with different time steps
samples=randn(N,1)+0.5*ones(N,1);
abs(sum((samples-sum(samples)/N).^2)/N-5)
q=0;
for p=[3 3.25 3.5 3.75 4 4.5 5 6 7 8 10]
%% initial condition
RCDOLMCsamples=samples;
RCADOLMCsamples=samples;
RCDOLMCW2=abs(sum((samples-sum(samples)/N).^2)/N-5);
RCADOLMCW2=abs(sum((samples-sum(samples)/N).^2)/N-5);
dt=Dt/p;
sqrdt = dt^(1/2);
x=RCADOLMCsamples;
RCADOLMCgini=(x-2).*exp(-(x-2).^2/2)+(x+2).*exp(-(x+2).^2/2);%%Derivative calculation
RCADOLMCgini=RCADOLMCgini./(exp(-(x-2).^2/2)+exp(-(x+2).^2/2));%%Derivative calculation   
%% first flux
RCADOLMCg=RCADOLMCgini;
m=0;
while m<(10^5*p)
W=randn(N,1);choice=(floor(rand(N,1)*d)==zeros(N,1));
%% RCD-OLMC
x=RCDOLMCsamples;
RCDOLMCflux=(x-2).*exp(-(x-2).^2/2)+(x+2).*exp(-(x+2).^2/2);%%Derivative calculation
RCDOLMCflux=RCDOLMCflux./(exp(-(x-2).^2/2)+exp(-(x+2).^2/2));%%Derivative calculation
RCDOLMCflux = RCDOLMCflux.*choice;
RCDOLMCsamples = RCDOLMCsamples - d*RCDOLMCflux * dt + sqrt(2)*sqrdt*W;
%% RCAD-OLMC
flux=(1-d)*RCADOLMCg;
x=RCADOLMCsamples;
RCADOLMCflux = (x-2).*exp(-(x-2).^2/2)+(x+2).*exp(-(x+2).^2/2);%%Derivative calculation
RCADOLMCflux = RCADOLMCflux./(exp(-(x-2).^2/2)+exp(-(x+2).^2/2));%%Derivative calculation
RCADOLMCflux = RCADOLMCflux.*choice;
RCADOLMCg = RCADOLMCg.*(1-choice)+RCADOLMCflux.*choice;
flux=flux+d*RCADOLMCg;
RCADOLMCsamples = RCADOLMCsamples - flux * dt + sqrt(2)*sqrdt*W;
m=m+1;
if (floor(m/10)==m/10)     
           m
           p    
     RCADOLMCW2=[RCADOLMCW2;abs(sum((RCADOLMCsamples(:,1)).^2)/N-5)];
     RCDOLMCW2=[RCDOLMCW2;abs(sum((RCDOLMCsamples(:,1)).^2)/N-5)];     
end
end
q=q+1;
RCADULMCError(q,1:length(RCADULMCW2))=RCADULMCW2';
RCDULMCError(q,1:length(RCADULMCW2))=RCDULMCW2';
RCADULMCFError(q)=RCADULMCW2(length(RCADULMCW2));
RCDULMCFError(q)=RCDULMCW2(length(RCDULMCW2));
save OLMCnonGaussianError
end





