%Prepared by Samyak Shah for EE227 Assignment 1 Question 6
%Date: September 08,2019
%______________________________________________
%Defining the constants given in the question
%Energy is taken in eV
q=1.602e-19;
m=0.25*9.11e-31;%Effective mass
hcut=1.055e-34;
eps0=8.854E-12;
epsr=4;
I0=q*q/hcut;
%System Dimensions and properties
W=1e-6;%W=Width
L=10e-9;%Length of active region
t=1.5e-9;%Oxide thickness
Cgate=epsr*eps0*W*L/t;
Csourc=0.05*Cgate;
Cdrain=0.05*Cgate;
Ceq=Cgate+Csourc+Cdrain;
U0=q/Ceq;
alphag=Cgate/Ceq;
alphad=Cdrain/Ceq;
kT=0.025; %Room temperature is considered
mu=0;
ep=0.2;
v=1e5; %Escape velocity
g1=hcut*v/(q*L);
g2=g1;
g=g1+g2;
%Energy grid
NE=501;
E=linspace(-1,1,NE);
dE=E(2)-E(1);
D0=m*q*W*L/(pi*hcut*hcut);%Density of states per eV
D=D0*[zeros(1,251) ones(1,250)];
%D=(2*g/(2*pi))./((E.ˆ2)+((g/2)ˆ2));% Lorentzian Density of states per eV
%D=D./(dE*sum(D));%Normalizing to one
%ReferenCeq number of electrons
f0=1./(1+exp((E+ep-mu)./kT));N0=2*dE*sum(D.*f0);ns=N0/(L*W*1e4);%/cmˆ2
%Bias
IV=61;VV=linspace(0,0.6,IV);
for iV=1:IV
%Vg=0.25;Vd=VV(iV);
Vd=0.5;Vg=VV(iV);
mu1=mu;mu2=mu1-Vd;UL=-(alphag*Vg)-(alphad*Vd);
U=0;%Self-consistent field
dU=1;
while dU>1e-6
f1=1./(1+exp((E+UL+U+ep-mu1)./kT));
f2=1./(1+exp((E+UL+U+ep-mu2)./kT));
N(iV)=dE*sum(D.*((f1.*g1/g)+(f2.*g2/g)));
Unew=U0*(N(iV)-N0);dU=abs(U-Unew);
U=U+0.1*(Unew-U);
end
I(iV)=dE*I0*(sum(D.*(f1-f2)))*g1*g2/g;
end
hold on
h=plot(VV,I,'b');

xlabel(' Voltage (V) --->')
ylabel('Current (A) ---> ')
grid on

Vrn=zeros(1,501);
for i=252:501
    Vrn(i)=sqrt(2*m*E(i)*1.6e-19);
end
g1n=hcut*Vrn/(q*L);
g2n=g1n;
gn=g1n+g2n;

for iV=1:IV
Vg=0.25;Vd=VV(iV);
%Vd=0.5;Vg=VV(iV);
mu1=mu;mu2=mu1-Vd;UL=-(alphag*Vg)-(alphad*Vd);
U=0;%Self-consistent field
dU=1;
while dU>1e-6
f1=1./(1+exp((E+UL+U+ep-mu1)./kT));
f2=1./(1+exp((E+UL+U+ep-mu2)./kT));
N(iV)=dE*sum(D.*((f1.*0.5)+(f2.*0.5)));
Unew=U0*(N(iV)-N0);dU=abs(U-Unew);
U=U+0.1*(Unew-U);
end
I(iV)=dE*I0*(sum(D.*(f1-f2).*g1n*0.5));
end
figure
hold on
h=plot(VV,I,'b');

xlabel(' Gate Voltage (in V)')
ylabel('Drain Current (in A) ')
grid on