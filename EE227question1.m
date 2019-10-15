%Prepeared by Samyak Shah for EE227 Assignment 1 Question 1
%Date: September 08,2019
%______________________________________________
%Plotting the wavefunction corresponding to the first allowed energy state
x=linspace(-2e-9,6e-9,1000);
y1=zeros(1,1000);
v0=zeros(1,1000);
k1=7.27e8;
alpha1=6.196e9;
%Defining the values of Piecewise wavefunction coefficients
A1=2315.58;
C1=19735.49;
D1=A1;
F1=1.343e14;
for i=1:1000
    if x(i)<=0
        y1(i)=A1*exp(alpha1*x(i));
        v0(i)=2.8e4;
    elseif x(i)<=4e-9
        y1(i)=C1*sin(k1*x(i))+D1*cos(k1*x(i));
        v0(i)=0;
    elseif x(i)>4e-9
        y1(i)=F1*exp(-alpha1*x(i));
        v0(i)=2.8e4;
    end
end
figure
plot(x,y1,'b','linewidth',1.3);hold on;
%Plotting the potential profile to show the potential well
title('Plot for \Psi_{1}');
xlabel('x (in m)');
ylabel('\Psi_{1}(x)');
plot(x,v0,'r','linewidth',1.5);

%Plotting the wavefunction corresponding to the second allowed energy state
k2=1.452e9;
alpha2=5.995e9;
A2=15498.7;
C2=63994.35;
D2=A2;
F2=-3.94e14;
y2=zeros(1,1000);
for i=1:1000
    if x(i)<=0
        y2(i)=A2*exp(alpha2*x(i));
    elseif x(i)<=4e-9
        y2(i)=C2*sin(k2*x(i))+D2*cos(k2*x(i));
    elseif x(i)>4e-9
        y2(i)=F2*exp(-alpha2*x(i));
    end
end
figure
plot(x,y2,'b','linewidth',1.3);hold on;
title('Plot for \Psi_{2}');
xlabel('x (in m)');
ylabel('\Psi_{2}(x)');
%Plotting the potential profile to show the potential well
plot(x,v0,'r','linewidth',1.5);

