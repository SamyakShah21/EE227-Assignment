%Prepeared by Samyak Shah for EE227 Assignment 1 Question 2
%Date: September 08,2019
%______________________________________________
%Defining the constants given in the question
m=9.1e-31;
hcut=1.054e-34;
v0=4.8e-19;
l=2e-9;
%Note: We were asked to draw the Tunneling probability for very small bias
%      The plot produced here is a more general one. We observe the initial
%      rise only in case of a very small bias
E=linspace(0,6.2e-19,1000);
alpha=zeros(1,1000);
T=zeros(1,1000);
%Finding tunneling probability for all values of E in the 'E' array
for i=1:1000
    alpha(i)=sqrt(2*m*(v0-E(i)))/hcut;
    T(i)=1/((1+0.25*(v0*v0/(E(i)*(v0-E(i))))*sinh(alpha(i)*l)*sinh(alpha(i)*l)));
end
figure
plot(E,T,'r','linewidth',1.3);
xlabel("E value (in J)");
ylabel('Tunneling Probability');
title('Tunneling probability v/s E');
