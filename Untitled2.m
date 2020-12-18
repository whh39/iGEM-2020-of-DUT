clc;clear;
b=5.3*10^(-6);%infection rate constant
k=4.0;%transition rate for infected cells to produce virus
a=3.8;%death rate of infected cells
p=1.3*10^(-1);%viral release rate per infected cell
c=3.8;%viral clearance rate
V=4*10^(8);
I1=0;
I2=0;
v=4.9;
T=0:1:100;
for idex=1:length(T)-1
    V(idex+1)=V(idex)-0.1*b*V(idex)*v(idex);
    I1(idex+1)=I1(idex)+0.1*b*V(idex)*v(idex)-0.1*k*I1(idex);
    I2(idex+1)=I2(idex)+0.1*k*I1(idex)-0.1*a*I2(idex);
    v(idex+1)=0.1*p*I2(idex)-0.1*c*v(idex);
end
plot(T/10,V,T/10,I1,T/10,I2,T/10,V+I1+I2);
legend('Target cells','Eclipse cells','Infected cells','Totality');
xlabel('Days');
ylabel('Cells');