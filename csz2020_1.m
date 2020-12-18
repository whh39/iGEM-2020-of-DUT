nbiCFUtot=100000;%the total bacterial load
S=1*10^(5);
R=1000;
CFUmax=1*10^(6);
i=8*10^(5);
q=0.01;%Fractional in-flow/out flow, h21
v=0.0108;%Resistant E. coli fraction in in-flow
b=0.004;%Plasmid transmission term, h-1
%CFUm;%the median effect bacterial concentration
kd=0.16;%Elimination of bacteria is by a first-order process according to a rate constant Kd
%VGmax=kd*(CFUm+CFUmax);%the maximal velocity of bacterial growth
%kg=VGmax/(CFUm+CFUtot);%the apparent rate constant for replication
e=0.2;
%kg=e*(1-(S+R)/CFUmax);
Emax=2;
MICsB=1;%Ceftiofur on LB
MICsD=8;%Ceftiofur on LB
C=8;
x=C;
N=S+R;
v1=0.01;
T=0:1:200;
for idex=1:length(T)-1
    E2=Emax*x(idex)^(1.5)/(MICsB^(1.5)+x(idex)^(1.5));
    E4=Emax*x(idex)^(1.5)/(MICsD^(1.5)+x(idex)^(1.5));
    if rem(idex,24)==0
        x(idex+1)=x(idex)*exp(-0.1)+8;
    else
        x(idex+1)=x(idex)*exp(-0.1);
    end
    S(idex+1)=S(idex)+S(idex)*(e*(1-((S(idex)+R(idex)+i)/(CFUmax)))*(1-E2))-b*R(idex)*S(idex)/(R(idex)+S(idex));%+(1-v)*q*(S(idex)+R(idex))-q*S(idex);
    if S(idex+1)>=0
        R(idex+1)=R(idex)+R(idex)*(e*(1-0.05)*(1-((S(idex)+R(idex)+i)/(CFUmax)))*(1-E4))+b*R(idex)*S(idex)/(R(idex)+S(idex));%+v*q*(S(idex)+R(idex))-q*R(idex);
    else
        S(idex+1)=0;
        %if R(idex)>=0;
        R(idex+1)=R(idex)+R(idex)*(e*(1-0.05)*(1-((S(idex)+R(idex))/(CFUmax)))*(1-E4))+b*R(idex)*S(idex)/(R(idex)+S(idex));
        %else
        %R(idex+1)=0;
        %end
    end
end
%%
subplot(2,3,6);
plot(T,x,'r','LineWidth', 1.5);
xlabel('Hours (h)');
ylabel('Concentration(¦Ì/ml)');
grid on
title('Ceftiofur concentration with intravenous administration per 24 hours')

subplot(2,3,[1 2 4 5]);
plot(T,R,T,S,T,S+R,'LineWidth', 1.5);
xlabel('Hours (h)');
ylabel('CFU/ml');
legend('R','S','T+S');
grid on
title('Pharmacokinetic-pharmacodynamic (PK-PD) model')
hold on
ylim([0 110000])
% for i = [25,49,73,98,122,146,169]
%     plot(T(i),R(i),'ro','MarkerFaceColor','b','color','b')
%     hold on
%     plot(T(i),S(i),'ro','MarkerFaceColor','r','color','r')
%     hold on
%     plot(T(i),R(i)+S(i),'ro','MarkerFaceColor','k','color','k')
%     hold on
% end





subplot(2,3,3);
plot(T,(R./(S+R)),'LineWidth', 1.5);
xlabel('Hours (h)');
grid on
title('Proporation of resistant bacterial')
ylabel('')