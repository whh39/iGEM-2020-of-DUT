clc
clear
%%
t1 = 0 : 0.01 : 0.48;
t2 = 0.48 : 0.01 : 4;
t3 = 0.48 : 0.01 : 1.94;
x1 = -3772.1 * t1 .^5 + 4193.6 * t1 .^4 - 1569.1 * t1.^3 + 405.65 * t1 .^2 - 5.3313 * t1 + 0.0935;
x2 =  1.5848 * t2 .^6 - 23.242 * t2 .^5 + 136.85 * t2 .^4 - 412.24 * t2 .^3 + 666.48* t2 .^2 - 522.88 * t2 + 181.54;
x3 = 233.43 * t3 .^6 - 1812.9 * t3 .^5 + 5736.3 * t3 .^4 - 9460.6 * t3 .^3 + 8571.4 * t3 .^2 - 4054.5 * t3 + 804.2;
m = 44.782220807727022;
x1(49) = m;
x2(1) = m;
x3(1) = m;
%%
plot(t1, x1,'LineWidth', 1.5,'color','k')
hold on
plot(t2, x2,'LineWidth', 1.5,'color','k')
hold on
plot(t3, x3,'--','LineWidth', 1.5,'color','k')
% xlim([0 4])
grid on
xlabel('Radius (nm)')
ylabel('Pore energy (units of kT)')

plot(0.48,x2(1),'ro','MarkerFaceColor','k','LineWidth', 1.5,'color','k')
plot(t2(40),x2(40),'ro','MarkerFaceColor','k','LineWidth', 1.5,'color','k')
%%
text(0.48,x2(1),'\leftarrow (r_*,E_*)','color','k','FontSize',18);
% text(0.48,x1(end),'o');
text(t2(40),x2(40),'\leftarrow (r_m,E_m)','color','k','FontSize',18);
% text(t2(40),x2(40),'o');
title('Energy of each pore corresponding to pore radius');


