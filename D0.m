clc
clear
%%
subplot(1,3,1)
size = [1.868,2.100,2.311,2.326,2.686,2.717,2.798,3.730,4.363,4.752,5.400,6.500,10.200,11.100,45.000,112.800,183.500,287.100];
D = [7.00,6.00,5.60,5.35,5.40,4.90,4.80,4.11,3.70,3.51,3.10,2.89,2.30,1.65,0.60,0.49,0.33,0.23];
scatter(size,D,'*');
plot(size,D,'*');
title('Plasmid diffusion coefficient data collection from previous article');
xlabel('Size (kbp)');
ylabel('D_0 \times 10^8');
legend('original data point')
grid on
%%
subplot(1,3,2)
size = [1.868,2.100,2.311,2.326,2.686,2.717,2.798,3.730,4.363,4.752,5.400,6.500,10.200,11.100];
D = [7.00,6.00,5.60,5.35,5.40,4.90,4.80,4.11,3.70,3.51,3.10,2.89,2.30,1.65];
scatter(size,D,'*');
plot(size,D,'*');
hold on
t=1:0.1:12;
plot(t,9.6174 * t .^-0.641,'LineWidth',1.5 )
hold on
title('Curve-fiting in low domain of plasmid size');
xlabel('Size (kbp)');
ylabel('D_0 \times 10^8');
legend('original point','curve-fitting');
grid on
%%
subplot(1,3,3)
size = [45.000,112.800,183.500,287.100];
D = [0.60,0.49,0.33,0.23];
scatter(size,D,'*');
plot(size,D,'*');
hold on
t=10:0.1:300;
plot(t,0.7351 .* exp(-0.004 .* t),'LineWidth',1.5 )
hold on
title('Curve-fiting in high domain of plasmid size');
xlabel('Size (kbp)');
ylabel('D_0 \times 10^8');
legend('original point','curve-fitting');
grid on