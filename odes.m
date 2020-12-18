function rk = odes(t, y, T, V0)
%     function rk = odes(t, y)
D0 = 2*10^-12 ; % DNA diffusion coefficient��ɢϵ��
R = 10*10^-6; % cell radiusϸ���뾶
H = 5*10^-9; % membrane thicknessϸ��Ĥ���h
Vcell = 3.14*10^-16; % cell volumeϸ�����
C0 = 2*10^-6; % Outside DNA concentraion����DNAŨ��DNAo
zeff = -27; % effective valence of DNA��Ч��
q0 = 1.602*10^-19; % elementary chargeԪ���e
kB = 1.381*10^-23; % Boltzman constant������������k
T0 = 310; % Absolute Temperature (37 C)���Ͼ����¶�T
g = 2; % Conductivity of the solution������Һ�絼��
f = 1.5; % form factor��״����
dVs = 1; % Potential Threshold
R0 = 100; % external resistance (series resistance of experimental setup)�ⲿ����Ri
Rm = 0.523; % Surface resistance of the membraneĤ�ϵ���
Cm = 9.5*10^-3; % Surface Capacitance of the membraneĤ�ϵ���
F_max = 0.7*10^-9; % Maximum electric force for Vm = 1V���糡��
rh = 0.97*10^-9; % Constant
rt = 0.31*10^-9; % Constant
a = 1*10^9; % creation rate coefficient��������alpha
Vep = 0.258; % characteristic voltage of electroporation�紩�׵�������ѹ
N0 = 1.5*10^9; % equillibrium pore denisty at Vm = 0
rp = 0.51*10^-9 ; %minimum radius of hydrophillic pores��ˮ��ˮ�׾���ֵr*
rm = 0.8*10^-9; % minimum energy radius at Vm = 0
D = 5*10^-14 ; % Difussion coefficient for pore radius�װ뾶����ϵ��
Beta = 1.4*10^-19; % Steric repulsion energy֬��Ĥ���ų���
gamma = 1.8*10^-11; %Edge Energy�׵ı�Ե����
s0 = 1*10^-3; %Tension of the bilayer without pores
s1 = 2*10^-2; %Tension of hydrocarbon-water interface
tau = 2.46; % electroporation constant��Ӿ����
kf = 2*10^-9; %coefficient
k0 = 1*10^10; %coefficient
q = (rm/rp)^2; % model paramater
A = 4*pi*R^2; % Total Surface Area
C = Cm*A; % total capacitance
R = Rm/A; % total resistance
Ip = 0 ; % initialize combined current through all pores
% K = 0; % initialize number of large pores in a given itteration
K1 = 0; % initialize total number of large pores
%%

%  format long

V0 = interp1(T, V0, t);

% dN/dt
rk(1) = a*exp((y(2)/Vep)^2)*(1 - (y(1)/(N0*exp(q*(y(2)/Vep)^2))));

%         Ip = (((y(2)/((H/(pi*g*(rm^2))) + (1/(2*g*(rm)))))*y(1)*A) + 2*pi*g*y(2)*sum((y(3:end).^2)./(2*pi + pi*y(3:end))));

% dVm/dt
rk(2) = ((V0/R0) - (((y(2)/((H/(pi*g*(rm^2))) +(1/(2*g*(rm)))))*y(1)*A) +...
    2*pi*g*y(2)*sum((y(3:end).^2)./(2*pi + pi*y(3:end)))) - (((1/R0) + (1/R))*y(2)))/C;

%         rk(2) = ((V0/R0) - (((y(2)/((H/(pi*g*(rm^2))) + (1/(2*g*(rm)))))*y(1)*A) + 2*pi*g*y(2)*sum((y(3:end).^2)./(2*pi + pi*y(3:end)),'all')) - (((1/R0) + (1/R))*y(2)))/C; % dVm/dt

for j = 3:K1 + 2
    %drj/dt
    
    rk(j) = (-D/kB*T0)*((-4*Beta*rp*(y(j))^-5) + 2*pi*gamma - 2*pi*2*s1 -...
        ((2*s1 - s0)/(1 - (((y(1))*A*pi*(rm^2) + (sum(((y(3:K1)).^2)))*pi)/A)))*...
        ((y(1))*A*pi*(rm^2) + (sum(((y(3:K1)).^2)))*pi)*y(j) + (F_max*((y(2))^2)/(1 + (rh/(y(j) + rt)))));
    
    % rk(j) = (-D/kB*T0)*((-4*Beta*rp*(y(j))^-5) + 2*pi*gamma - 2*pi*2*s1 - ((2*s1 - s0)/(1 - (((y(1))*A*pi*(rm^2) + (sum(((y(3:K1)).^2),'all'))*pi)/A)))*((y(1))*A*pi*(rm^2) + (sum(((y(3:K1)).^2),'all'))*pi)*y(j) + (F_max*((y(2))^2)/(1 + (rh/(y(j) + rt)))));
    
end

for v = 3:num+2
    if y(v) == 0
        rk(v) = 0; % drj/dt = 0 for pores that have not been created, or that have been absorbed
    end
end

rk = rk(:);

end