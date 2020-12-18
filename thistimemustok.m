clc,clear
%%
%������Χ��45-287kp, R2=0.9887
%D0
kp = 9;
if kp>=45
    D0 = 0.7351 * exp(-0.004 * kp) * 10^-8;
else
    D0 = 9.6174 * kp^-0.641 *10^-8;
end
D0 = D0 / 7;
%%
%��ѹVm
U = 2500;%��ѹ(V)
d = 0.2 * 10^-2;%����(m)
Eext = U / d;%��ӵ糡ǿ��(V/m)
t = 0 : 10^-8 : 4 * 10^-5;%ʱ�䷶Χ
r = 1.5 * 10^-6;%ϸ���뾶
Cm = 6.5 * 10^-3;%ϸ��Ĥ����
rou_i = 5.1^-1;%ϸ�����ʵĵ�����
rou_e = (1 / 12.55) * 10^4;%������Һ�ĵ�����(g)
g = 2;
h = 7 * 10^-9; % membrane thicknessϸ��Ĥ���7nm
tau = r * Cm * (rou_i + 0.5 * rou_e);%ϸ����Ĥ��λ����ӵ���̬���ӳ�ʱ��
Ri = 1 / (2 * g * r);
Rp = h / (pi * g * r * r);


Vm = (3/8) * Eext * r * (1 - exp(-t ./ tau));
% Vm = (1 + Ri / Rp) * Vp;

subplot(2,2,1);
plot(Vm)
title('Vm')
xlabel('t(ps)')
ylabel('��Ĥ��ѹ(V)')
grid on
%%
%Ĥ���濪���ܶ�
Vep = 0.258;%������Ĥ��ѹ
rp = 0.51 * 10^-9 ; %minimum radius of hydrophillic pores
rm = 0.8 * 10^-9; % minimum energy radius at Vm = 0
q = (rm / rp)^2;%�����Ǹ�����ֵ
a = 1 * 10^9; % creation rate coefficient
N0 = 1.5 * 10^9; % equillibrium pore denisty at Vm = 0

Nep = N0 * exp(q * (Vm ./ Vep).^2);
k1 = a * exp((Vm ./ Vep).^2);
N = k1 ./ (1 + t ./ Nep);

subplot(2,2,2);
plot(N)
title('���ܶ�')
xlabel('t(ps)')
ylabel('N')
grid on
%%
%�׾�
num = 4 * pi * r^2 * N;%���׸���
D = 5 * 10^-14 ; % Difussion coefficient for pore radius
T = 293.15;%�����¶�
Beta = 1.4 * 10^-19; % Steric repulsion energy
s1 = 2 * 10^-2;
s0 = 1 * 10^-3;
A = 1.26 * 10^-9;
rh = 0.97 * 10^-9;
rt = 0.31 * 10^-9;
Fmax = 0.7 * 10^-9;
gamma = 1.8 * 10^-11;
k = 1.381 * 10^-23;%������������

k1 = (4 * num .* Beta * rp^4 * D .* t) ./ k * T;
k2 = 2 * num * pi * gamma .* t;
k3 = 4 * num * pi * s1 .* t;
k4 = num .* Vm .* Fmax .* t;
k5 = num * pi / A;
k6 = num * pi * (2 * s1 - s0) .* t;

temp = size(t);
root = zeros(1,temp(2));

for i = 1:temp(2)
    f = @(r) k1(i)/r^5 + k2(i) - k3(i)*r + k4(i)/(1+rh/(r+rt)) +...
        k6(i)*( (2*r)/(k5(i)*r*r-1)^2 - (4*k5(i)*r^3)/(k5(i)*r*r-1)^3) - r;    %������ĺ���������ֱ���ڱ������޸ĺ������Ϊ��������
    U = 1 * 10^-6;%�Ͻ�
    L = 1 * 10^-10;%�½�
    while U - L > 0.0000000001    %�趨һ��������򾫶ȣ�Ȼ������ж�
        root(i) = (U + L) / 2;    %���������������������ʱ�����ö��ַ����¹滮�������
        if f(root(i)) == 0
            break;    %rǡ��Ϊ�������ֱ������ѭ��
        end
        if f(root(i)) * f(U) < 0    %�������ڶ����жϸ����ڵ�����
            L = root(i);
        else
            U = root(i);
        end
    end
    %root����������ֵ
end
R = root;

subplot(2,2,3);
plot(R)
title('ƽ���װ뾶(m)')
% ylim([0 1*10^-7])
xlabel('t(ps)')
ylabel('r')
grid on


%%
S = 4 * pi * R.^2 .* num;%�������
DNAo = 1.3 * 10^-11;%��ҺDNAŨ��
Zeff = 1;%DNA��Ч�۾���ֵ����ǰ��,���������϶�����һ�������������Ϊ������ֵ����ֱ��Ӱ��DNA��Ũ�ȡ�����zeff=20����ô
%����Ӱ���ž���ԭ����ֵ��10~20����Ҳ����һ���㡣
k = 1.381 * 10^-23;%������������
e = 1.6 * 10^-19;%Ԫ���
Vcell = 3 * pi * 10^-18; %ϸ���������������


k1 = D0/(h*Vcell);
k2 = (DNAo * Zeff * e) / (k * T);
DNAi = (k1 * S .* t .* (DNAo + k2 * Vm)) ./ (1 + k1 * S .* t);
% DNAi = (k1 * S .* t .* (DNAo )) ./ (1 + k1 * S .* t);

subplot(2,2,4)
plot(DNAi)
title('DNAת��Ũ��(mol/L)')
xlabel('t(ps)')
ylabel('DNAi')
grid on
