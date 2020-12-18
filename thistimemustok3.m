clc,clear

for z1 = 1 : 94
    for z2 = 1 : 401
        
        z1
        
        kp = 1.8 + 0.1 * (z1 - 1);
        D0 = 9.6174 * kp^-0.641 *10^-8;
        D0 = D0 / 7;
       %%
        %电压Vm
        % U = 1500; %电压(V)
        U = 1000 + (z2 - 1) * 5;
        d = 0.2 * 10^-2;%板间距(m)
        Eext = U / d;%外加电场强度(V/m)
        t = 0:10^-8:4*10^-5;%时间范围
        r = 1.5 * 10^-6;%细胞半径
        Cm = 6.5 * 10^-3;%细胞膜电容
        rou_i = 5.1^-1;%细胞基质的电阻率
        rou_e = (1/12.55) * 10^4;%胞外溶液的电阻率(g)
        g = 2;
        h = 7 * 10^-9; % membrane thickness细胞膜厚度7nm
        tau = r * Cm * (rou_i + 0.5 * rou_e);%细胞跨膜电位从零加到稳态的延迟时间
        Ri = 1 / (2 * g * r);
        Rp = h / (pi * g * r * r);
        
        
        Vm = (3/8) * Eext * r * (1 - exp(-t./tau));
        % Vm = (1 + Ri / Rp) * Vp;
        
        %%
        %膜表面开孔密度
        Vep = 0.258;%特征跨膜电压
        rp = 0.51 * 10^-9 ; %minimum radius of hydrophillic pores
        rm = 0.8 * 10^-9; % minimum energy radius at Vm = 0
        q = (rm / rp)^2;%反正是个特征值
        a = 1 * 10^9; % creation rate coefficient
        N0 = 1.5 * 10^9; % equillibrium pore denisty at Vm = 0
        
        Nep = N0 * exp(q * (Vm./Vep).^2);
        k1 = a * exp((Vm./Vep).^2);
        N = k1 ./ (1 + t./Nep);
        
        %%
        %孔径
        num = 4 * pi * r^2 * N;%开孔个数
        D = 5 * 10^-14 ; % Difussion coefficient for pore radius
        T = 293.15;%开氏温度
        Beta = 1.4 * 10^-19; % Steric repulsion energy
        s1 = 2 * 10^-2;
        s0 = 1 * 10^-3;
        A = 1.26 * 10^-9;
        rh = 0.97 * 10^-9;
        rt = 0.31 * 10^-9;
        Fmax = 0.7 * 10^-9;
        gamma = 1.8 * 10^-11;
        k = 1.381 * 10^-23;%玻尔兹曼常数
        
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
                k6(i)*( (2*r)/(k5(i)*r*r-1)^2 - (4*k5(i)*r^3)/(k5(i)*r*r-1)^3) - r;    %求给定的函数，可以直接在本行中修改后面代码为其他函数
            U = 1 * 10^-6;%上界
            L = 1 * 10^-10;%下界
            while U - L > 0.0000000001    %设定一个求根区域精度，然后进行判断
                root(i) = (U + L) / 2;    %当根的区间大于所给精度时，利用二分法重新规划求根区间
                if f(root(i)) == 0
                    break;    %r恰好为所求根，直接跳出循环
                end
                if f(root(i)) * f(U) < 0    %用零点存在定理判断根所在的区域
                    L = root(i);
                else
                    U = root(i);
                end
            end
            %root输出所求根的值
        end
        R = root;
        
        %%
        S = 4 * pi * R.^2 .* num;%开孔面积
        DNAo = 1.3 * 10^-6;%溶液DNA浓度
        Zeff = 1;%DNA有效价绝对值
        k = 1.381 * 10^-23;%玻尔兹曼常数
        e = 1.6 * 10^-19;%元电荷
        Vcell = 3 * pi * 10^-18; %细胞体积
        
        k1 = D0/(h*Vcell);
        k2 = (DNAo * Zeff * e) / (k * T);
        DNAi = (k1 * S .* t .* (DNAo + k2 * Vm)) ./ (1 + k1 * S .* t);
        
        %%
        DNA(z1,z2) = DNAi(end);
    end
end

ans=(DNA/(1.77*10^-7))';