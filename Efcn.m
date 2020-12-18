function [check,isterminal,direction] = Efcn(t,y)

isterminal = ones(num+1, 1);
direction = -1*ones(num+1, 1);  %only when decreasing from positive to negative
check = ones(num+1 , 1); 


% check(1) = 1 - y(1)*(4*pi*(10*10^-6)^2); % checks if N*A is greater than 1
check(1) = 1 - y(1)*(4*pi*(10*10^-6)^2); % checks if N*A is greater than 1
% check(1) = y(1)*(4*pi*(7*10^-6)^2) - 1*(1+eps);

check(2:end) = 0.8*10^-9 - y(3:end); % Absorb pores that have fallen below rm
for m = 3:y(end)
    if y(m)== 0 
        check(m-1) = 1;
    end 
end 



end