function [k] = K(C, TNAW, p)
    
if p.TNAW_const == 1
    a = ((2-log(2))/2) + 0.01*(p.TNAW - 10);

else
    a = ((2-log(2))/2) + 0.01*(TNAW - 10);
end
    k = exp(2*(C-a));
end