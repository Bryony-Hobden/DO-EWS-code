function [outputArg1] = Theta(x)
%my own heaviside function 
%if x<=0 then Theta(x)=0
%if x>0 then Theta(x)=1

if x<= 0
    outputArg1  = 0;
else 
    outputArg1 = 1;
end