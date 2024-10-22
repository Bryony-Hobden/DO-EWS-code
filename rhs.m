function[y_out] = rhs(y, t, starfuncs, p)

    %starfuncs = [tgstar; tnawstar; iastar; igstar];
    %tgstar(t) = starfuncs(1,t)
    %tnawstar(t) = starfuncs(2,t)
    %iastar = starfuncs(3,t)
    %igstar = starfuncs(4,t)

    % t just i in DO_solve

    % function for RHS of equation 

    %threshold for Tnaw, if Tnaw> tau, ice sheet destabilized and C=0
    
    pp = [5.56, 0.01, 4.32, 11.49];

    %y(1) = tg, y(2) = tnaw, y(3) = phi, y(4)= Ia, y(5) = lambda, y(6) = C

    %Tg
    y_out1 = (- pp(1)*(y(1)- starfuncs(1,end+1-t))+ 10*pp(1)*y(2)*Theta(y(2))*(1-y(6))); 
    %Tnaw
    y_out2 = ( -pp(1)*(y(2) - starfuncs(2,end+1-t)) - 10*pp(1)*y(2)*Theta(y(2))*(1 - y(6)) - pp(1)*y(3));
    %phi
    y_out3 = (y(3) - (y(3))^3 - pp(1)*(y(2) - p.tau/2));
    %Ia
    y_out4 = -(y(4) - (starfuncs(3,end+1-t) - y(3))) ;
    %lambda
    y_out5 = (1/p.epsilon)*(-y(6) -(8/135)*(y(5))^3 +(4/45)*(y(5))^2 +(16/45)*y(5)+ (11/27))*(2*(y(5)-0.5)*(y(5)-1.5)+y(6));
    %C
    y_out6 =  H(y(5), p, y(6))*(y(5) - K(y(6), y(2), p));
    
    %final rhs output- column vector 
    y_out = [y_out1; y_out2; y_out3; y_out4; y_out5; y_out6];
  
end
    
