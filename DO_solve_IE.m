function [t,y,starfuncs, Ig, Cdot] = DO_solve_IE(y0, p) 
%y(1) = tg, y(2) = tnaw, y(3) = phi, y(4)= Ia, y(5) = lambda, y(6) = C

%define number of integration steps, depends on dt as we span 120Kyrs
    N = round(p.years/p.dt);

    % Interpolate for integration later
    xq = 1:p.dt:121; %query points
    
    scaledinterpB = interp1(1:121, p.scaledbenthic, xq);
    
    interpB = interp1(1:121, p.filteredbenthic, xq);
    
    % Tg* time series
    tgdiff= -37.39--47.39;
    tgstar = (scaledinterpB*tgdiff) + (-37.83-(tgdiff*scaledinterpB(1)));
    
    %Ig* time series
    igdiff = -38.29--45.43;
    igstar = (scaledinterpB*igdiff) + (-38.57-(igdiff*scaledinterpB(1)));
    
    % Tnaw* time series
    tnawdiff = 16-6;
    tnawstar = (scaledinterpB*tnawdiff) - (scaledinterpB(1) *tnawdiff) + 15.68;
    
    % Ia* time series
    iadiff = -31.07--42.68;
    iastar = (scaledinterpB*iadiff) + (-31.25 -(iadiff*scaledinterpB(1)));
    starfuncs = [tgstar; tnawstar; iastar; igstar];
    
    %defining required input for rhs
    %y(1) = tg, y(2) = tnaw, y(3) = phi, y(4)= Ia, y(5) = Ig, y(6) = lambda, y(7) = C
    % constant f*(t) function

    %initialise variable vector and time vector
    y = zeros(length(y0), N+1);
    t = zeros(1, N+1);
    Cdot = zeros(1, N+1);

    %put in initial conditions for each of the variables
    y(:,1)= y0;
    Cdot(1)=0;
    t(1) = 120;

    %initialise IG
    Ig = zeros(1, N+1);
    Ig(1) = p.Ig_0;
    
    for i = 1 : N
        
        % delta W_i
        dwi = randn(1);

        %current delta O18 value 
        p.delta = interpB(end+1-i);        
                     
        %implementing Heun's method of forward Euler

        k1 = rhs(y(:,i), i, starfuncs, p);
        Cdot(i+1) = k1(6);
        k2 = rhs(y(:,i) + p.dt*k1, i+1, starfuncs, p);

        y(:,i+1) = y(:,i) + (1/2)*p.dt*(k1+k2);

        lambda_tilda_next = y(5,i) + k1(5)*p.dt + p.sigma* (sqrt(p.dt)/sqrt(p.epsilon))*dwi;
        foo = y(:,i)+ p.dt*k1;
        foo(5) = lambda_tilda_next;

        RHS = rhs(foo, i+1, starfuncs, p);

        y(5, i+1)= y(5,i) + (1/2)*p.dt*(k1(5) + RHS(5)) + p.sigma*(sqrt(p.dt)/sqrt(p.epsilon))*dwi;

        %update Ig

        Ig(i+1) = starfuncs(4, end+1-i) + 0.01*y(1,i) + 4.32*(1-y(6,i)) - p.beta*y(5,i);

        %update T
   
        t(i+1) = t(i) - p.dt;
        
    end
end