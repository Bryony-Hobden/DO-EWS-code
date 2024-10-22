function [x, W]= bry_func(compiled_ts, sampling_rate, N_win)
%input in section of timeseries of time on top and isotope ratio on bottom,
%the smapling rate in kyrs and the number of data points you want in your
%sliding window

%output is detrended timeseries, time in kyrs, var, lag1, alpha, upper
%alpha, lower alpha
    dt = 0.0001;
    N_tot = length(compiled_ts);

    N_sub = sampling_rate /dt;
    dtau = N_sub*dt;

    y_n = zeros(1, floor(N_tot/N_sub));
    t_n = zeros(1, floor(N_tot/N_sub));
   

    for i = 1:length(y_n)
        y_n(i) =1/N_sub * sum(compiled_ts(2,((i-1)*N_sub +1):((i-1)*N_sub + N_sub)));
        t_n(i) = compiled_ts(1,1) - (i-1)*sampling_rate;
        %
    end

    W = N_win*dtau;

    y_detrend = y_n;

    [var, lag1] = ews_func(N_win,y_n);

    alpha = -log(max(lag1,0))/dtau;

    sd_lag1 = sqrt((2*alpha*dtau)/N_win);

    alpha_u = -log(max(lag1-sd_lag1,0))/dtau;
    alpha_l = -log(max(lag1+sd_lag1,0))/dtau;

    x = [y_detrend; t_n; var; lag1; alpha; alpha_u; alpha_l; y_n];
end 

