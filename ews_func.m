function [variance, lag1]= ews_func(w, x)

    %creating ’empty’ arrays to collect variance and
    %lag−1 autocorrelation
    %variance = zeros (1, length ( x(1,:) )-w);
    variance = zeros (1, length ( x(1,:) ));

    %lag1 = zeros (1, length ( x(1,:) )-w);
    lag1 = zeros (1, length ( x(1,:) ));

    for i = w:length(x(1,:))

        xw = x(1, i-w+1: i);

        %xw = detrend(xw, 1);
        xw = xw -mean(xw);

        variance (i)= var(xw);

        [acf, lags] = autocorr(xw);

        lag1(i) = acf(2);


    end
    variance(1, 1:w-1) = NaN;
    lag1(1, 1:w-1) = NaN;
    
end 