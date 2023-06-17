function [causality, pval] = mdmfa(x, y, N)

    % Input: x and y are time series of length N
    % Output: if x causes y, then causality is 1. 0 otherwise. 
    %         a p value is also given.

    Qmax = 5;
    delta_q = 1e-1;
    Q = -Qmax:delta_q:Qmax;

    scmin = 8;
    scmax = floor(N / 4); 
    scres = 20; 

    exponents = linspace(log2(scmin), log2(scmax), scres);
    Scales = round(2 .^ exponents);
	
	detrend_fluctuations = true;
	detrend_order = 1;

    mf_params = struct(...
                'detrend', detrend_fluctuations, ...
                'detrend_order', detrend_order,...
                'scale_min', scmin,...
                'scale_max', scmax,...
                'scales', Scales,...
                'Q', Q);

    delta = 1;
    x = (x - min(x)) ./ (max(x) - min(x));
    y = (y - min(y)) ./ (max(y) - min(y));
    yt = [y(delta+1:end), y(1:delta)];

    out_space = {yt};
    condition_space = {y};
    [alpha_y, f_alpha_y, tau_q_y] = conditional_spectrum(out_space, condition_space,...
                                mf_params);
                            
    out_space = {yt};
    condition_space = {y, x};
    [alpha_xy, f_alpha_xy, tau_q_xy] = conditional_spectrum(out_space, condition_space,...
                                mf_params);
                            
    width_y = spectrum_width1(alpha_y, f_alpha_y, f_alpha_xy, 20);
    width_xy = spectrum_width1(alpha_xy, f_alpha_xy, f_alpha_y, 20);

    if (width_y > width_xy)
        causality = 1;
    else
        causality = 0;
    end

    pval = (width_y - width_xy) / width_y;

end
