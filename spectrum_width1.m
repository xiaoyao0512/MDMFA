function width = spectrum_width1(alpha, falpha, falpha2, n)

	% n: number of points to be selected in the spectrum
    % usage: spectrum_width1(alpha_y, f_alpha_y, 5)
    
    
	% find the "appropriate" minimum of the spectrum between two spectrums
	min_val = min_spectrum(falpha, falpha2);
    width = 0;
    
	if (min_val == 0)  
        return
    end
	
	% find the interval
    intval = min_val / n;
	
	% find the width for each interval
    for i = 1:n-1
       width = width + spectrum_width_delta(alpha, falpha, i*intval);
    end
	width = width / (n-1);
    %width = spectrum_width_delta(alpha, falpha, min_val * (n-1)/n);
end
