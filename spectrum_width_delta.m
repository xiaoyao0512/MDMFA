function width = spectrum_width_delta(alpha, falpha, delta)

	% delta: fine grained y interval used in the interpolation (number)
	
	% split the spectrum into two halves from the maximum value
    [max_val, max_pos] = max(falpha);
	right = 1:(max_pos-1);
	left = (max_pos+1):length(falpha);
	alpha_left = alpha(left);
	alpha_right = alpha(right);
	falpha_left = falpha(left);
	falpha_right = falpha(right);
	
	% reference y value
	y_ref = max_val - delta;
	
	% find the x value at the reference y value
	x_ref_left = ref_point(alpha_left, falpha_left, y_ref);
	x_ref_right = ref_point(alpha_right, falpha_right, y_ref);
	width = abs(x_ref_right - x_ref_left);
end
