function x_ref = ref_point(alpha, falpha, y_ref)
	
	% two_points: Return x_ref from the curve at the y_ref point

	% pick 2 points to fit a polynomial
	
	% find the x value that is closest to y_ref
	falpha_diff = abs(falpha - y_ref);
	[minVal, minInd] = min(falpha_diff);
	
	% find the next closest one
	sec_min = 0;
    if (minInd == 1)
        sec_min = length(alpha);
        x_ref = 0;
    elseif minInd == length(alpha)
        sec_min = 1;
        x_ref = 0;
    elseif (abs(y_ref - falpha(minInd - 1)) > abs(y_ref - falpha(minInd + 1)))
		sec_min = minInd + 1;
	else
		sec_min = minInd - 1;
    end
    
	
	alpha_p = [alpha(minInd), alpha(sec_min)];
	falpha_p = [falpha(minInd), falpha(sec_min)];
	
	% linear fitting
	p = polyfit(alpha_p, falpha_p, 1);
	
	% find x_ref from the model and y_ref
	x_ref = (y_ref - p(2)) / p(1);

end