function min_intv = min_spectrum(falpha1, falpha2)

	len1 = length(falpha1);
	len2 = length(falpha2);
	% find leftmost or rightmost points;
	% select the one with higher value in the y axis
    min_idx1 = 0;
    if (falpha1(1) > falpha1(len1))
        min_idx1 = 1;
    else
        min_idx1 = len1;
    end
    min_idx2 = 0;
    if (falpha2(1) > falpha2(len2))
        min_idx2 = 1;
    else
        min_idx2 = len2;
    end
	
	%intval1 = max(falpha1) - falpha1(min_idx1);
	%intval2 = max(falpha2) - falpha2(min_idx2);
    intval1 = max(falpha1) - min(falpha1);
    intval2 = max(falpha2) - min(falpha2);
	min_intv = min(intval1, intval2);
	
    if (intval1 > 10*intval2)
        fprintf('too small\n');
        min_intv = 0;%max(intval1, intval2);
    end
end