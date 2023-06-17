
function [Fxy, Hq, tq, hq, Dq] = MFDFA(x, Scales, Q)

% Hq:           q-order Hurst exponent
% tq:           q-order mass exponent 
% hq:           q-order singularity exponent
% Dq:           q-order dimension 
% Fq:           q-order scaling function
for ns = 1:length(Scales)
    I = 1:Scales(ns):length(x)-Scales(ns)+1;
    % I = 1:length(x)-Scales(ns)+1;
    for n = 1:length(I)
        Xn = x(I(n):I(n)+Scales(ns)-1);
        Fxyn(n) = F_xy(Xn);
    end
    Fxyn = Fxyn ./ sum(Fxyn);
    for j = 1:length(Q)
        q = Q(j);
        %if q == 0
        %    Fxy(j, ns) = exp(mean(log(abs(Fxyn))));
                %    Fyx(j, ns) = exp(mean(log(abs(Fxyn))));
        %else
        Fxy(j,ns) = sum( (abs(Fxyn).^q) );
        %end
    end
    clear Fxyn
end
for nq=1:length(Q),
    C = polyfit(log2(Scales),log2(Fxy(nq,:)),1);
    tq(nq) = C(1);
end
Hq = (tq + 1) ./ Q;
hq = diff(tq)./(Q(2)-Q(1));
Dq = (Q(1:end-1).*hq)-tq(1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fxyn = F_xy(Xn)
% Calculate the detrended covariance of the n-th box
%ResX = DetrendedResiduals(Xn); 
%ResY = DetrendedResiduals(Yn); 
Fxyn = sum(Xn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Res = DetrendedResiduals(z)
% P = POLYFIT(X,Y,N) finds the coefficients of a polynomial P(X) of
% degree N that fits the data Y best in a least-squares sense. P is a
% row vector of length N+1 containing the polynomial coefficients in
% descending powers, P(1)*XˆN + P(2)*Xˆ(N-1) +...+ P(N)*X + P(N+1).
t = reshape([1:length(z)],size(z));
P = polyfit(t, z, 1);
trend_z = polyval(P,t);
Res = z - trend_z;