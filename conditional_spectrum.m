function [alpha, f_alpha, tau_q] = conditional_spectrum(out, condition, params)

L = length(out);
Lc = length(condition);
if L ~= 1
    error('The output space should be 1D') % future extension for multiple dimensions
end

Fq = conditional_partition_function(out, condition, params);

Fq_cond = joint_partition_function(condition, params);

Q = params.Q;
scales = params.scales;
[~, ~, tau_q_cond] = get_tau_q(Fq_cond, Q, Lc, scales);

Fq_avg = get_averaged_partition(Fq, tau_q_cond, Lc, params);
[alpha, f_alpha, tau_q] = get_tau_q(Fq_avg, Q, L, scales);


function Fq_avg = get_averaged_partition(Fq, tau_q_cond, Lc, params)


Q = params.Q;
dq = Q(2)-Q(1);
len_Q = length(Q);
scales = params.scales;
ns = length(scales);

Fq_avg = zeros(len_Q, ns); % support for 1D out space for now
switch Lc

    case 1

        alpha = gradient(tau_q_cond, dq);
        d2tau2_dq2 = gradient(alpha, dq);
       
        f_alpha = alpha .* reshape(Q, size(alpha)) - tau_q_cond;

        for si = 1:ns
            s = scales(si);
            prod_ = s.^-f_alpha .* abs(d2tau2_dq2);
            if Lc>1
                prod_ = reshape(prod_, [1, size(prod_)]); % adding extra dimension
            end
            prod_ = repmat(prod_, len_Q, 1);
            normalization  = sum(prod_, 2) * dq;
            Fq_avg(:,si) = sum(Fq(:,:,si) .* prod_, 2) * dq ./ normalization;
            

        end

    case 2

        [alpha_1, alpha_2] = gradient(tau_q_cond, dq, dq);
        [d2tau2_dq1_dq1, d2tau2_dq2_dq1] = gradient(alpha_1, dq, dq);
        [d2tau2_dq1_dq2, d2tau2_dq2_dq2] = gradient(alpha_2, dq, dq);

        J = d2tau2_dq1_dq1 .* d2tau2_dq2_dq2 - d2tau2_dq1_dq2 .* d2tau2_dq2_dq1;
        [Q1, Q2] = meshgrid(Q, Q);
        f_alpha = alpha_1 .* Q1 + alpha_2.*Q2 - tau_q_cond;

        for si = 1:ns
            s = scales(si);
            prod_ = s.^-f_alpha .* abs(J);
            prod_ = reshape(prod_, [1, size(prod_)]); % adding extra dimension
            prod_ = repmat(prod_, len_Q, 1);
            normalization = sum(prod_ * dq * dq, [2,3]);
            Fq_avg(:,si) = sum(Fq(:,:,:,si) .* prod_, [2,3]) * dq * dq ./normalization;
        end
end


function [alpha, f_alpha, tau_q] = get_tau_q(Fq, Q, L, scales)

len_Q = length(Q);
if L>1
    tau_q = zeros(ones(1, L)*len_Q);
else
    tau_q = zeros(1, len_Q);
end

coord_grid = get_Q_grid(1:len_Q, L);
for l = 1:L
    coord_grid{l} = coord_grid{l}(:);
end
indices_ = fliplr(cat(2, coord_grid{:}));
for i = 1:size(indices_,1)
    ind = num2cell(indices_(i,:));
    C = polyfit(log(scales), log(Fq(ind{:},:)), 1);
    tau_q(ind{:})=C(1);
end

dq = Q(2)-Q(1);
switch L
    case 1
        alpha = gradient(tau_q, dq);
        f_alpha = alpha .* reshape(Q, size(alpha)) - tau_q;

    case 2
        [alpha_1, alpha_2] = gradient(tau_q, dq, dq);
        [Q1, Q2] = meshgrid(Q, Q);
        f_alpha = alpha_1 .* Q1 + alpha_2.*Q2 - tau_q;
        alpha = {alpha_1, alpha_2};
end


function Q_grid = get_Q_grid(Q, L)

Q_ = cell(1, L);
for l = 1:L
    Q_{l} = Q;
end
Q_grid = cell(size(Q_));
[Q_grid{:}] = ndgrid(Q_{:});


function Fq = conditional_partition_function(out_space, cond_space, params)

L = length(out_space);
Lc = length(cond_space);
Lt = L+Lc;
total_space = [out_space, cond_space];

Q = params.Q;
scales = params.scales;
ns = length(scales);
len_Q = length(Q);
Fq = zeros([ones(1, Lt)*len_Q, ns]);

log_fluctuations = cell(Lt, ns);
for l = 1:Lt % for each signal in the out space
    for si = 1:ns
        log_fluctuations{l, si} = log(get_fluctuations(total_space{l}, scales(si), params));
    end
end

switch Lc

    case 1
        [Q1_grid, Q2_grid] = meshgrid(Q, -Q);

        for si = 1:ns
            num_intervals = length(log_fluctuations{1,si});
            for ii = 1:num_intervals
                prod_ = log_fluctuations{1,si}(ii) * Q1_grid + ...
                        log_fluctuations{2,si}(ii) * Q2_grid;
                Fq(:,:,si) = Fq(:,:,si) + exp(prod_);
            end
        end

    case 2
        [Q1_grid, Q2_grid, Q3_grid] = meshgrid(Q, -Q, -Q);
        for si = 1:ns
            num_intervals = length(log_fluctuations{1,si});
            for ii = 1:num_intervals
                prod_ = log_fluctuations{1,si}(ii) * Q1_grid + ...
                        log_fluctuations{2,si}(ii) * Q2_grid + ...
                        log_fluctuations{3,si}(ii) * Q3_grid;
                Fq(:,:,:,si) = Fq(:,:,:,si) + exp(prod_);
            end
        end
end


function Fq = joint_partition_function(out_space, params)

L = length(out_space);

Q = params.Q;
scales = params.scales;
ns = length(scales);
len_Q = length(Q);
Fq = zeros([ones(1, L)*len_Q, ns]);

fluctuations = cell(L, ns);
for l = 1:L % for each signal in the out space
    for si = 1:ns
        fluctuations{l, si} = get_fluctuations(out_space{l}, scales(si), params);
    end
end

num_Q = length(Q);
Q_ = cell(1, L);
for l = 1:L
    Q_{l} = Q;
end
Q_grid = cell(size(Q_));
[Q_grid{:}] = ndgrid(Q_{:});

for si = 1:ns
    num_intervals = length(fluctuations{1,si});
    for ii = 1:num_intervals
        prod_ = ones(size(Q_grid{1}));
        for l = 1:L
            prod_ = times(prod_, fluctuations{l, si}(ii).^Q_grid{l});
        end
        Fq(1+(si-1)*(num_Q^L): si*num_Q^L) = Fq(1+(si-1)*(num_Q^L): si*num_Q^L) ...
            + reshape(prod_, 1, num_Q^L);
    end
end


function fluc = get_fluctuations(x_, scale, params)

N = length(x_);
num_intervals = floor(N/scale);
intervals = [(1:scale:num_intervals*scale)', (scale:scale:num_intervals*scale)'];

fluc = zeros(num_intervals,1);
for ii = 1:num_intervals
    x_i = x_(intervals(ii,1):intervals(ii,2));
    if params.detrend
        C = polyfit(intervals(ii,1):intervals(ii,2), x_i, params.detrend_order);
        fluc(ii) = sqrt(mean((x_i - polyval(C, intervals(ii,1):intervals(ii,2))).^2));
    else
        fluc(ii) = (x_i(end) - x_i(1)).^2;
    end
end