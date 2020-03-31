function [est_ent, est_tc, est_dtc, est_oinfo, entropy_time, tc_time, dtc_time, oinfo_time, temp_time] = information_metrics_calc(data, bias_corr)
% Estimation for entropy('est_ent'), total correlation('est_tc'),
% dual total correlation ('est_dtc') and O-info ('est_oinfo'),
% calculated from 'data' (T samples x N dimensionmatrix) 
% using empirical copulas (see data2gaussian method).
% 
% INPUT
% data = T samples x N variables matrix
%
% OUTPUT
% est_ent = Estimated total system entropy from 'data'
% est_tc = Estimated total correlation from 'data'
% est_dtc = Estimated dual total correlation from 'data'
% est_oinfo = Estimated oinfo from 'data'
% entropy_time = Time spent calculating system entropy
% tc_time = Time spent calculating total correlation
% dtc_time = Time spent calculating dual total correlation
% oinfo_time = Time spent calculating O-info
% temp_time = Time spent calculating the marginal entropies
%

% Obtaining dimensions for 'data'
[T,N] = size(data);
all_min_1=(arrayfun(@(x) setdiff(1:N,x),1:N,'uni',0)');

% Generating Biases for the information metrics
if bias_corr
    biascorrN = gaussian_ent_biascorr(N,T);
    biascorrNmin1 = gaussian_ent_biascorr(N-1,T);
    biascorr_1 = gaussian_ent_biascorr(1,T);
else
    biascorrN = 0;
    biascorrNmin1 = 0;
    biascorr_1 = 0;
end


% Transforming data to Gaussian Copulas
[~,est_covmat] = data2gaussian(data);

% Computing estimated measures for multi-variate gaussian

tic;
detmv = log(det(est_covmat));
detmv_min_1=log(cellfun(@(x) det(est_covmat(x,x)),all_min_1));
single_vars = diag(est_covmat);
% estimating info measures
est_var_ents= 0.5.*log(((2*pi*exp(1))).*single_vars) - biascorr_1;
est_ent_min_one = 0.5.*(log((2*pi*exp(1)).^(N-1)) + detmv_min_1) - biascorrNmin1;
temp_time = toc;

% Output Calculation

tic;
est_ent = 0.5.*(log(((2*pi*exp(1)).^N)) + detmv) - biascorrN;
entropy_time = toc;

tic;
est_tc = sum(est_var_ents) - est_ent;
tc_time = toc;

tic;
est_dtc = sum(est_ent_min_one) - (N-1).*est_ent;
dtc_time = toc;

tic;
est_oinfo = est_tc - est_dtc;
oinfo_time = toc;

entropy_time = entropy_time ; %+ temp_time;
tc_time = tc_time ; %+ temp_time;
dtc_time = dtc_time; % + temp_time;
oinfo_time = oinfo_time;% + (2*temp_time);