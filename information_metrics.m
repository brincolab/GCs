function results = information_metrics(data, infodynamics_path, octave_path, varargin)
%
%
% Inputs
% - data = matrix to obtain information metrics
% - infodynamics_path = path to JIDT infodynamics.jar 
% - octave_path = patho to JIDT octave bindings
% - varargin = Unused, for future implementation of flags
% 
% Outputs 
% - results = Results obtained from information_metrics_calc and information_metrics_calc_JIDT methods 
%

results = containers.Map;

%Check if covariance matrix is present on the workspace
flag_covmat = false;
if exist('true_covmat','var')
    flag_covmat = true;
end

%Calculate the teoretical information metrics if covmat is present
if flag_covmat
    [T,N] = size(data);
	all_min_1=(arrayfun(@(x) setdiff(1:N,x),1:N,'uni',0)');
	detmv = log(det(true_covmat));
	detmv_min_1=log(cellfun(@(x) det(true_covmat(x,x)),all_min_1));
	single_vars = diag(true_covmat);
	
	teor_var_ents = 0.5.*log(((2*pi*exp(1))).*single_vars);
	teor_ent_min_one = 0.5.*(log((2*pi*exp(1)).^(N-1)) + detmv_min_1);
	
	teor_ent = 0.5.*(log(((2*pi*exp(1)).^N)) + detmv);
	teor_tc = sum(teor_var_ents) - teor_ent;
	teor_dtc = sum(teor_ent_min_one) - (N-1).*teor_ent;
	teor_oinfo = teor_tc - teor_dtc;
	
	results('teor_ent') = teor_ent;
	results('teor_tc') = teor_tc;
	results('teor_dtc') = teor_dtc;
	results('teor_oinfo') = teor_oinfo;
end

%Calculate Gaussian Copula information metrics of the data
[est_ent, est_tc, est_dtc, est_oinfo, entropy_time_gc, tc_time_gc, dtc_time_gc, oinfo_time_gc, temp_time] = information_metrics_calc(data, bias_flag);

%Calculate JIDT gaussian information metrics of the data
c_type = "gaussian";
[JIDT_ent_gaussian, JIDT_tc_gaussian, JIDT_dtc_gaussian, JIDT_oinfo_gaussian, entropy_time_gaussian, tc_time_gaussian, dtc_time_gaussian, oinfo_time_gaussian] = information_metrics_calc_JIDT(data, infodynamics_path, octave_path,c_type);

%Calculate JIDT kozachencko-kraskov information metrics of the data
c_type = "kozachenko_kraskov";
[JIDT_ent_kozachenko, JIDT_tc_kraskov, JIDT_dtc_kraskov, JIDT_oinfo_kraskov, entropy_time_kozachenko, tc_time_kraskov, dtc_time_kraskov, oinfo_time_kraskov] = information_metrics_calc_JIDT(data, infodynamics_path, octave_path,c_type);

% Save the results
results('est_ent') = est_ent;
results('est_tc') = est_tc;
results('est_dtc') = est_dtc;
results('est_oinfo') = est_oinfo;
results('est_ent_time') = entropy_time_gc;
results('temp_time') = temp_time;
results('est_tc_time') = est_tc;
results('est_dtc_time') = est_dtc;
results('est_oinfo_time') = est_oinfo;


results('JIDT_ent_gaussian') = JIDT_ent_gaussian;
results('JIDT_tc_gaussian') = JIDT_tc_gaussian;
results('JIDT_dtc_gaussian') = JIDT_dtc_gaussian;
results('JIDT_oinfo_gaussian') = JIDT_oinfo_gaussian;

results('JIDT_ent_gaussian_time') = entropy_time_gaussian;
results('JIDT_tc_gaussian_time') = tc_time_gaussian;
results('JIDT_dtc_gaussian_time') = dtc_time_gaussian;
results('JIDT_oinfo_gaussian_time') = oinfo_time_gaussian;

results('JIDT_ent_kozachenko') = JIDT_ent_kozachenko;
results('JIDT_tc_kraskov') = JIDT_tc_kraskov;
results('JIDT_dtc_kraskov') = JIDT_dtc_kraskov;
results('JIDT_oinfo_kraskov') = JIDT_oinfo_kraskov;

results('JIDT_ent_kozachenko_time') = entropy_time_kozachenko;
results('JIDT_tc_kraskov_time') = tc_time_kraskov;
results('JIDT_dtc_kraskov_time') = dtc_time_kraskov;
results('JIDT_oinfo_kraskov_time') = oinfo_time_kraskov;

if flag_covmat
	 results('GC_ent_error') = teor_ent - est_ent;
	 results('GC_tc_error') = teor_tc - est_tc;
	 results('GC_dtc_error') = teor_dtc - est_dtc;
	 results('GC_oinfo_error') = teor_oinfo - est_oinfo;
	 
	 results('JIDT_gaussian_ent_error') = teor_ent - JIDT_ent_gaussian;
	 results('JIDT_gaussian_ent_error') = teor_tc - JIDT_tc_gaussian;
	 results('JIDT_gaussian_ent_error') = teor_dtc - JIDT_dtc_gaussian;
	 results('JIDT_gaussian_ent_error') = teor_oinfo - JIDT_oinfo_gaussian;
	 
	 results('JIDT_kozachenko_ent_error') = teor_ent - JIDT_ent_kozachenko;
	 results('JIDT_kozachenko_ent_error') = teor_tc - JIDT_tc_kozachenko;
	 results('JIDT_kraskov_ent_error') = teor_dtc - JIDT_dtc_kraskov;
	 results('JIDT_kraskov_ent_error') = teor_oinfo - JIDT_oinfo_kraskov;
end