function [est_ent, est_tc, est_dtc, est_oinfo, entropy_time, tc_time, dtc_time, oinfo_time] = information_metrics_calc_JIDT(data, infodynamics_path, octave_path,calc_type)
% Estimation for entropy('est_ent'), total correlation('est_tc'), 
% dual total correlation ('est_dtc') and O-info ('est_oinfo'),
% calculated from 'data' (T samples x N dimensionmatrix) 
% using empirical copulas (see data2gaussian method).
% 
% INPUT
% data = T samples x N variables matrix
% infodynamics_path = Location for JIDT "infodynamics.jar"
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

% Avoid reloading infodynamics.jar and matlab wrappers
if isempty(contains(javaclasspath('-dynamic'), 'infodynamics.jar')) || not(contains(javaclasspath('-dynamic'), 'infodynamics.jar'))
	javaaddpath(infodynamics_path)
end

if not(contains(path,octave_path))
	addpath(octave_path);
end

% obtain calculators from JIDT
if strcmp(calc_type, "gaussian")
	calc_ent = javaObject('infodynamics.measures.continuous.gaussian.EntropyCalculatorMultiVariateGaussian');
    calc_tc = javaObject('infodynamics.measures.continuous.gaussian.MultiInfoCalculatorGaussian');
	calc_dtc = javaObject('infodynamics.measures.continuous.gaussian.DualTotalCorrelationCalculatorGaussian');
	%calc_oinfo = javaObject('infodynamics.measures.continuous.gaussian.OInfoCalculatorGaussian');
elseif strcmp(calc_type, "kozachenko_kraskov")
	calc_ent = javaObject('infodynamics.measures.continuous.kozachenko.EntropyCalculatorMultiVariateKozachenko');
    calc_tc = javaObject('infodynamics.measures.continuous.kraskov.MultiInfoCalculatorKraskov1');
	calc_dtc = javaObject('infodynamics.measures.continuous.kraskov.DualTotalCorrelationCalculatorKraskov');
	%calc_oinfo = javaObject('infodynamics.measures.continuous.kraskov.OInfoCalculatorKraskov');
else
    ME = MException('information_metrics_calc_JIDT:calc_typeNotSupported', ...
        'calc_type %s not supported',calc_type);
    throw(ME)
end
% Obtaining dimensions for 'data'

[T,N] = size(data);

variable = octaveToJavaDoubleArray(data(:,:));

% System Entropy Estimation

tic;
calc_ent.initialise(N);
calc_ent.setObservations(variable);
est_ent = calc_ent.computeAverageLocalOfObservations();
entropy_time = toc;

% Marginal Entropy Estimation

%calc_ent.initialise(1);
%ent_marg = 0;

tic;
calc_tc.initialise(N);
calc_tc.setObservations(variable);
est_tc = calc_tc.computeAverageLocalOfObservations();
tc_time = toc;
%for i=1:N
%	variable_marg = octaveToJavaDoubleArray(data(:,i));
%	calc_ent.setObservations(variable_marg);
%	ent_marg = ent_marg + calc_ent.computeAverageLocalOfObservations();
%end

% Total Correlation Estimation
%est_tc = ent_marg - est_ent;

tic
calc_dtc.initialise(N);
%calc_oinfo.initialise(N);

calc_dtc.setObservations(variable);
%calc_oinfo.setObservations(variable);
% Output Calculation

est_dtc = calc_dtc.computeAverageLocalOfObservations(); 
dtc_time = toc;
tic;
est_oinfo = est_tc - est_dtc; %calc_oinfo.computeAverageLocalOfObservations(); 
oinfo_time = toc;