%%%%%%%
%
% Parsing matlab to latex table
%
%%%%%%%
disp(['ESTIMACION  GC TC ============'])

for i=1:4
    dat = squeeze(est_tc(:,i,:)); mean_dat = dat; %std_dat = std(dat,1);
    sub_mean = mean_dat([5 12 20 34]);
    %sub_std = std_dat([5 12 20 34]);
    disp(['$ ', num2str((i * 2) +1), ' $ & $', num2str(sub_mean(1)), '$', ' & ', ...
    '$' , num2str(sub_mean(2)), '$', ' & ', ...
    '$' ,num2str(sub_mean(3)), '$', ' & ', ...
    '$' ,num2str(sub_mean(4)), '$',' \\ '])
end

disp(['ESTIMACION  JIDT GAUSSIAN TC ============'])
for i=1:4
    dat = squeeze(gaussian_tc(:,i,:)); mean_dat = dat; %std_dat = std(dat,1);
    sub_mean = mean_dat([5 12 20 34]);
    %sub_std = std_dat([5 12 20 34]);
    disp(['$ ', num2str((i * 2) +1), ' $ & $', num2str(sub_mean(1)), '$', ' & ', ...
    '$' , num2str(sub_mean(2)), '$', ' & ', ...
    '$' ,num2str(sub_mean(3)), '$', ' & ', ...
    '$' ,num2str(sub_mean(4)), '$',' \\ '])
end

disp(['ESTIMACION  JIDT kraskov TC ============'])
for i=1:4
    dat = squeeze(kraskov_tc(:,i,:)); mean_dat = dat; %std_dat = std(dat,1);
    sub_mean = mean_dat([5 12 20 34]);
    %sub_std = std_dat([5 12 20 34]);
    disp(['$ ', num2str((i * 2) +1), ' $ & $', num2str(sub_mean(1)), '$', ' & ', ...
    '$' , num2str(sub_mean(2)), '$', ' & ', ...
    '$' ,num2str(sub_mean(3)), '$', ' & ', ...
    '$' ,num2str(sub_mean(4)), '$',' \\ '])
end

disp(['ESTIMACION GC DTC ============'])
for i=1:4
    dat = squeeze(est_dtc(:,i,:)); mean_dat = dat; % mean(dat,1); std_dat = std(dat,1);
    sub_mean = mean_dat([5 12 20 34]);
    %sub_std = std_dat([5 12 20 34]);
    disp(['$ ' ,num2str((i * 2) +1), ' $ & $', num2str(sub_mean(1)), '$', ' & ', ...
    '$', num2str(sub_mean(2)), '$', ' & ', ...
    '$' ,num2str(sub_mean(3)), '$', ' & ', ...
    '$' ,num2str(sub_mean(4)), '$',' \\ '])
end

disp(['ESTIMACION JIDT GAUSSIAN DTC ============'])
for i=1:4
    dat = squeeze(gaussian_dtc(:,i,:)); mean_dat = dat; %mean(dat,1); std_dat = std(dat,1);
    sub_mean = mean_dat([5 12 20 34]);
    %sub_std = std_dat([5 12 20 34]);
    disp(['$ ' ,num2str((i * 2) +1), ' $ & $', num2str(sub_mean(1)), '$', ' & ', ...
    '$', num2str(sub_mean(2)), '$', ' & ', ...
    '$' ,num2str(sub_mean(3)), '$', ' & ', ...
    '$' ,num2str(sub_mean(4)), '$',' \\ '])
end

disp(['ESTIMACION JIDT kraskov DTC ============'])
for i=1:4
    dat = squeeze(kraskov_dtc(:,i,:)); mean_dat = dat; %mean(dat,1); std_dat = std(dat,1);
    sub_mean = mean_dat([5 12 20 34]);
    %sub_std = std_dat([5 12 20 34]);
    disp(['$ ' ,num2str((i * 2) +1), ' $ & $', num2str(sub_mean(1)), '$', ' & ', ...
    '$', num2str(sub_mean(2)), '$', ' & ', ...
    '$' ,num2str(sub_mean(3)), '$', ' & ', ...
    '$' ,num2str(sub_mean(4)), '$',' \\ '])
end