% Plot the absolute value between teoretical and estimated information metrics
%

curves = nc;
mean_diffs = zeros(curves,3,34);

for i=1:curves
    dat0 = squeeze(teor_tc(:,i,:));
    dat1 = squeeze(est_tc(:,i,:)); mean_dat1 = mean(dat1,1);
    dat2 = squeeze(gaussian_tc(:,i,:)); mean_dat2 = mean(dat2,1);
    dat3 = squeeze(kraskov_tc(:,i,:)); mean_dat3 = mean(dat3,1);
    mean_diffs(i,1,:) = abs(dat0 - dat1);
    mean_diffs(i,2,:) = abs(dat0 - dat2);
    mean_diffs(i,3,:) = abs(dat0 - dat3);
end


figure('Name','TC Estimation','NumberTitle','on');
ax1 = subplot(3,1,1);
plot(ax1, Ts, squeeze(mean_diffs(:,1,:)))
ylim(ax1, [-0.1 0.5])
title(ax1,'abs((TC Teorico) - (Gaussian Copula TC))')

ax2 = subplot(3,1,2);
plot(ax2, Ts, squeeze(mean_diffs(:,2,:)))
ylim(ax2, [-0.1 0.5])
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31');
title(ax2,'abs((TC Teorico)  - (Gaussian JIDT TC))')

ax3 = subplot(3,1,3);
plot(ax3, Ts, squeeze(mean_diffs(:,3,:)))
ylim(ax3, [-0.1 0.5])
title(ax3,'abs((TC Teorico) - (kraskov JIDT TC))')
%%%%
% ZOOM
%%%%
figure('Name','TC Estimation Zoom','NumberTitle','on');
ax1 = subplot(3,1,1);
plot(ax1, Ts(1:10), squeeze(mean_diffs(:,1,1:10)))
ylim(ax1, [-0.1 0.5])
title(ax1,'abs((TC Teorico) - (Gaussian Copula TC))')

ax2 = subplot(3,1,2);
plot(ax2, Ts(1:10), squeeze(mean_diffs(:,2,1:10)))
ylim(ax2, [-0.1 0.5])
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31');
title(ax2,'abs((TC Teorico)  - (Gaussian JIDT TC))')

ax3 = subplot(3,1,3);
plot(ax3, Ts(1:10), squeeze(mean_diffs(:,3,1:10)))
ylim(ax3, [-0.1 0.5])
title(ax3,'abs((TC Teorico) - (kraskov JIDT TC))')

%%%%
% DTC
%%%%
mean_diffs = zeros(curves,3,34);

for i=1:curves
    dat0 = squeeze(teor_dtc(:,i,:));
    dat1 = squeeze(est_dtc(:,i,:)); mean_dat1 = mean(dat1,1);
    dat2 = squeeze(gaussian_dtc(:,i,:)); mean_dat2 = mean(dat2,1);
    dat3 = squeeze(kraskov_dtc(:,i,:)); mean_dat3 = mean(dat3,1);
    mean_diffs(i,1,:) = abs(dat0 - dat1);
    mean_diffs(i,2,:) = abs(dat0 - dat2);
    mean_diffs(i,3,:) = abs(dat0 - dat3);
end


figure('Name','DTC Estimation','NumberTitle','on');
ax1 = subplot(3,1,1);
plot(ax1, Ts, squeeze(mean_diffs(:,1,:)))
ylim(ax1, [-0.1 0.5])
title(ax1,'abs((DTC Teorico) - (Gaussian Copula DTC))')

ax2 = subplot(3,1,2);
plot(ax2, Ts, squeeze(mean_diffs(:,2,:)))
ylim(ax2, [-0.1 0.5])
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31');
title(ax2,'abs((DTC Teorico)  - (Gaussian JIDT DTC))')

ax3 = subplot(3,1,3);
plot(ax3, Ts, squeeze(mean_diffs(:,3,:)))
ylim(ax3, [-0.1 0.5])
title(ax3,'abs((DTC Teorico) - (kraskov JIDT DTC))')

figure('Name','DTC Estimation Zoom','NumberTitle','on');
ax1 = subplot(3,1,1);
plot(ax1, Ts(1:10), squeeze(mean_diffs(:,1,1:10)))
ylim(ax1, [-0.1 0.5])
title(ax1,'abs((DTC Teorico) - (Gaussian Copula DTC))')

ax2 = subplot(3,1,2);
plot(ax2, Ts(1:10), squeeze(mean_diffs(:,2,1:10)))
ylim(ax2, [-0.1 0.5])
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31');
title(ax2,'abs((DTC Teorico)  - (Gaussian JIDT DTC))')

ax3 = subplot(3,1,3);
plot(ax3, Ts(1:10), squeeze(mean_diffs(:,3,1:10)))
ylim(ax3, [-0.1 0.5])
title(ax3,'abs((DTC Teorico) - (kraskov JIDT DTC))')