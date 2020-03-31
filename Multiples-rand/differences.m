% Plot the absolute value between the different estimated information metrics

curves = nc;
mean_diffs = zeros(curves,3,34);

for i=1:curves
    dat1 = squeeze(est_tc(:,i,:)); mean_dat1 = mean(dat1,1);
    dat2 = squeeze(gaussian_tc(:,i,:)); mean_dat2 = mean(dat2,1);
    dat3 = squeeze(kraskov_tc(:,i,:)); mean_dat3 = mean(dat3,1);
    mean_diffs(i,1,:) = abs(dat1 - dat2);
    mean_diffs(i,2,:) = abs(dat1 - dat3);
    mean_diffs(i,3,:) = abs(dat3 - dat2);
end


figure('Name','TC Estimation','NumberTitle','on');
ax1 = subplot(3,1,1);
plot(ax1, Ts, squeeze(mean_diffs(:,1,:)))
ylim(ax1, [-0.1 0.5])
title(ax1,'abs(mean(Gaussian Copula TC) - mean(Gaussian JIDT TC))')

ax2 = subplot(3,1,2);
plot(ax2, Ts, squeeze(mean_diffs(:,2,:)))
ylim(ax2, [-0.1 0.5])
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31');
title(ax2,'abs(mean(Gaussian Copula TC) - mean(kraskov JIDT TC))')

ax3 = subplot(3,1,3);
plot(ax3, Ts, squeeze(mean_diffs(:,3,:)))
ylim(ax3, [-0.1 0.5])
title(ax3,'abs(mean(Gaussian JIDT TC) - mean(kraskov JIDT TC))')


mean_diffs = zeros(curves,3,34);

for i=1:curves
    dat1 = squeeze(est_dtc(:,i,:)); mean_dat1 = mean(dat1,1);
    dat2 = squeeze(gaussian_dtc(:,i,:)); mean_dat2 = mean(dat2,1);
    dat3 = squeeze(kraskov_dtc(:,i,:)); mean_dat3 = mean(dat3,1);
    mean_diffs(i,1,:) = abs(dat1 - dat2);
    
    mean_diffs(i,2,:) = abs(dat1 - dat3);
    mean_diffs(i,3,:) = abs(dat3 - dat2);
end


figure('Name','DTC Estimation','NumberTitle','on');
ax1 = subplot(3,1,1);
plot(ax1, Ts, squeeze(mean_diffs(:,1,:)))
ylim(ax1, [-0.1 0.5])
title(ax1,'abs(mean(Gaussian Copula DTC) - mean(Gaussian JIDT DTC))')

ax2 = subplot(3,1,2);
plot(ax2, Ts, squeeze(mean_diffs(:,2,:)))
ylim(ax2, [-0.1 0.5])
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31');
title(ax2,'abs(mean(Gaussian Copula DTC) - mean(kraskov JIDT DTC))')

ax3 = subplot(3,1,3);
plot(ax3, Ts, squeeze(mean_diffs(:,3,:)))
ylim(ax3, [-0.1 0.5])
title(ax3,'abs(mean(Gaussian JIDT DTC) - mean(kraskov JIDT DTC))')



figure('Name','TC Estimation Zoom','NumberTitle','on');
ax1 = subplot(3,1,1);
plot(ax1, Ts(1:10), squeeze(mean_diffs(:,1,1:10)))
ylim(ax1, [-0.1 0.5])
title(ax1,'abs(mean(Gaussian Copula TC) - mean(Gaussian JIDT TC))')

ax2 = subplot(3,1,2);
plot(ax2, Ts(1:10), squeeze(mean_diffs(:,2,1:10)))
ylim(ax2, [-0.1 0.5])
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31');
title(ax2,'abs(mean(Gaussian Copula TC) - mean(kraskov JIDT TC))')

ax3 = subplot(3,1,3);
plot(ax3, Ts(1:10), squeeze(mean_diffs(:,3,1:10)))
ylim(ax3, [-0.1 0.5])
title(ax3,'abs(mean(Gaussian JIDT TC) - mean(kraskov JIDT TC))')

figure('Name','DTC Estimation Zoom','NumberTitle','on');
ax1 = subplot(3,1,1);
plot(ax1, Ts(1:10), squeeze(mean_diffs(:,1,1:10)))
ylim(ax1, [-0.1 0.5])
title(ax1,'abs(mean(Gaussian Copula DTC) - mean(Gaussian JIDT DTC))')

ax2 = subplot(3,1,2);
plot(ax2, Ts(1:10), squeeze(mean_diffs(:,2,1:10)))
ylim(ax2, [-0.1 0.5])
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31');
title(ax2,'abs(mean(Gaussian Copula DTC) - mean(kraskov JIDT DTC))')

ax3 = subplot(3,1,3);
plot(ax3, Ts(1:10), squeeze(mean_diffs(:,3,1:10)))
ylim(ax3, [-0.1 0.5])
title(ax3,'abs(mean(Gaussian JIDT DTC) - mean(kraskov JIDT DTC))')