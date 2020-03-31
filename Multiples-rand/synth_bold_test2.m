addpath("../");
load synth_bold_sig_bp_filtered.mat;

Ts = [100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 5990];
Cols = [3 5 7 9];
%Ts = [100 200 300 400];
javapath = '/home/srodriguez/rodrigoJIDT/RodrigoJIDT/infodynamics.jar';
octavepath = '/home/srodriguez/repositorios/jidt/demos/octave';

nts=length(Ts);
nc = length(Cols);
n_runs = 10;
kozachenko_ent = zeros(n_runs,nc,nts);
gaussian_ent = zeros(n_runs,nc,nts);

kraskov_tc = zeros(n_runs,nc,nts);
gaussian_tc = zeros(n_runs,nc,nts);

kraskov_dtc = zeros(n_runs,nc,nts);
gaussian_dtc = zeros(n_runs,nc,nts);

kraskov_oinfo = zeros(n_runs,nc,nts);
gaussian_oinfo = zeros(n_runs,nc,nts);
[~,n_cols] = size(bold_sig_1);
for r=1:n_runs
    for c=1:nc
        for t=1:nts
            T = Ts(t);
            C = Cols(c);
            randcolums = randperm(n_cols,C);
            data = bold_sig_1(1:T,randcolums);
            res = information_metrics(data,javapath,octavepath);

            gaussian_ent(r,c,t) = res('JIDT_ent_gaussian'); %tmp_ent_g;
            gaussian_tc(r,c,t) = res('JIDT_tc_gaussian'); %tmp_tc_g;
            gaussian_dtc(r,c,t) = res('JIDT_dtc_gaussian'); %tmp_dtc_g;
            gaussian_oinfo(r,c,t) = res('JIDT_oinfo_gaussian'); %tmp_oinfo_g;

            kozachenko_ent(r,c,t) = res('JIDT_ent_kozachenko'); %tmp_ent_k;
            kraskov_tc(r,c,t) = res('JIDT_tc_kraskov'); %tmp_tc_k;
            kraskov_dtc(r,c,t) = res('JIDT_dtc_kraskov'); %tmp_dtc_k;
            kraskov_oinfo(r,c,t) = res('JIDT_oinfo_kraskov'); %

            est_ent(r,c,t) = res('est_ent'); %tmp_ent;
            est_tc(r,c,t) = res('est_tc'); %tmp_tc;
            est_dtc(r,c,t) = res('est_dtc'); %tmp_dtc;
            est_oinfo(r,c,t) = res('est_oinfo'); %tmp_oinfo;

        end
        disp(['Done for C=', num2str(C)])
    end
    disp(['=============== Done run = ', num2str(r), ' ================='])
end

cols = hsv(nc);

cleg = arrayfun(@(x) ['C=',num2str(x)],Cols,'uni',0);
pp = cell(nc,1);

figure('Name','TC Estimation','NumberTitle','on');
for c=1:nc
    ax1 = subplot(3,1,1);
    dat = squeeze(est_tc(:,c,:));
    shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    %plot(ax1, Ts, est_tc)
    title(ax1,'TC GC Estimation')
    ylabel(ax1,'Nats')
    %ylim(ax1, [0 0.2])
    %legend(ax1, '3', '4', '5', '6', '7', '8', '9', '10');


    % Bottom plot
    ax2 = subplot(3,1,2);;
    dat = squeeze(gaussian_tc(:,c,:));
    shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    %plot(ax2, Ts, gaussian_tc )
    title(ax2,'TC Gaussian-JIDT Estimation')
    ylabel(ax2,'Nats')
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')


    % Bottom plot
    ax3 = subplot(3,1,3);;
    dat = squeeze(kraskov_tc(:,c,:));
    pp{c}=shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    %plot(ax3, Ts, kraskov_tc )
    title(ax3,'TC Kraskov Estimation')
    ylabel(ax3,'Nats')
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')
    
end
legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
    cleg,'location','best')
    
savefig("figures/TC-5990-filtered.fig")

figure('Name','DTC Estimation','NumberTitle','on');
for c=1:nc
    ax1 = subplot(3,1,1);
    dat = squeeze(est_dtc(:,c,:));
    shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    %plot(ax1, Ts, est_tc)
    title(ax1,'DTC GC Estimation')
    ylabel(ax1,'Nats')
    %ylim(ax1, [0 0.2])
    %legend(ax1, '3', '4', '5', '6', '7', '8', '9', '10');


    % Bottom plot
    ax2 = subplot(3,1,2);;
    dat = squeeze(gaussian_dtc(:,c,:));
    shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    %plot(ax2, Ts, gaussian_dtc )
    title(ax2,'DTC Gaussian-JIDT Estimation')
    ylabel(ax2,'Nats')
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')


    % Bottom plot
    ax3 = subplot(3,1,3);;
    dat = squeeze(kraskov_dtc(:,c,:));
    pp{c}=shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    %plot(ax3, Ts, kraskov_tc )
    title(ax3,'DTC Kraskov Estimation')
    ylabel(ax3,'Nats')
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')
    
end
legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
    cleg,'location','best')
savefig("figures/DTC-5990-filtered.fig")