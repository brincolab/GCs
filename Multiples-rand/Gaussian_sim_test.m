clear;
addpath("../");
%load longer_synth_bold_sig_bp_filtered.mat;

Ts = [100 200 300 400 500 600 700 800 900 1000 1500 ...
      2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 ...
      10000 10500 11000 11500 12000 12500 12990];
Cols = [3 5 7 9]; % 11 13 15 17 19 21 23 25 27 29 31];
%Ts = [100 200 300 400];
javapath = '/home/srodriguez/rodrigoJIDT/RodrigoJIDT/infodynamics.jar';
octavepath = '/home/srodriguez/repositorios/jidt/demos/octave';

nts=length(Ts);
nc = length(Cols);
n_runs = 1;
kozachenko_ent = zeros(n_runs,nc,nts);
gaussian_ent = zeros(n_runs,nc,nts);

teor_tc = zeros(n_runs,nc,nts);
teor_dtc = zeros(n_runs,nc,nts);

est_tc = zeros(n_runs,nc,nts);
kraskov_tc = zeros(n_runs,nc,nts);
gaussian_tc = zeros(n_runs,nc,nts);

est_tc_time = zeros(n_runs,nc,nts);
kraskov_tc_time = zeros(n_runs,nc,nts);
gaussian_tc_time = zeros(n_runs,nc,nts);


est_dtc = zeros(n_runs,nc,nts);
kraskov_dtc = zeros(n_runs,nc,nts);
gaussian_dtc = zeros(n_runs,nc,nts);

est_dtc_time = zeros(n_runs,nc,nts);
temp_times = zeros(n_runs,nc,nts);
kraskov_dtc_time = zeros(n_runs,nc,nts);
gaussian_dtc_time = zeros(n_runs,nc,nts);

kraskov_oinfo = zeros(n_runs,nc,nts);
gaussian_oinfo = zeros(n_runs,nc,nts);

%bold_sig = bold_sig_2.';
%[~,n_cols] = size(bold_sig);
for r=1:1
    for c=1:nc
        C = Cols(c);
        %N = Ns(n);
        mu = zeros(C,1);
         % Covariance matrices
        covmat = (-sqrt(C)) + 2.*sqrt(C).*rand(C);            
        covmat = 0.5.*(covmat + covmat');
        covmat(1:C+1:end) = (C);
        all_min_1=(arrayfun(@(x) setdiff(1:C,x),1:C,'uni',0)');
        all_dat = mvnrnd(mu,covmat,12990);
        for t=1:nts
            
            T = Ts(t);
            data = all_dat(1:T,:);
            
            %data = bold_sig(1:T,randcolums);
            res = information_metrics(data,javapath,octavepath);
            
            %Teoricos
            detmv = log(det(covmat));
            detmv_min_1=log(cellfun(@(x) det(covmat(x,x)),all_min_1));
            single_vars = diag(covmat);
            
            teor_var_ents = 0.5.*log(((2*pi*exp(1))).*single_vars);
            teor_ent = 0.5.*(log(((2*pi*exp(1)).^C)) + detmv);
            teor_ent_min_one= 0.5.*(log((2*pi*exp(1)).^(C-1)) + detmv_min_1);
            
            teor_tc(r,c,t) = sum(teor_var_ents) - teor_ent;
            teor_dtc(r,c,t) = sum(teor_ent_min_one) - (C-1).*teor_ent;
            
            %Calculados
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
            
            gaussian_tc_time(r,c,t) = res('JIDT_tc_gaussian_time'); %tmp_tc_g;
            gaussian_dtc_time(r,c,t) = res('JIDT_dtc_gaussian_time'); %tmp_dtc_g;
            
            kraskov_tc_time(r,c,t) = res('JIDT_tc_kraskov_time'); %tmp_tc_k;
            kraskov_dtc_time(r,c,t) = res('JIDT_dtc_kraskov_time'); %tmp_dtc_k;
            
            est_tc_time(r,c,t) = res('est_tc_time'); %tmp_tc;
            est_dtc_time(r,c,t) = res('est_dtc_time'); %tmp_dtc;
            temp_times(r,c,t) = res('temp_time'); %tmp_dtc;
        end
        disp(['Done for C=', num2str(C)])
    end
    disp(['=============== Done run = ', num2str(r), ' ================='])
end

cols = hsv(nc);

cleg = arrayfun(@(x) ['C=',num2str(x)],Cols,'uni',0);
pp = cell(nc,1);

figure('Name','TC Estimation','NumberTitle','on');
ax1 = subplot(4,1,1);
ax2 = subplot(4,1,2);
ax3 = subplot(4,1,3);
ax4 = subplot(4,1,4);
title(ax1,'TC Teorico');
ylabel(ax1,'Nats');
title(ax2,'TC GC Estimation');
ylabel(ax2,'Nats');
title(ax3,'TC Gaussian-JIDT Estimation');
ylabel(ax3,'Nats');
title(ax4,'TC Kraskov Estimation')
ylabel(ax4,'Nats')
%hold on;
hold (ax1,'on');
hold (ax2,'on');
hold (ax3,'on');
hold (ax4,'on');
for c=1:nc
    dat = squeeze(teor_tc(:,c,:));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax1, Ts,dat)

    dat = squeeze(est_tc(:,c,:));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax2, Ts,dat)
    
    
    %ylim(ax1, [0 0.2])
    %;


    % Bottom plot
    
    dat = squeeze(gaussian_tc(:,c,:));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax3, Ts, dat )
    
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')


    % Bottom plot
    
    dat = squeeze(kraskov_tc(:,c,:));
    %pp{c}=shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax4, Ts, dat )
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')
    
end
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31')
hold (ax1,'off');
hold (ax2,'off');
hold (ax3,'off');
hold (ax4,'off');

%legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
%    cleg,'location','best')
savefig("figures/TC-Gaussian-sim.fig")

figure('Name','TC Estimation Zoom','NumberTitle','on');
ax1 = subplot(4,1,1);
ax2 = subplot(4,1,2);
ax3 = subplot(4,1,3);
ax4 = subplot(4,1,4);
title(ax1,'TC Teorico');
ylabel(ax1,'Nats');
title(ax2,'TC GC Estimation');
ylabel(ax2,'Nats');
title(ax3,'TC Gaussian-JIDT Estimation');
ylabel(ax3,'Nats');
title(ax4,'TC Kraskov Estimation')
ylabel(ax4,'Nats')


%hold on;
hold (ax1,'on');
hold (ax2,'on');
hold (ax3,'on');
hold (ax4,'on');
for c=1:nc
    dat = squeeze(teor_tc(:,c,:));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax1, Ts,dat)
    dat = squeeze(est_tc(:,c,1:10));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax2, Ts(1:10),dat)
    
    
    %ylim(ax1, [0 0.2])
    %;


    % Bottom plot
    
    dat = squeeze(gaussian_tc(:,c,1:10));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax3, Ts(1:10), dat )
    
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')


    % Bottom plot
    
    dat = squeeze(kraskov_tc(:,c,1:10));
    %pp{c}=shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax4, Ts(1:10), dat )
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')
    
end
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31')
hold (ax1,'off');
hold (ax2,'off');
hold (ax3,'off');
hold (ax4,'off');

%legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
%    cleg,'location','best')
savefig("figures/TC-Gaussian-sim-zoom.fig")

figure('Name','DTC Estimation','NumberTitle','on');
ax1 = subplot(4,1,1);
ax2 = subplot(4,1,2);;
ax3 = subplot(4,1,3);;
ax4 = subplot(4,1,4);;

title(ax1,'DTC Teorico')
title(ax2,'DTC GC Estimation')
title(ax3,'DTC Gaussian-JIDT Estimation')
title(ax4,'DTC Kraskov Estimation')

ylabel(ax1,'Nats')
ylabel(ax2,'Nats')
ylabel(ax3,'Nats')
ylabel(ax4,'Nats')

hold (ax1,'on');
hold (ax2,'on');
hold (ax3,'on');
hold (ax4,'on');
for c=1:nc
    dat = squeeze(teor_dtc(:,c,:));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax1, Ts, dat)
    
    dat = squeeze(est_dtc(:,c,:));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax2, Ts, dat)
    
    
    %ylim(ax1, [0 0.2])
    %legend(ax1, '3', '4', '5', '6', '7', '8', '9', '10');


    % Bottom plot
    
    dat = squeeze(gaussian_dtc(:,c,:));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax3, Ts, dat )
    
    
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')


    % Bottom plot
    
    dat = squeeze(kraskov_dtc(:,c,:));
    %pp{c}=shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax4, Ts, dat )
    
    %ylim(ax4, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')
    
end
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31')
hold (ax1,'off');
hold (ax2,'off');
hold (ax3,'off');
hold (ax4,'off');



%legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
%    cleg,'location','best')
savefig("figures/DTC-Gaussian-sim.fig")

figure('Name','DTC Estimation Zoom','NumberTitle','on');
ax1 = subplot(4,1,1);
ax2 = subplot(4,1,2);;
ax3 = subplot(4,1,3);;
ax4 = subplot(4,1,4);;

title(ax1,'DTC Teorico')
title(ax2,'DTC GC Estimation')
title(ax3,'DTC Gaussian-JIDT Estimation')
title(ax4,'DTC Kraskov Estimation')

ylabel(ax1,'Nats')
ylabel(ax2,'Nats')
ylabel(ax3,'Nats')
ylabel(ax4,'Nats')
hold (ax1,'on');
hold (ax2,'on');
hold (ax3,'on');
hold (ax4,'on');
for c=1:nc
    
    dat = squeeze(teor_dtc(:,c,1:10));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax1, Ts(1:10), dat)
    
    dat = squeeze(est_dtc(:,c,1:10));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax2, Ts(1:10), dat)
    
    
    %ylim(ax1, [0 0.2])
    %legend(ax1, '3', '4', '5', '6', '7', '8', '9', '10');


    % Bottom plot
    
    dat = squeeze(gaussian_dtc(:,c,1:10));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax3, Ts(1:10), dat )
    
    
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')


    % Bottom plot
    
    dat = squeeze(kraskov_dtc(:,c,1:10));
    %pp{c}=shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax4, Ts(1:10), dat )
    
    %ylim(ax2, [0 0.2])
    %legend(ax4, 'TC', 'DTC', 'O-Info')
    
end
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31')
hold (ax1,'off');
hold (ax2,'off');
hold (ax3,'off');
hold (ax4,'off');



%legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
%    cleg,'location','best')
savefig("figures/DTC-Gaussian-sim-zoom.fig")

figure('Name','Time for TC Estimation','NumberTitle','on');
ax1 = subplot(4,1,1);
ax2 = subplot(4,1,2);
ax3 = subplot(4,1,3);
ax4 = subplot(4,1,4);
title(ax1,'Time TC GC Estimation');
ylabel(ax1,'Seconds');
title(ax2,'Entropies calc)');
ylabel(ax2,'Seconds');
title(ax3,'Time TC Gaussian-JIDT Estimation');
ylabel(ax3,'Seconds');
title(ax4,'Time TC Kraskov Estimation')
ylabel(ax4,'Seconds')

%hold on;
hold (ax1,'on');
hold (ax2,'on');
hold (ax3,'on');
hold (ax4,'on');
for c=1:nc
    dat = squeeze(est_tc_time(:,c,:));
    plot(ax1, Ts,dat);
    
    dat = squeeze(temp_times(:,c,:));
    plot(ax2, Ts,dat);

    % Bottom plot
    
    dat = squeeze(gaussian_tc_time(:,c,:));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax3, Ts, dat )
    
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')


    % Bottom plot
    
    dat = squeeze(kraskov_tc_time(:,c,:));
    %pp{c}=shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax4, Ts, dat )
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')
    
end
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31')
hold (ax1,'off');
hold (ax2,'off');
hold (ax3,'off');
hold (ax4,'off');

%legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
%    cleg,'location','best')
savefig("figures/TC-12990-filtered_time.fig")

figure('Name','DTC Estimation','NumberTitle','on');
ax1 = subplot(4,1,1);
ax2 = subplot(4,1,2);;
ax3 = subplot(4,1,3);;
ax4 = subplot(4,1,4);;

title(ax1,'Time DTC GC Estimation')
title(ax2,'Entropies calc')
title(ax3,'Time DTC Gaussian-JIDT Estimation')
title(ax4,'Time DTC Kraskov Estimation')

ylabel(ax1,'Seconds')
ylabel(ax2,'Seconds')
ylabel(ax3,'Seconds')
ylabel(ax4,'Seconds')
hold (ax1,'on');
hold (ax2,'on');
hold (ax3,'on');
hold (ax4,'on');
for c=1:nc
    
    dat = squeeze(est_dtc_time(:,c,:));
    plot(ax1, Ts, dat)
    
    dat = squeeze(temp_times(:,c,:));
    plot(ax2, Ts,dat);
    
    % Bottom plot
    
    dat = squeeze(gaussian_dtc_time(:,c,:));    
    plot(ax3, Ts, dat )

    % Bottom plot
    
    dat = squeeze(kraskov_dtc_time(:,c,:));
    %pp{c}=shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax4, Ts, dat )
    
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')
    
end
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31')
hold (ax1,'off');
hold (ax2,'off');
hold (ax3,'off');
hold (ax4,'off');


%legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
%    cleg,'location','best')
savefig("figures/DTC-12990-filtered_time.fig")