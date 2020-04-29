clear;
addpath("../");
%load longer_synth_bold_sig_bp_filtered.mat;

%Ts = 0:1000:60000;

Ts = [10,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,15000,20000,30000,40000,50000,59990];
%Ts(1) = 100;
Cols = [3];

load bold_60000.mat
bold = bold(:,11:60000);
javapath = '/home/srodriguez/rodrigoJIDT/RodrigoJIDT/infodynamics.jar';
octavepath = '/home/srodriguez/repositorios/jidt/demos/octave';

nts=length(Ts);
nc = length(Cols);
n_runs = 1;
kozachenko_ent = zeros(n_runs,nc,nts);
gaussian_ent = zeros(n_runs,nc,nts);

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

nsurr = 250;
surrogates_gc_tc = zeros(n_runs,nc,nts,nsurr);
surrogates_jg_tc = zeros(n_runs,nc,nts,nsurr);
surrogates_jk_tc = zeros(n_runs,nc,nts,nsurr);

surrogates_gc_dtc = zeros(n_runs,nc,nts,nsurr);
surrogates_jg_dtc = zeros(n_runs,nc,nts,nsurr);
surrogates_jk_dtc = zeros(n_runs,nc,nts,nsurr);

bold = bold';
[~,n_cols] = size(bold);
C_s = Cols(1);
randcolums = randperm(n_cols,C_s);
save randcolumns.mat randcolums

%load randcolumns.mat
full_data = bold(:,randcolums);

for r=1:1
    for c=1:nc
        for t=1:nts
            T = Ts(t);
            data = full_data(1:T,:);
            shuffled_dat = data;
            disp(['Surrogates for C=', num2str(C_s),', T=', num2str(T)])
            for i=1:nsurr
                for j=1:size(shuffled_dat,2)
                    shuffled_dat(:,j) = shuffled_dat(randperm(size(shuffled_dat, 1)), j);
                end
                res_surr = information_metrics(shuffled_dat,javapath,octavepath);
                surrogates_gc_tc(r,c,t,i) = res_surr('est_tc');
                surrogates_jg_tc(r,c,t,i) = res_surr('JIDT_tc_gaussian');
                surrogates_jk_tc(r,c,t,i) = res_surr('JIDT_tc_kraskov');
                surrogates_gc_dtc(r,c,t,i) = res_surr('est_dtc');
                surrogates_jg_dtc(r,c,t,i) = res_surr('JIDT_dtc_gaussian');
                surrogates_jk_dtc(r,c,t,i) = res_surr('JIDT_dtc_kraskov');
                
            end

            %for i=1:nsurr
            %    for j=1:size(shuffled_dat,2)
            %        shuffled_dat(:,j) = shuffled_dat(randperm(size(shuffled_dat, 1)), j);
            %    end
            %    res_surr = information_metrics(shuffled_dat,javapath,octavepath);
            %    surrogates_gc_dtc(r,c,t,i) = res_surr('est_dtc');
            %    surrogates_jg_dtc(r,c,t,i) = res_surr('JIDT_dtc_gaussian');
            %    surrogates_jk_dtc(r,c,t,i) = res_surr('JIDT_dtc_kraskov');
            %end
            disp(['Surrogates Done'])
            res = information_metrics(data,javapath,octavepath);
            
            gaussian_ent(r,c,t) = res('JIDT_ent_gaussian'); %tmp_ent_g;
            gaussian_tc(r,c,t) = res('JIDT_tc_gaussian') - mean(surrogates_jg_tc(r,c,t,:)); %tmp_tc_g;
            gaussian_dtc(r,c,t) = res('JIDT_dtc_gaussian') - mean(surrogates_jg_dtc(r,c,t,:)); %tmp_dtc_g;
            gaussian_oinfo(r,c,t) = res('JIDT_oinfo_gaussian'); %tmp_oinfo_g;

            kozachenko_ent(r,c,t) = res('JIDT_ent_kozachenko'); %tmp_ent_k;
            kraskov_tc(r,c,t) = res('JIDT_tc_kraskov') - mean(surrogates_jk_tc(r,c,t,:)); %tmp_tc_k;
            kraskov_dtc(r,c,t) = res('JIDT_dtc_kraskov') - mean(surrogates_jk_dtc(r,c,:)); %tmp_dtc_k;
            kraskov_oinfo(r,c,t) = res('JIDT_oinfo_kraskov'); %

            est_ent(r,c,t) = res('est_ent'); %tmp_ent;
            est_tc(r,c,t) = res('est_tc') - mean(surrogates_gc_tc(r,c,t,:)); %tmp_tc;
            est_dtc(r,c,t) = res('est_dtc') - mean(surrogates_gc_dtc(r,c,t,:)); %tmp_dtc;
            est_oinfo(r,c,t) = res('est_oinfo'); %tmp_oinfo;
            
            gaussian_tc_time(r,c,t) = res('JIDT_tc_gaussian_time'); %tmp_tc_g;
            gaussian_dtc_time(r,c,t) = res('JIDT_dtc_gaussian_time'); %tmp_dtc_g;
            
            kraskov_tc_time(r,c,t) = res('JIDT_tc_kraskov_time'); %tmp_tc_k;synth_bold
            kraskov_dtc_time(r,c,t) = res('JIDT_dtc_kraskov_time'); %tmp_dtc_k;
            
            est_tc_time(r,c,t) = res('est_tc_time'); %tmp_tc;
            est_dtc_time(r,c,t) = res('est_dtc_time'); %tmp_dtc;
            temp_times(r,c,t) = res('temp_time'); %tmp_dtc;
            save results/raw/synth_bold_60k_calcs.mat Ts gaussian_tc gaussian_dtc kraskov_tc kraskov_dtc est_tc est_dtc
            save results/raw/synth_bold_60k_times.mat gaussian_tc_time gaussian_dtc_time kraskov_tc_time kraskov_dtc_time est_tc_time est_dtc_time temp_times
            save results/raw/synth_bold_60k_surrogates.mat surrogates_jg_tc surrogates_jg_dtc surrogates_jk_tc surrogates_jk_dtc surrogates_gc_tc surrogates_gc_dtc
        end
        disp(['Done for C=', num2str(c)])
    end
    disp(['=============== Done run = ', num2str(r), ' ================='])
end

cols = hsv(nc);

cleg = arrayfun(@(x) ['C=',num2str(x)],Cols,'uni',0);
pp = cell(nc,1);

figure('Name','TC Estimation','NumberTitle','on');
ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);
ax3 = subplot(3,1,3);
title(ax1,'TC GC Estimation');
ylabel(ax1,'Nats');
title(ax2,'TC Gaussian-JIDT Estimation');
ylabel(ax2,'Nats');
title(ax3,'TC Kraskov Estimation')
ylabel(ax3,'Nats')
%hold on;
hold (ax1,'on');
hold (ax2,'on');
hold (ax3,'on');
for c=1:nc
    dat = squeeze(est_tc(:,c,:));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax1, Ts,dat)
    
    
    %ylim(ax1, [0 0.2])
    %;


    % Bottom plot
    
    dat = squeeze(gaussian_tc(:,c,:));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax2, Ts, dat )
    
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')


    % Bottom plot
    
    dat = squeeze(kraskov_tc(:,c,:));
    %pp{c}=shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax3, Ts, dat )
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')
    
end
legend(ax2, '3', '5', '7', '9', '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31')
hold (ax1,'off');
hold (ax2,'off');
hold (ax3,'off');

%legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
%    cleg,'location','best')
%savefig("figures/TC-12990-filtered.fig")

figure('Name','TC Estimation Zoom','NumberTitle','on');
ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);
ax3 = subplot(3,1,3);
title(ax1,'TC GC Estimation');
ylabel(ax1,'Nats');
title(ax2,'TC Gaussian-JIDT Estimation');
ylabel(ax2,'Nats');
title(ax3,'TC Kraskov Estimation')
ylabel(ax3,'Nats')
%hold on;
hold (ax1,'on');
hold (ax2,'on');
hold (ax3,'on');
for c=1:nc
    dat = squeeze(est_tc(:,c,1:10));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax1, Ts(1:10),dat)
    
    
    %ylim(ax1, [0 0.2])
    %;


    % Bottom plot
    
    dat = squeeze(gaussian_tc(:,c,1:10));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax2, Ts(1:10), dat )
    
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')


    % Bottom plot
    
    dat = squeeze(kraskov_tc(:,c,1:10));
    %pp{c}=shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax3, Ts(1:10), dat )
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')
    
end
legend(ax2, '3', '5', '7', '9'); %, '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31')
hold (ax1,'off');
hold (ax2,'off');
hold (ax3,'off');

%legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
%    cleg,'location','best')
%savefig("figures/TC-12990-filtered-zoom.fig")

figure('Name','DTC Estimation','NumberTitle','on');
ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);;
ax3 = subplot(3,1,3);;

title(ax1,'DTC GC Estimation')
title(ax2,'DTC Gaussian-JIDT Estimation')
title(ax3,'DTC Kraskov Estimation')

ylabel(ax1,'Nats')
ylabel(ax2,'Nats')
ylabel(ax3,'Nats')
hold (ax1,'on');
hold (ax2,'on');
hold (ax3,'on');
for c=1:nc
    
    dat = squeeze(est_dtc(:,c,:));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax1, Ts, dat)
    
    
    %ylim(ax1, [0 0.2])
    %legend(ax1, '3', '4', '5', '6', '7', '8', '9', '10');


    % Bottom plot
    
    dat = squeeze(gaussian_dtc(:,c,:));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax2, Ts, dat )
    
    
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')


    % Bottom plot
    
    dat = squeeze(kraskov_dtc(:,c,:));
    %pp{c}=shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax3, Ts, dat )
    
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')
    
end
legend(ax2, '3', '5', '7', '9', '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31')
hold (ax1,'off');
hold (ax2,'off');
hold (ax3,'off');



%legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
%    cleg,'location','best')
%savefig("figures/DTC-12990-filtered.fig")

figure('Name','DTC Estimation Zoom','NumberTitle','on');
ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);;
ax3 = subplot(3,1,3);;

title(ax1,'DTC GC Estimation')
title(ax2,'DTC Gaussian-JIDT Estimation')
title(ax3,'DTC Kraskov Estimation')

ylabel(ax1,'Nats')
ylabel(ax2,'Nats')
ylabel(ax3,'Nats')
hold (ax1,'on');
hold (ax2,'on');
hold (ax3,'on');
for c=1:nc
    
    dat = squeeze(est_dtc(:,c,1:10));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax1, Ts(1:10), dat)
    
    
    %ylim(ax1, [0 0.2])
    %legend(ax1, '3', '4', '5', '6', '7', '8', '9', '10');


    % Bottom plot
    
    dat = squeeze(gaussian_dtc(:,c,1:10));
    %shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    
    plot(ax2, Ts(1:10), dat )
    
    
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')


    % Bottom plot
    
    dat = squeeze(kraskov_dtc(:,c,1:10));
    %pp{c}=shadedErrorBar(Ts,mean(dat,1),std(dat,1),'lineprops',{'color',cols(c,:),'linewidth',2});
    plot(ax3, Ts(1:10), dat )
    
    %ylim(ax2, [0 0.2])
    %legend(ax3, 'TC', 'DTC', 'O-Info')
    
end
legend(ax2, '3', '5', '7', '9', '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31')
hold (ax1,'off');
hold (ax2,'off');
hold (ax3,'off');



%legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
%    cleg,'location','best')
%savefig("figures/DTC-12990-filtered-Zoom.fig")

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
legend(ax2, '3', '5', '7', '9', '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31')
hold (ax1,'off');
hold (ax2,'off');
hold (ax3,'off');
hold (ax4,'off');

%legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
%    cleg,'location','best')
%savefig("figures/TC-12990-filtered_time.fig")

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
legend(ax2, '3', '5', '7', '9', '11', '13', '15', '17', '19', '21', '23', '25', '27', '29', '31')
hold (ax1,'off');
hold (ax2,'off');
hold (ax3,'off');
hold (ax4,'off');


%legend([pp{1}.mainLine,pp{2}.mainLine,pp{3}.mainLine,pp{4}.mainLine],...
%    cleg,'location','best')
%savefig("figures/DTC-12990-filtered_time.fig")