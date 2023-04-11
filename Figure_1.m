close all
clear all
clc


load('L:\nttk-data2\palva\Madrid\Full Cohort\Figures Codes\Standalone figure codes\DFA_fEI_estimates.mat')
for k=1:3
    for i=1:size(temp_exp{k},3)
        avg_exp_subj{k}(i,:) = nanmean(temp_exp{k}(:,:,i));
    end
    for i=1:size(temp_fEI{k},3)
        avg_fEI_subj{k}(i,:) = nanmean(temp_fEI{k}(:,:,i));
    end
end


figure
subplot(211)
color_code = {[0.3 0.44 1],[0 0.5 0],[1 0.44 0]};
color_code1 = {[1 0 1],[0 0.5 0.5],[1 0.44 0.5]};
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
p = 95;

%%%%%%%%%%%%%%%%%%%  Plot average DFA exponents %%%%%%%%%%%%%%%%%%%%%%%%

semilogx(MorletBank,mean(avg_exp_subj{1}),'color',color_code{1},'linewidth',2)
hold on
semilogx(MorletBank,mean(avg_exp_subj{2}),'color',color_code{2},'linewidth',2)
semilogx(MorletBank,mean(avg_exp_subj{3}),'color',color_code{3},'linewidth',2)
for i=1:3
    hold on
    for j=1:size(avg_exp_subj{1, 1},2)
        surgt_means = [];
        for t=1:10000
            surgt_idx = [];
            surgt_data =[];
            surgt_idx = randi(size(avg_exp_subj{i},1),size(avg_exp_subj{i},1),1);
            surgt_data = avg_exp_subj{i}(surgt_idx,j);
            surgt_means(t) = mean(surgt_data);
        end
        CI = CIFcn(surgt_means,p);
        Conf_interval{i}(1,j) = CI(1);
        Conf_interval{i}(2,j) = CI(2);
    end
    curve1 = Conf_interval{i}(1,:);
    curve2 = Conf_interval{i}(2,:);
    t = MorletBank;
    time = [t, fliplr(t)];
    inBetween = [curve1, fliplr(curve2)];
    pl=patch(time, inBetween,color_code{i},'edgecolor','none');
    alpha(0.15)
    set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');   
end
semilogx(MorletBank,mean(avg_exp_subj{1}),'color',color_code{1},'linewidth',2,'HandleVisibility','off')
hold on
semilogx(MorletBank,mean(avg_exp_subj{2}),'color',color_code{2},'linewidth',2,'HandleVisibility','off')
semilogx(MorletBank,mean(avg_exp_subj{3}),'color',color_code{3},'linewidth',2,'HandleVisibility','off')

%%%%%%%%%%%%%%%%%%%  Panel b: Between groups stats %%%%%%%%%%%%%%%%%%%%%%%%

grouplabel = [repmat({'cntrl'},size(avg_exp_subj{1},1),1); repmat({'scd'},size(avg_exp_subj{2},1),1) ...
    ; repmat({'mci'},size(avg_exp_subj{3},1),1)];
tbl_true={};
sub_true_distribution =[];
sub_ovrall_sig_05 = zeros(1,size(avg_exp_subj{1},2));
sub_p_true = [];
exp_per_sub = [];
for k=1:3
    temp_hold = [];
    temp_hold = avg_exp_subj{k};
    exp_per_sub = [exp_per_sub;  temp_hold];
end

for m=1:size(exp_per_sub,2)
    [sub_p_true(m), tbl_true{m}] = kruskalwallis(exp_per_sub(:,m),grouplabel,'off');
    sub_true_distribution(m) = tbl_true{m}{2,5};                   
end

%%%%%%%%%%%% FDR correction %%%%%%%%%%%%%%%%%%%%%%%

    [rank rnk_num]= sort(sub_p_true);
    adj_pval = ([1:size(rank,2)]./size(rank,2)).*0.1;
    k = find(adj_pval <=0.05,1,'last')
    rank(k:end) = 1;
    ylim([0.5 0.8])
    ytick = yticks;
    xloc = find(sub_p_true <= 0.05);
    
    adj_p_align(rnk_num) = adj_pval;
    adj_p_align(adj_p_align >= 0.05) = 1;
    round(adj_p_align',4)
    
    if ~isempty(xloc)
        plot(MorletBank(xloc),ytick(end)-0.02,'ks', 'MarkerSize',7);
        adj_p_align(rnk_num) = rank;
        xloc = find(adj_p_align <= 0.05);
        plot(MorletBank(xloc),ytick(end)-0.02,'rd', 'MarkerFaceColor', 'r','MarkerSize',5);
        
        st = xloc(1)    
        lst = find(adj_p_align <=0.05,1,'last')
    else
        st = 1;    
        lst = size(MorletBank,2);
    end

    legend('NC','SCD','MCI')%,'Bonferroni','permutation')
    legend boxoff
    xlim([2.2 90])
    xticks([3 5 10 20 30 50 90]);
    xticklabels([3 5 10 20 30 50 90])
    box off
    ylabel('DFA exponent')

    
%%%%%%%%%%%%%%%%%%%%%%%% Posthoc analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
diff_sub_exp ={};
for k=2:3
    for i=1:size(avg_exp_subj{k})
        diff_sub_exp{k}(i,:) = avg_exp_subj{k}(i,:)-nanmean(avg_exp_subj{1});
    end
end
% figure
subplot(212)
semilogx(MorletBank,mean(diff_sub_exp{2}),'color',color_code1{2},'linewidth',2)
hold on
semilogx(MorletBank,mean(diff_sub_exp{3}),'color',color_code1{3},'linewidth',2)
semilogx(MorletBank,zeros(1,32),'--','color',[0.7 0.7 0.7],'linewidth',1,'HandleVisibility','off')

for i=2:3
    hold on
    for j=1:size(diff_sub_exp{2},2)
        surgt_means = [];
        for t=1:10000
            surgt_idx = [];
            surgt_data =[];
            surgt_idx = randi(size(diff_sub_exp{i},1),size(diff_sub_exp{i},1),1);
            surgt_data = diff_sub_exp{i}(surgt_idx,j);
            surgt_means(t) = mean(surgt_data);
        end
        CI = CIFcn(surgt_means,p);
        Conf_interval{i}(1,j) = CI(1);
        Conf_interval{i}(2,j) = CI(2);
    end
    curve1 = Conf_interval{i}(1,:);
    curve2 = Conf_interval{i}(2,:);
%     t = MorletBank;
    t = MorletBank;
    time = [t, fliplr(t)];
    inBetween = [curve1, fliplr(curve2)];
    pl=patch(time, inBetween,[color_code1{i}],'edgecolor','none');
    alpha(0.15)
    set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
   
end
    
    h = gca; % Get axis to modify
    h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
    h.XAxis.MinorTickValues = [4 6 7 8 9 40 60 70 80];
    grouplabel12 = [ repmat({'cntrl'},size(avg_exp_subj{1},1),1) ...
        ; repmat({'scd'},size(avg_exp_subj{2},1),1)];
    
    sub_true_distribution_dif12=zeros(size(MorletBank)); sub_p_true_dif12=ones(size(MorletBank)); tbl_true_dif12={};
    
    exp_per_sub = [];
    for k=1:2
        temp_hold = [];
        temp_hold = avg_exp_subj{k};
        exp_per_sub = [exp_per_sub;  temp_hold];
    end
    for m=1:size(exp_per_sub,2)
        [sub_p_true_dif12(m), tbl_true_dif12{m}] = kruskalwallis(exp_per_sub(:,m),grouplabel12,'off');
        sub_true_distribution_dif12(m) = tbl_true_dif12{m}{2,5};                   
    end
    
    [rank12 rnk_num12]= sort(sub_p_true_dif12);
    ypos = mean(diff_sub_exp{2});
    xloc = find(sub_p_true_dif12 <= 0.05);
    yloc = ypos(xloc);
    
    adj_pval = ([1:size(rank12,2)]./size(rank12,2)).*0.1;
    k = find(adj_pval <=0.05,1,'last')
    rank12(k+1:end) = 1;    
    for i=1:size(MorletBank,2)
        if rnk_num12(i) >= st && rnk_num12(i) <= lst
            adj_p_align12(rnk_num12(i)) = rank12(i);
        else
            adj_p_align12(rnk_num12(i)) = 1;
        end
    end
    
% % %     adj_p_align12(adj_p_align12 > 0.05)=1;
% % %     xloc = find(adj_p_align12 <= 0.05);
% % %     adj_pval(rnk_num12) = adj_pval;
% % %     adj_p_align12(xloc) = adj_pval(xloc);
    
    yloc = ypos(xloc);
    semilogx(MorletBank(xloc),yloc,'rd', 'MarkerFaceColor', 'r','MarkerSize',5,'HandleVisibility','off'); 
    
    grouplabel13 = [ repmat({'cntrl'},size(avg_exp_subj{1},1),1) ...
        ; repmat({'mci'},size(avg_exp_subj{3},1),1)];    
    sub_true_distribution_dif13=zeros(size(MorletBank)); sub_p_true_dif13=ones(size(MorletBank)); tbl_true_dif13={};
    
    exp_per_sub = [];
    for k=1:2:3
        temp_hold = [];
        temp_hold = avg_exp_subj{k};
        exp_per_sub = [exp_per_sub;  temp_hold];
    end
    for m=1:size(exp_per_sub,2)
        [sub_p_true_dif13(m), tbl_true_dif13{m}] = kruskalwallis(exp_per_sub(:,m),grouplabel13,'off');
        sub_true_distribution_dif13(m) = tbl_true_dif13{m}{2,5};                   
    end

    [rank13 rnk_num13]= sort(sub_p_true_dif13);
    ypos = mean(diff_sub_exp{3});
    xloc = find(sub_p_true_dif13 <= 0.05);
    yloc = ypos(xloc);
    
    adj_pval = ([1:size(rank13,2)]./size(rank13,2)).*0.1
    k = find(adj_pval <=0.05,1,'last')
    rank13(k+1:end) = 1;
    for i=1:size(MorletBank,2)
        if rnk_num13(i) >= st && rnk_num13(i) <= lst
            adj_p_align13(rnk_num13(i)) = rank13(i);
        else
            adj_p_align13(rnk_num13(i)) = 1;
        end
    end
% % %     adj_p_align13(adj_p_align13 > 0.05)=1;
% % %     xloc = find(adj_p_align13 <= 0.05);
% % %     adj_pval(rnk_num13) = adj_pval;
% % %     adj_p_align13(xloc) = adj_pval(xloc);
    
    yloc = ypos(xloc);
    semilogx(MorletBank(xloc),yloc,'rd', 'MarkerFaceColor', 'r','MarkerSize',5,'HandleVisibility','off');
    
    
    diff_sub_exp ={};
    for i=1:size(avg_exp_subj{3})
        diff_sub_exp{1}(i,:) = avg_exp_subj{3}(i,:)-nanmean(avg_exp_subj{2});
    end
    hold on
    semilogx(MorletBank,mean(diff_sub_exp{1}),'color',[1 0 1],'linewidth',2)
    for i=1:1
        hold on
        for j=1:size(diff_sub_exp{1},2)
            surgt_means = [];
            for t=1:10000
                surgt_idx = [];
                surgt_data =[];
                surgt_idx = randi(size(diff_sub_exp{i},1),size(diff_sub_exp{i},1),1);
                surgt_data = diff_sub_exp{i}(surgt_idx,j);
                surgt_means(t) = mean(surgt_data);
            end
            CI = CIFcn(surgt_means,p);
            Conf_interval{i}(1,j) = CI(1);
            Conf_interval{i}(2,j) = CI(2);
        end
        curve1 = Conf_interval{i}(1,:);
        curve2 = Conf_interval{i}(2,:);
%         t = MorletBank;
        t = MorletBank;
        time = [t, fliplr(t)];
        inBetween = [curve1, fliplr(curve2)];
        pl=patch(time, inBetween,[1 0 1],'edgecolor','none');
        alpha(0.10)
        set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

    end
    
    h = gca; % Get axis to modify
    h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
    h.XAxis.MinorTickValues = [4 6 7 8 9 40 60 70 80];
    grouplabe23 = [ repmat({'scd'},size(avg_exp_subj{2},1),1) ...
        ; repmat({'mci'},size(avg_exp_subj{3},1),1)];

    sub_true_distribution_dif23=zeros(size(MorletBank)); sub_p_true_dif23=ones(size(MorletBank)); tbl_true_dif23={};
    
    exp_per_sub = [];
    for k=2:3
        temp_hold = [];
        temp_hold = avg_exp_subj{k};
        exp_per_sub = [exp_per_sub;  temp_hold];
    end
    for m=1:size(exp_per_sub,2)
        [sub_p_true_dif23(m), tbl_true_dif23{m}] = kruskalwallis(exp_per_sub(:,m),grouplabe23,'off');
        sub_true_distribution_dif23(m) = tbl_true_dif23{m}{2,5};                   
    end

    [rank23 rnk_num23]= sort(sub_p_true_dif23);
    ypos = mean(diff_sub_exp{1});
    xloc = find(sub_p_true_dif23 <= 0.05);
    yloc = ypos(xloc);
    
    adj_pval = ([1:size(rank23,2)]./size(rank23,2)).*0.1
    k = find(adj_pval <=0.05,1,'last')
    rank23(k+1:end) = 1;
    adj_p_align23 = sub_p_true_dif23;
    for i=1:size(MorletBank,2)
        if rnk_num23(i) >= st && rnk_num23(i) <= lst
%             adj_p_align23(rnk_num23(i)) = rank23(i);
        else
            adj_p_align23(rnk_num23(i)) = 1;
        end
    end 
    
% % %     adj_p_align23(adj_p_align23 > 0.05)=1;
% % %     xloc = find(adj_p_align23 <= 0.05);
% % %     adj_pval(rnk_num23) = adj_pval;
% % %     adj_p_align23(xloc) = adj_pval(xloc);
    
    yloc = ypos(xloc);
    semilogx(MorletBank(xloc),yloc,'rd', 'MarkerFaceColor', 'r','MarkerSize',5,'HandleVisibility','off');

    xlim([2.2 90])
    xticks([3 5 10 20 30 50 90]);
    xticklabels([3 5 10 20 30 50 90])
    
    xlabel('Frequency (Hz)')
    ylabel('DFA differences')
    legend('SCD-NC','MCI-NC','MCI-SCD')
    legend boxoff
    box off
    

%%%%%%%%%%%%%%%%%% DFA density and tree plots %%%%%%%%%%%%%%%%%%% 


sig_freq=find(adj_p_align23 <= 0.05);
cntrl_sig_freq = avg_exp_subj{1}(:,sig_freq(1:end-1));
SCD_sig_freq = avg_exp_subj{2}(:,sig_freq(1:end-1));
MCI_sig_freq = avg_exp_subj{3}(:,sig_freq(1:end-1));


figure
subplot(131)
yaxis_val = raincloud_yaxis(cntrl_sig_freq(:));
raincloud_plot(round(cntrl_sig_freq(:),2),yaxis_val, color_code{1});
xlim([0.45 1.1])
camroll(90)
x1 = get(gca,'XLim');
yl = get(gca,'YLim');
set(gca,'YLim',[-yl(2)-5 yl(2)+3]);
subplot(132)
yaxis_val = raincloud_yaxis(SCD_sig_freq(:));
raincloud_plot(round(SCD_sig_freq(:),2),yaxis_val, color_code{2});
camroll(90)
set(gca,'XLim',x1);
yl = get(gca,'YLim');
set(gca,'YLim',[-yl(2)-1 yl(2)-0.5]);
subplot(133)
yaxis_val = raincloud_yaxis(MCI_sig_freq(:));
raincloud_plot(round(MCI_sig_freq(:),2),yaxis_val, color_code{3});
camroll(90)
set(gca,'XLim',x1);
yl = get(gca,'YLim');
set(gca,'YLim',[-yl(2)-5 yl(2)+3]);
xlabel('Scaling exponent (\alpha)','interpreter','Tex','FontSize',11)

diff_scd_sig = SCD_sig_freq(:) - mean(cntrl_sig_freq(:));
diff_mci_sig = MCI_sig_freq(:) - mean(cntrl_sig_freq(:));
diff_mciscd_sig = MCI_sig_freq(:) - mean(SCD_sig_freq(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(132)
raincloud_density(round(SCD_sig_freq(:),2),color_code{2});
plot(median(round(SCD_sig_freq(:),2)),-1,'.k','MarkerSize',30)
%%%%%%%%%
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
p = 95;

line([median(round(SCD_sig_freq(:),2))+std(round(SCD_sig_freq(:),2)) median(round(SCD_sig_freq(:),2))],[-1,-1],'color','k','linewidth',2)
line([median(round(SCD_sig_freq(:),2))-std(round(SCD_sig_freq(:),2)) median(round(SCD_sig_freq(:),2))],[-1,-1],'color','k','linewidth',2)
%%%%%%%%%            
ylim([-10 10])
xlim([0.45 1.1])
camroll(90)
x1 = get(gca,'XLim');
box off


subplot(133)
h = raincloud_density(round(MCI_sig_freq(:),2),color_code{3});
plot(median(round(MCI_sig_freq(:),2)),-1,'.k','MarkerSize',30)

line([median(round(MCI_sig_freq(:),2))+std(round(MCI_sig_freq(:),2)) median(round(MCI_sig_freq(:),2))],[-1,-1],'color','k','linewidth',2)
line([median(round(MCI_sig_freq(:),2))-std(round(MCI_sig_freq(:),2)) median(round(MCI_sig_freq(:),2))],[-1,-1],'color','k','linewidth',2)
%%%%%
ylim([-10 10])
set(gca,'XLim',x1);
camroll(90)
box off


subplot(131)
h = raincloud_density(round(cntrl_sig_freq(:),2),color_code{1});
plot(median(round(cntrl_sig_freq(:),2)),-1,'.k','MarkerSize',30)

line([median(round(cntrl_sig_freq(:),2))+std(round(cntrl_sig_freq(:),2)) median(round(cntrl_sig_freq(:),2))],[-1,-1],'color','k','linewidth',2)
line([median(round(cntrl_sig_freq(:),2))-std(round(cntrl_sig_freq(:),2)) median(round(cntrl_sig_freq(:),2))],[-1,-1],'color','k','linewidth',2)
%%%%%
ylim([-10 10])
set(gca,'XLim',x1);
camroll(90)
box off
subplot(133)
xlabel('Scaling exponent (\alpha)','interpreter','Tex','FontSize',11)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(132)
yaxis_val = raincloud_yaxis(diff_scd_sig(:));
raincloud_plot(round(diff_scd_sig(:),2),yaxis_val, color_code1{2});
xlim([-0.21 0.31])
camroll(90)
x1 = get(gca,'XLim');
yl = get(gca,'YLim');
set(gca,'YLim',[-yl(2)-1.5 yl(2)-0.5]);
subplot(133)
yaxis_val = raincloud_yaxis(diff_mci_sig(:));
raincloud_plot(round(diff_mci_sig(:),2),yaxis_val, color_code1{3});
camroll(90)
set(gca,'XLim',x1);
yl = get(gca,'YLim');
set(gca,'YLim',[-yl(2)-5 yl(2)+2.5]);
subplot(131)
yaxis_val = raincloud_yaxis(diff_mciscd_sig(:));
raincloud_plot(round(diff_mciscd_sig(:),2),yaxis_val, color_code1{1});
camroll(90)
set(gca,'XLim',x1);
yl = get(gca,'YLim');
set(gca,'YLim',[-yl(2)-4.5 yl(2)+3]);
xlabel('Scaling exponent (\alpha)','interpreter','Tex','FontSize',11)


figure
subplot(132)
raincloud_density(diff_scd_sig,color_code1{2});
plot(median(diff_scd_sig),-1,'.k','MarkerSize',30)

line([median(diff_scd_sig)+std(diff_scd_sig) median(diff_scd_sig)],[-1,-1],'color','k','linewidth',2)
line([median(diff_scd_sig)-std(diff_scd_sig) median(diff_scd_sig)],[-1,-1],'color','k','linewidth',2)
%%%%%%%%%            
ylim([-10 10])
xlim([-0.21 0.31])
camroll(90)
x1 = get(gca,'XLim');
box off

subplot(133)
h = raincloud_density(diff_mci_sig,color_code1{3});
plot(median(diff_mci_sig),-1,'.k','MarkerSize',30)

line([median(diff_mci_sig)+std(diff_mci_sig) median(diff_mci_sig)],[-1,-1],'color','k','linewidth',2)
line([median(diff_mci_sig)-std(diff_mci_sig) median(diff_mci_sig)],[-1,-1],'color','k','linewidth',2)
%%%%%
ylim([-10 10])
set(gca,'XLim',x1);
camroll(90)
box off

subplot(131)
h = raincloud_density(diff_mciscd_sig,color_code1{1});
plot(median(diff_mciscd_sig),-1,'.k','MarkerSize',30)

line([median(diff_mciscd_sig)+std(diff_mciscd_sig) median(diff_mciscd_sig)],[-1,-1],'color','k','linewidth',2)
line([median(diff_mciscd_sig)-std(diff_mciscd_sig) median(diff_mciscd_sig)],[-1,-1],'color','k','linewidth',2)
%%%%%
ylim([-10 10])
set(gca,'XLim',x1);
camroll(90)
box off
subplot(133)
xlabel('DFA differences')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

















