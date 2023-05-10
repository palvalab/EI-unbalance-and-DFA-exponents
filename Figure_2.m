%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%                   Full Cohort Results                                 %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
rsdata = ['L:\nttk-data2\palva\Madrid\Full Cohort\Materials for Javed et al/Data'];
addpath ( sprintf ( '%s/Figure 2/',  rsdata) );



load([rsdata '\DFA_fEI_estimates.mat'])
for k=1:3
    for i=1:size(temp_exp{k},3)
        avg_exp_subj{k}(i,:) = nanmean(temp_exp{k}(:,:,i));
    end
    for i=1:size(temp_fEI{k},3)
        avg_fEI_subj{k}(i,:) = nanmean(temp_fEI{k}(:,:,i));
    end
end

temp_exp1 = temp_exp;    
temp_exp  =  temp_fEI;
avg_exp_subj1 = avg_exp_subj;
avg_exp_subj = avg_fEI_subj;
    
figure%('units','normalized','outerposition',[0 0 1 1])
subplot(211)
color_code = {[0.3 0.44 1],[0 0.5 0],[1 0.44 0]}; %{'b';'m';'g';'r'};%[0.47 0.47 1;1 0.53 0.1;0.56 0.82 0.56;1 0.48 0.48];
color_code1 = {[1 0 1],[0 0.5 0.5],[1 0.44 0.5]};
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
p = 95;


semilogx(MorletBank,mean(avg_exp_subj{1}),'color',color_code{1},'linewidth',2)
hold on
semilogx(MorletBank,mean(avg_exp_subj{2}),'color',color_code{2},'linewidth',2)
semilogx(MorletBank,mean(avg_exp_subj{3}),'color',color_code{3},'linewidth',2)
% plot(mean(avg_exp_subj{4}),'r','linewidth',2)


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
    adj_pval = ([1:size(rank,2)]./size(rank,2)).*0.1
    k = find(adj_pval <=0.05,1,'last')
    rank(k+1:end) = 1;
    rank(1:k) = adj_pval(1:k);
    
    ylim([0.75 1.05])
    yticks([0.75 0.8 0.85 0.9 0.95 1 1.05])
    ytick = yticks;
    xloc = find(sub_p_true <= 0.05);
    
    
    if ~isempty(xloc)
        plot(MorletBank(xloc),ytick(end)-0.02,'ks', 'MarkerSize',7);
        adj_p_align(rnk_num) = rank;
        xloc = find(adj_p_align <= 0.05);
        plot(MorletBank(xloc),ytick(end)-0.02,'rd', 'MarkerFaceColor', 'r','MarkerSize',5);
        
        st = find(adj_p_align <=0.05,1,'first')    
        lst = find(adj_p_align <=0.05,1,'last')
        indxs = xloc;
    else
        st = 1;    
        lst = size(MorletBank,2);
    end
   
    legend('NC','SCD','MCI')%,'Bonferroni','permutation')
    ylabel('Scaling exponent (\alpha)','interpreter','Tex','FontSize',14)
    legend boxoff
    xlim([2.2 90])
    xticks([3 5 10 20 30 50 90]);
    xticklabels([3 5 10 20 30 50 90])
    box off
% % % MorletBank1=roundn(MorletBank,-1);
% % % xticklabels([MorletBank1(tick1)])
ylabel('fE/I')

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
    
    adj_pval = ([1:size(rank12,2)]./size(rank12,2)).*0.1
    k = find(adj_pval <=0.05,1,'last')
    rank12(k+1:end) = 1;
    
    for i=1:size(MorletBank,2)
        if rnk_num12(i) >= st && rnk_num12(i) <= lst
            adj_p_align12(rnk_num12(i)) = rank12(i);
        else
            adj_p_align12(rnk_num12(i)) = 1;
        end
    end
    xloc = find(adj_p_align12 <= 0.05);
    tot = size(xloc,2);
    i=1;
    while i <= tot
        if isempty(find(xloc(i)==indxs))
            xloc(i)=[];
            i=i-1;
            tot = tot-1;
        end
        i=i+1;
    end       
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
    xloc = find(adj_p_align13 <= 0.05);
    tot = size(xloc,2);
    i=1;
    while i <= tot
        if isempty(find(xloc(i)==indxs))
            xloc(i)=[];
            i=i-1;
            tot = tot-1;
        end
        i=i+1;
    end 
    yloc = ypos(xloc);
    semilogx(MorletBank(xloc),yloc,'rd', 'MarkerFaceColor', 'r','MarkerSize',5,'HandleVisibility','off');
    diff_sub_exp ={};
    for i=1:size(avg_exp_subj{3})
        diff_sub_exp{1}(i,:) = avg_exp_subj{3}(i,:)-nanmean(avg_exp_subj{2});
    end

%     subplot(313)
    hold on
%     semilogx(zeros(1,90),'--','color',[0.7 0.7 0.7],'linewidth',2,'HandleVisibility','off')
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
    for i=1:size(MorletBank,2)
        if rnk_num23(i) >= st && rnk_num23(i) <= lst
            adj_p_align23(rnk_num23(i)) = rank23(i);
        else
            adj_p_align23(rnk_num23(i)) = 1;
        end
    end

    xloc = find(adj_p_align23 <= 0.05);
    tot = size(xloc,2);
    i=1;
    while i <= tot
        if isempty(find(xloc(i)==indxs))
            xloc(i)=[];
            i=i-1;
            tot = tot-1;
        end
        i=i+1;
    end 
    yloc = ypos(xloc);
    semilogx(MorletBank(xloc),yloc,'rd', 'MarkerFaceColor', 'r','MarkerSize',5,'HandleVisibility','off');
    xlim([2.2 90])
    xticks([3 5 10 20 30 50 90]);
    xticklabels([3 5 10 20 30 50 90])
    
    xlabel('Frequency (Hz)')
    ylabel('fE/I differences ','interpreter','Tex','FontSize',12)
    legend('SCD-NC','MCI-NC','MCI-SCD')
    legend boxoff
    box off
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
% % %     
% % %     
% % % display('       Cntrl vs SCD')
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cntrl vs SCD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % figure
% % % grouplabel12 = [repmat({'cntrl'},size(temp_exp{1},3),1); repmat({'scd'},size(temp_exp{2},3),1)]; %...
% % % %     ; repmat({'mci'},size(temp_exp{3},3),1)];
% % % 
% % % 
% % % p_true12 = [];tbl_true12={}; h_true12=[];
% % % true_distribution12 =[];
% % % 
% % % for i=1:size(temp_exp{1},1)
% % %     x=[];
% % %     exp_per_parcel = [];
% % %     for k=1:2
% % %         temp_hold = [];
% % %         temp_hold(:,:) = temp_exp{k}(i,:,:);
% % %         exp_per_parcel = [exp_per_parcel;  temp_hold'];
% % %     end
% % %     for m=1:size(exp_per_parcel,2)
% % %         [p_true12(i,m),h_true12(i,m),tbl_true12{i,m}] = ranksum(exp_per_parcel(1:size(temp_exp{1},3),m),exp_per_parcel(size(temp_exp{1},3)+1:end,m));
% % %         true_distribution12(i,m) = tbl_true12{i,m}.zval;
% % %     end    
% % % end
% % % p_true12_1d= p_true12(:)';
% % % % % % save('-v7.3','L:\nttk-data2\palva\Madrid\Full Cohort\Figures Codes\Between group fEI stats\Full Cohort\wilcoxn_p_true12_1d','true_distribution12','MorletBank','p_true12_1d');
% % % uncorr_sig_05 = zeros(size(temp_exp{k},1),size(exp_per_parcel,2));
% % % uncorr_sig_05(find(p_true12 <= 0.05))=1;
% % % imagesc(true_distribution12)
% % % hold on
% % % contour(logical(uncorr_sig_05),1,'linecolor','k');
% % % hcb = colorbar;
% % % Bonferroni_sig_05 = zeros(size(temp_exp{k},1),size(exp_per_parcel,2));
% % % Bonferroni_sig_05(find(p_true12 <= (0.05/(400*32))))=1;
% % % contour(logical(Bonferroni_sig_05),1,'linecolor','b');
% % % 
% % % display('       Cntrl vs MCI')
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cntrl vs MCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % figure
% % % grouplabel13 = [repmat({'cntrl'},size(temp_exp{1},3),1); repmat({'mci'},size(temp_exp{3},3),1)];
% % % 
% % % 
% % % p_true13 = [];tbl_true13={};
% % % true_distribution13 =[];
% % % 
% % % for i=1:size(temp_exp{1},1)
% % %     x=[];
% % %     exp_per_parcel = [];
% % %     for k=1:2:3
% % %         temp_hold = [];
% % %         temp_hold(:,:) = temp_exp{k}(i,:,:);
% % %         exp_per_parcel = [exp_per_parcel;  temp_hold'];
% % %     end
% % %     for m=1:size(exp_per_parcel,2)
% % %         [p_true13(i,m),h_true13(i,m),tbl_true13{i,m}] = ranksum(exp_per_parcel(1:size(temp_exp{1},3),m),exp_per_parcel(size(temp_exp{1},3)+1:end,m));
% % %         true_distribution13(i,m) = tbl_true13{i,m}.zval;
% % %     end    
% % % end
% % % p_true13_1d= p_true13(:)';
% % % % % % save('-v7.3','L:\nttk-data2\palva\Madrid\Full Cohort\Figures Codes\Between group fEI stats\Full Cohort\wilcoxn_p_true13_1d','true_distribution13','MorletBank','p_true13_1d');
% % % 
% % % uncorr_sig_05 = zeros(size(temp_exp{1},1),size(exp_per_parcel,2));
% % % uncorr_sig_05(find(p_true13 <= 0.05))=1;
% % % imagesc(true_distribution13)
% % % hold on
% % % contour(logical(uncorr_sig_05),1,'linecolor','k');
% % % hcb = colorbar;
% % % Bonferroni_sig_05 = zeros(size(temp_exp{1},1),size(exp_per_parcel,2));
% % % Bonferroni_sig_05(find(p_true13 <= (0.05/(400*32))))=1;
% % % contour(logical(Bonferroni_sig_05),1,'linecolor','b');
% % % display('       SCD vs MCI')
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCD vs MCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % figure
% % % grouplabel23 = [repmat({'scd'},size(temp_exp{2},3),1); repmat({'mci'},size(temp_exp{3},3),1)];
% % % 
% % % 
% % % p_true23 = [];tbl_true23={};
% % % true_distribution23 =[];
% % % 
% % % for i=1:size(temp_exp{1},1)
% % %     x=[];
% % %     exp_per_parcel = [];
% % %     for k=2:3
% % %         temp_hold = [];
% % %         temp_hold(:,:) = temp_exp{k}(i,:,:);
% % %         exp_per_parcel = [exp_per_parcel;  temp_hold'];
% % %     end
% % %     for m=1:size(exp_per_parcel,2)
% % %         [p_true23(i,m),h_true23(i,m),tbl_true23{i,m}] = ranksum(exp_per_parcel(1:size(temp_exp{2},3),m),exp_per_parcel(size(temp_exp{2},3)+1:end,m));
% % %         true_distribution23(i,m) = tbl_true23{i,m}.zval;
% % %     end    
% % % end
% % % p_true23_1d= p_true23(:)';
% % % % % % save('-v7.3','L:\nttk-data2\palva\Madrid\Full Cohort\Figures Codes\Between group fEI stats\Full Cohort\wilcoxn_p_true23_1d','true_distribution23','MorletBank','p_true23_1d');
% % % uncorr_sig_05 = zeros(size(temp_exp{1},1),size(exp_per_parcel,2));
% % % uncorr_sig_05(find(p_true23 <= 0.05))=1;
% % % imagesc(true_distribution23)
% % % hold on
% % % contour(logical(uncorr_sig_05),1,'linecolor','k');
% % % hcb = colorbar;
% % % Bonferroni_sig_05 = zeros(size(temp_exp{1},1),size(exp_per_parcel,2));
% % % Bonferroni_sig_05(find(p_true23 <= (0.05/(400*32))))=1;
% % % contour(logical(Bonferroni_sig_05),1,'linecolor','b');
% % % 
% % % 
% % % 
% close all
clearvars -except rsdata
clc 
load([rsdata '\Figure 2\wilcoxn_p_true12_1d.mat'])
color_code1 = {[1 0 1],[0 0.5 0.5],[1 0.44 0.5]};
Q = 0.2;
p_true12    = reshape(p_true12_1d, size(true_distribution12));
uncorr_sig_05 = zeros(size(true_distribution12));
uncorr_sig_05(find(p_true12 <= 0.05))=1;
    [rank12 rnk_num12]= sort(p_true12_1d);
    k = find(rank12 <=0.05,1,'last')
    k_reduction = k-round(k*Q)
    rank12(k_reduction+1:end) = 1;
    adj_p_align12(rnk_num12) = rank12;
    
    fdr_adj_pval    = reshape(adj_p_align12, size(true_distribution12));
    fdr_sig_05 = zeros(size(true_distribution12));
    fdr_sig_05(find(fdr_adj_pval <= 0.05))=1;
    
    SCD_fdr_sig_05 = fdr_sig_05;
    numel(find(adj_p_align12 <= 0.05))
    numel(find(p_true12_1d <= 0.05))

uncorr_sig_prcl12 = [];
fdr_sig_prcl12 = [];

for m=1:size(true_distribution12,2)
    uncorr_sig_prcl12(m) = numel(find(uncorr_sig_05(:,m)==1))/400;
    fdr_sig_prcl12(m) = numel(find(fdr_sig_05(:,m)==1))/400;
end

figure(12)
I2=semilogx(MorletBank,fdr_sig_prcl12,'-','color',color_code1{2},'linewidth',2);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cntrl vs MCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except parcel_perm_distribution FDR_SCD pi0_SCD I1 I2 Q12 SCD_fdr_sig_05 Q color_code1 rsdata
load([rsdata '\Figure 2\wilcoxn_p_true13_1d.mat'])
p_true13    = reshape(p_true13_1d, size(true_distribution13)); 
uncorr_sig_05 = zeros(size(true_distribution13));
uncorr_sig_05(find(p_true13 <= 0.05))=1;
[rank13 rnk_num13]= sort(p_true13_1d);
    k = find(rank13 <=0.05,1,'last')
    k_reduction = k-round(k*Q)
    rank13(k_reduction+1:end) = 1;
    adj_p_align13(rnk_num13) = rank13;
    
    fdr_adj_pval    = reshape(adj_p_align13, size(true_distribution13));
    fdr_sig_05 = zeros(size(true_distribution13));
    fdr_sig_05(find(fdr_adj_pval <= 0.05))=1;
    
    MCI_fdr_sig_05 = fdr_sig_05;
    numel(find(adj_p_align13 <= 0.05))
    numel(find(p_true13_1d <= 0.05))

uncorr_sig_prcl13 = [];
fdr_sig_prcl13 = [];
for m=1:size(true_distribution13,2)
    uncorr_sig_prcl13(m) = numel(find(uncorr_sig_05(:,m)==1))/400;
    fdr_sig_prcl13(m) = numel(find(fdr_sig_05(:,m)==1))/400;
end

h1 = figure(12)
hold on
I4=semilogx(MorletBank,fdr_sig_prcl13,'-','color',color_code1{3},'linewidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCD vs MCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except parcel_perm_distribution FDR_SCD pi0_SCD I1 I2 FDR_MCI pi0_MCI I3 I4 Q12 Q13 SCD_fdr_sig_05 MCI_fdr_sig_05 Q color_code1 rsdata
load([rsdata '\Figure 2\wilcoxn_p_true23_1d.mat'])
p_true23    = reshape(p_true23_1d, size(true_distribution23));
uncorr_sig_05 = zeros(size(true_distribution23));
uncorr_sig_05(find(p_true23 <= 0.05))=1;
    [rank23 rnk_num23]= sort(p_true23_1d);
%     adj_pval = ([1:size(rank,2)]./size(rank,2)).*0.1;
    k = find(rank23 <=0.05,1,'last')
    k_reduction = k-round(k*Q)
    rank23(k_reduction+1:end) = 1;
    adj_p_align23(rnk_num23) = rank23;
    
    fdr_adj_pval    = reshape(adj_p_align23, size(true_distribution23));
    fdr_sig_05 = zeros(size(true_distribution23));
    fdr_sig_05(find(fdr_adj_pval <= 0.05))=1;
    
    SM_fdr_sig_05 = fdr_sig_05;
    numel(find(adj_p_align23 <= 0.05))
    numel(find(p_true23_1d <= 0.05))

uncorr_sig_prcl23 = [];
fdr_sig_prcl23 = [];
for m=1:size(true_distribution23,2)
    uncorr_sig_prcl23(m) = numel(find(uncorr_sig_05(:,m)==1))/400;
    fdr_sig_prcl23(m) = numel(find(fdr_sig_05(:,m)==1))/400;
end

h1 = figure(12)
hold on
I6=semilogx(MorletBank,fdr_sig_prcl23,'-','color',color_code1{1},'linewidth',2);
hl = legend()
str=['Q=', sprintf('%.2f', Q*100),' %'];
str1=['Q=', sprintf('%.2f', Q*100),' %'];
str2=['Q=', sprintf('%.2f', Q*100),' %'];
legend_str = {'SCD-NC','MCI-NC','MCI-SCD'};
set(hl, 'string', legend_str);
legend boxoff
box off
xlim([2.2 90])
xticks([2.2 5 10 20 30 50 90]);
xticklabels([2 5 10 20 30 50 90])
h = gca; % Get axis to modify
h.XAxis.MinorTick = 'on';
ylabel('% of significant (p < 0.05) parcels')
xlabel('Frequency (Hz)')

% % % save('-v7.3','L:\nttk-data2\palva\Madrid\Full Cohort\Figures Codes\Between group fEI stats\Full Cohort\wilcoxn_FDR_Sig_05','MCI_fdr_sig_05','SCD_fdr_sig_05','SM_fdr_sig_05');



clearvars -except rsdata
load([rsdata '\DFA_fEI_estimates.mat'])

load([rsdata '\Figure 2\py_aligned_network_for_inflated_surface_plot.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([rsdata '\Figure 2\wilcoxn_FDR_Sig_05'])
bnd = [{12:15},{16:22}];

display('SCD vs NC')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cntrl vs SCD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig_per_freq_scd=[];
for i=1:size(bnd{1},2)
    sig_per_freq_scd(1,i) = numel(find(SCD_fdr_sig_05(:,bnd{1}(i))==1))/400;
end
[val,idx]=max(sig_per_freq_scd)

for i=1:size(bnd{1},2)
    SCD_fdr_sig_05(:,bnd{1}(i)) = SCD_fdr_sig_05(:,bnd{1}(i)).*SCD_fdr_sig_05(:,bnd{1}(idx));
end
sig_per_freq_scd=[];
for i=1:size(bnd{2},2)
    sig_per_freq_scd(1,i) = numel(find(SCD_fdr_sig_05(:,bnd{2}(i))==1))/400;
end
[val,idx]=max(sig_per_freq_scd)

for i=1:size(bnd{2},2)
    SCD_fdr_sig_05(:,bnd{2}(i)) = SCD_fdr_sig_05(:,bnd{2}(i)).*SCD_fdr_sig_05(:,bnd{2}(idx));
end

avg_nc(:,:) = nanmean(temp_exp{1},3);
for i=1:size(temp_exp{2},3)
    diff12(:,:) = temp_exp{2}(:,:,i)-avg_nc;
    SCD_fdr_sig_05(find(SCD_fdr_sig_05==0))=NaN;
    prcl_diff12(:,:,i) = diff12.*SCD_fdr_sig_05;
end

display('MCI vs NC')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cntrl vs MCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig_per_freq_MCI=[];
for i=1:size(bnd{1},2)
    sig_per_freq_MCI(1,i) = numel(find(MCI_fdr_sig_05(:,bnd{1}(i))==1))/4;
end
[val,idx]=max(sig_per_freq_MCI)

for i=1:size(bnd{1},2)
    MCI_fdr_sig_05(:,bnd{1}(i)) = MCI_fdr_sig_05(:,bnd{1}(i)).*MCI_fdr_sig_05(:,bnd{1}(idx));
end
sig_per_freq_MCI=[];
for i=1:size(bnd{2},2)
    sig_per_freq_MCI(1,i) = numel(find(MCI_fdr_sig_05(:,bnd{2}(i))==1))/400;
end
[val,idx]=max(sig_per_freq_MCI)

for i=1:size(bnd{2},2)
    MCI_fdr_sig_05(:,bnd{2}(i)) = MCI_fdr_sig_05(:,bnd{2}(i)).*MCI_fdr_sig_05(:,bnd{2}(idx));
end

for i=1:size(temp_exp{3},3)
    diff13(:,:) = temp_exp{3}(:,:,i)-avg_nc;
    MCI_fdr_sig_05(find(MCI_fdr_sig_05==0))=NaN;
    prcl_diff13(:,:,i) = diff13.*MCI_fdr_sig_05;    
end

display('       SCD vs MCI')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCD vs MCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig_per_freq_SM=[];
for i=1:size(bnd{1},2)
    sig_per_freq_SM(1,i) = numel(find(SM_fdr_sig_05(:,bnd{1}(i))==1))/400;
end
[val,idx]=max(sig_per_freq_SM)

for i=1:size(bnd{1},2)
    SM_fdr_sig_05(:,bnd{1}(i)) = SM_fdr_sig_05(:,bnd{1}(i)).*SM_fdr_sig_05(:,bnd{1}(idx));
end
sig_per_freq_SM=[];
for i=1:size(bnd{2},2)
    sig_per_freq_SM(1,i) = numel(find(SM_fdr_sig_05(:,bnd{2}(i))==1))/400;
end
[val,idx]=max(sig_per_freq_SM)

for i=1:size(bnd{2},2)
    SM_fdr_sig_05(:,bnd{2}(i)) = SM_fdr_sig_05(:,bnd{2}(i)).*SM_fdr_sig_05(:,bnd{2}(idx));
end

avg_scd(:,:) = nanmean(temp_exp{2},3);
for i=1:size(temp_exp{3},3)
    diff23(:,:) = temp_exp{3}(:,:,i)-avg_scd;
    SM_fdr_sig_05(find(SM_fdr_sig_05==0))=NaN;
    prcl_diff23(:,:,i) = diff23.*SM_fdr_sig_05;    
end


alpha_diff12 = nanmean(nanmean(prcl_diff12(:,bnd{1},:),3),2);
alpha_diff12(find(isnan(alpha_diff12))) = 0;
alpha_diff12 = alpha_diff12(py_aligned_network);
fname = ['L:\nttk-data2\palva\Madrid\Full Cohort\Figures Codes\Paper-figures codes\Inflated_Surf_data\fEI\fEIresults_alpha_diff12_max_freq.csv' ];
%         csvwrite(fname,alpha_diff12)

alpha_diff13 = nanmean(nanmean(prcl_diff13(:,bnd{1},:),3),2);
alpha_diff13(find(isnan(alpha_diff13))) = 0;
alpha_diff13 = alpha_diff13(py_aligned_network);
fname = ['L:\nttk-data2\palva\Madrid\Full Cohort\Figures Codes\Paper-figures codes\Inflated_Surf_data\fEI\fEIresults_alpha_diff13_max_freq.csv' ];
%         csvwrite(fname,alpha_diff13)

alpha_diff23 = nanmean(nanmean(prcl_diff23(:,bnd{1},:),3),2);
alpha_diff23(find(isnan(alpha_diff23))) = 0;
alpha_diff23 = alpha_diff23(py_aligned_network);
fname = ['L:\nttk-data2\palva\Madrid\Full Cohort\Figures Codes\Paper-figures codes\Inflated_Surf_data\fEI\fEIresults_alpha_diff23_max_freq.csv' ];
%         csvwrite(fname,alpha_diff23)


beta_diff12 = nanmean(nanmean(prcl_diff12(:,bnd{2},:),3),2);
beta_diff12(find(isnan(beta_diff12))) = 0;
beta_diff12 = beta_diff12(py_aligned_network);
fname = ['L:\nttk-data2\palva\Madrid\Full Cohort\Figures Codes\Paper-figures codes\Inflated_Surf_data\fEI\fEIresults_beta_diff12_max_freq.csv' ];
%         csvwrite(fname,beta_diff12)

beta_diff13 = nanmean(nanmean(prcl_diff13(:,bnd{2},:),3),2);
beta_diff13(find(isnan(beta_diff13))) = 0;
beta_diff13 = beta_diff13(py_aligned_network);
fname = ['L:\nttk-data2\palva\Madrid\Full Cohort\Figures Codes\Paper-figures codes\Inflated_Surf_data\fEI\fEIresults_beta_diff13_max_freq.csv' ];
%         csvwrite(fname,beta_diff13)

beta_diff23 = nanmean(nanmean(prcl_diff23(:,bnd{2},:),3),2);
beta_diff23(find(isnan(beta_diff23))) = 0;
beta_diff23 = beta_diff23(py_aligned_network);
fname = ['L:\nttk-data2\palva\Madrid\Full Cohort\Figures Codes\Paper-figures codes\Inflated_Surf_data\fEI\fEIresults_beta_diff23_max_freq.csv' ];
%         csvwrite(fname,beta_diff23)


%%
%Parcel Based Correlation between DFA and fE/I

band ={[1:5],[6:11],[12:15],[16:22]};
for bnd = 3:4
    band_prcls_exp = [];
    band_prcls_fEI = [];
    for k=1:3
        temp=[];
        temp(:,:) = nanmean(temp_exp{k}(:,band{bnd},:),2);
        band_prcls_exp = [band_prcls_exp temp];
        
        temp=[];
        temp(:,:) = nanmean(temp_fEI{k}(:,band{bnd},:),2);
        band_prcls_fEI = [band_prcls_fEI temp];
        
        for i=1:size(band_prcls_fEI,1)
            [Corr_fEI_DFA_pearson_prcls{1}{bnd-2}(i) pval_prcls{bnd-2}(i,1)] = corr(band_prcls_exp(i,:)',band_prcls_fEI(i,:)','Type','pearson');
        end
    end
end

% % for i=1:size(Corr_fEI_DFA_pearson_prcls{1},2)
% %     Corr_fEI_DFA_pearson_prcls{1}{i}(isnan(Corr_fEI_DFA_pearson_prcls{1}{i})) = [];
% %     pval_prcls{i}(isnan(pval_prcls{i})) = [];
% % end

sig=0.05;
Q = 0.05;
[fdr_sig_05_DFA_fEI]=fdr_correct_sig05([pval_prcls{1},pval_prcls{2}],sig,Q); numel(find(fdr_sig_05_DFA_fEI))
fdr_sig_05_DFA_fEI(fdr_sig_05_DFA_fEI==0) = NaN;

% % % cases = {'Alpha band', 'Beta band'};
% % % Origem = cases;
% % %     figure, pirateplot(Corr_fEI_DFA_pearson_prcls{1},Origem,0.95);
% % %     ylabel('correlation coefficient');
% % % 
% % %     
    data(:,1) = Corr_fEI_DFA_pearson_prcls{1}{1};
    data(:,2) = Corr_fEI_DFA_pearson_prcls{1}{2};
% % %     
% % %     figure
% % %     vcolor =[0 0 0.562500000000000;0.500000000000000 0 0];
% % %     for i=1:2        
% % %         x=data(:,i);
% % %         vec_X =[1,2];
% % %         xScat{i} = vec_X(i)*ones(1,length(x)) + 0.1*(rand(1,length(x))-0.5);
% % %         hold on, scatter(xScat{i},x,floor(40/2),'filled',...
% % %             'MarkerFaceColor',vcolor(i,:), 'MarkerEdgeColor',[0,0,0],...
% % %             'MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15);        
% % %     end    
% % %     boxplot(data, vec_X);
% % %     xticklabels(Origem)
    
data = data.*fdr_sig_05_DFA_fEI;
for i=1:size(data,2)
            temp=[]; temp = data(:,i);
            temp(isnan(temp)) = [];
            data1{i} = temp;            
end
            

    figure
    vcolor =[0.5 0.3 0.562500000000000;0.500000000000000 0.4 0.1];
    for i=1:size(data,2)        
        x=data1{i};
        subplot(1,2,i)

        yaxis_val = raincloud_yaxis(x(:));
        raincloud_plot_color(round(x(:),2),yaxis_val, vcolor(i,:));        
        xlim([-0.5 0])
        camroll(90)
        x1 = get(gca,'XLim');
        yl = get(gca,'YLim');
        set(gca,'YLim',[-yl(2)-5 yl(2)+3]);
    end
    
 %%

f3 = figure
color_code = {[0.3 0.44 1],[0 0.5 0],[1 0.44 0]};
band ={[1:5],[6:11],[12:15],[16:22]};
for bnd = 3:4
    band_prcls_exp = [];
    band_prcls_fEI = [];
    for k=1:3
        temp=[];
        temp(:,:) = nanmean(temp_exp{k}(:,band{bnd},:),2);
        band_prcls_exp = [band_prcls_exp temp];
        
        temp=[];
        temp(:,:) = nanmean(temp_fEI{k}(:,band{bnd},:),2);
        band_prcls_fEI = [band_prcls_fEI temp];
    end
    band_prcls_fEI(find(isnan(band_prcls_fEI))) = 1;

    for i=1:size(band_prcls_exp,1)    
        y1=[];x1=[];
        x1 = band_prcls_exp(i,:);
        y1 = band_prcls_fEI(i,:);

        [rho_prcl{1}{bnd-2}(i),p_rho_prcl{bnd-2}(i)] = quad_fit_all(y1(:),x1(:),color_code{k},2,'fE/I','DFA-exp',k-1);
        quad_fit_all(y1(:),x1(:),'r',2,'fE/I','DFA-exp',k);
    end    
end
close(f3)

sig=0.05;
Q = 0.1;
[fdr_sig_05_quad_DFA_fEI]=fdr_correct_sig05([p_rho_prcl{1};p_rho_prcl{2}]',sig,Q); numel(find(fdr_sig_05_quad_DFA_fEI))
fdr_sig_05_quad_DFA_fEI(fdr_sig_05_quad_DFA_fEI==0) = NaN;

data = [];
data(:,1) = rho_prcl{1}{1};
data(:,2) = rho_prcl{1}{2};

data = data.*fdr_sig_05_quad_DFA_fEI;
for i=1:size(data,2)
            temp=[]; temp = data(:,i);
            temp(isnan(temp)) = [];
            data1{i} = temp;            
end

    figure
    vcolor =[0.7 0.5 0.762500000000000;0.900000000000000 0.4 0.5];
    for i=1:size(data1,2)        
        x=data1{i};
        subplot(1,2,i)

        yaxis_val = raincloud_yaxis(x(:));
        raincloud_plot_color(round(x(:),2),yaxis_val, vcolor(i,:));        
        xlim([0 0.3])
        camroll(90)
        x1 = get(gca,'XLim');
        yl = get(gca,'YLim');
        set(gca,'YLim',[-yl(2)-5 yl(2)+3]);
    end







