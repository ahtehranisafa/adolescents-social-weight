clc
clearvars
close all
Num_Subject =38; % number of subjects
%
ML_Rho  = NaN (Num_Subject,3);
ML_Beta = NaN(Num_Subject,3);
ML_NLL  = NaN(Num_Subject,3);
Prediction_Performance_Time = NaN(Num_Subject,4);
Prediction_Performance  = NaN(Num_Subject,1);
Tom2USD =4.2;
%Tom2USD =1;
peer_choice_averse =[]; peer_choice_seeker = [];
for Subject  = 1:Num_Subject
    load(['SS',num2str(100+Subject),'.mat']);
    u_self1(:,2) =  u_self1(:,2)/Tom2USD;
    u_self1(:,5) =  u_self1(:,5)/Tom2USD;
    u_self3(:,2) =  u_self3(:,2)/Tom2USD;
    u_self3(:,5) =  u_self3(:,5)/Tom2USD;
    u_prediction(:,2)=u_prediction(:,2)/Tom2USD;
    u_prediction(:,5)=u_prediction(:,5)/Tom2USD;
    %% Prediction Analysis
    Prediction_Performance(Subject,1) = mean(y_prediction(6:35,1)==u_prediction(6:35,1));
    Prediction_Performance_Time(Subject,:) = [mean(y_prediction(1:5,1)==u_prediction(1:5,1))...
        mean(y_prediction(6:15,1)==u_prediction(6:15,1))...
        mean(y_prediction(16:25,1)==u_prediction(16:25,1))...
        mean(y_prediction(26:35,1)==u_prediction(26:35,1))];
    %% Model Fitting
    % self1
    [ML_Rho(Subject,1),ML_Beta(Subject,1),ML_model_free(Subject,1),ML_NLL(Subject,1)]...
        = ML_fitting(u_self1);
    % observee
    [ML_Rho(Subject,2),ML_Beta(Subject,2),ML_model_free(Subject,2), ML_NLL(Subject,2)]...
        = ML_fitting(u_prediction);
    %ML_Rho(Subject,2) = round(ML_Rho(Subject,2),1);
    %self3
    [ML_Rho(Subject,3),ML_Beta(Subject,3),ML_model_free(Subject,3),ML_NLL(Subject,3)]...
        = ML_fitting(u_self3);
    clear u_self1 u_self3 u_prediction y_prediction
    clear data_choice data_mag data_prob U_sure
end             
%% Exlusion Criteria
% subject #15 :choose sure option in check trials
% Subject#13 and #21 Prediction Performance is less than %60
Min_Accepted_Pre = 0.633;
idx_Exclusion1 = find(Prediction_Performance<=Min_Accepted_Pre);
idx_Exclusion2 = find(ML_Beta(:,1)==0);
idx_Exclusion= [15 ; idx_Exclusion1; idx_Exclusion2]';
%% ML results: After Exclusion
ML_results.Rho = ML_Rho;
ML_results.Rho(idx_Exclusion,:) =NaN;
ML_results.Beta = ML_Beta;
ML_results.model_free = ML_model_free;
ML_results.model_free(idx_Exclusion,:) =NaN;
ML_results.NLL = ML_NLL;
ML_results.NLL(idx_Exclusion,:) =NaN;
Prediction_Performance(idx_Exclusion,:) =NaN;
clear ML_Rho ML_Beta  ML_model_free ML_NLL
%% risk attitude profile
Rho_mean = nanmean(ML_results.Rho,1);
Rho_sd   = nanstd (ML_results.Rho,1);
Rho_median = nanmedian(ML_results.Rho,1);
Rho_min = min(ML_results.Rho);
Rho_max = max(ML_results.Rho);
C1 = cell(5,2); C1{1,1} = 'mean';C1{2,1} = 'sd'; C1{3,1} = 'median';
C1{4,1} = 'min'; C1{5,1} = 'max';
C1{1,2} = Rho_mean(1); C1{2,2} = Rho_sd(1);
C1{3,2} = Rho_median(1);C1{4,2} = Rho_min(1);C1{5,2} = Rho_max(1);
Rho_profile = C1 ;
%% baseline risk attitude 
%  Comparing Risk attitudes in session 1 between participants who observed
%  risk-averse peers and participants who observed risk-seekers
%  group1:participants who observed risk-averse peers 
%  group2: participants who observed risk-seeking peers
Rho_group1 = ML_results.Rho (ML_results.model_free(:,2)<0.5,:);
Rho_group2 = ML_results.Rho (ML_results.model_free(:,2)>0.5,:);
[StatTest_baseline_g1vsg2.h,StatTest_baseline_g1vsg2.p,...
StatTest_baseline_g1vsg2.ci,StatTest_baseline_g1vsg2.stat]  = ...
    ttest2 (Rho_group1(:,1),Rho_group2(:,1));
%% Contagion
Contagion_group1 = Rho_group1(:,1) - Rho_group1(:,3);
Contagion_group2 = Rho_group2(:,3) - Rho_group2(:,1);
Contagion_all = [Contagion_group1;Contagion_group2];
% contagion_all profile
Contagion_all_mean = nanmean(Contagion_all);
Contagion_all_sd   = nanstd (Contagion_all);
Contagion_all_median = nanmedian(Contagion_all);
Contagion_all_min = min(Contagion_all);
Contagion_all_max = max(Contagion_all);
C2 = cell(5,2); C2{1,1} = 'mean';C2{2,1} = 'sd'; C2{3,1} = 'median';
C2{4,1} = 'min'; C2{5,1} = 'max';
C2{1,2} = Contagion_all_mean; C2{2,2} = Contagion_all_sd;
C2{3,2} = Contagion_all_median;C2{4,2} = Contagion_all_min;C2{5,2} = Contagion_all_max;
Contagion_all_profile = C2 ;
% one-sided ttest for group1: H1: Contagion is positive
[StatTest_Contagion_group1.h,StatTest_Contagion_group1.p,...
StatTest_Contagion_group1.ci,StatTest_Contagion_group1.stat]  = ... 
ttest(Contagion_group1,0,'Tail','right');% one-sided ttest: H1: Contagion is positive
% one-sided ttest for group2: H1: Contagion is positive
[StatTest_Contagion_group2.h,StatTest_Contagion_group2.p,...
StatTest_Contagion_group2.ci,StatTest_Contagion_group2.stat]  = ... 
ttest(Contagion_group2,0,'Tail','right');
% compare contagion  between two groups
[StatTest_contagion_g1vsg2.h,StatTest_contagion_g1vsg2.p,...
StatTest_contagion_g1vsg2.ci,StatTest_contagion_g1vsg2.stat]  = ... 
ttest2 (Contagion_group1,Contagion_group2); 
% one-sided ttest for all participants  contagion  H1: Contagion is positive
[StatTest_Contagion_all.h,StatTest_Contagion_all.p,...
StatTest_Contagion_all.ci,StatTest_Contagion_all.stat]  = ... 
ttest(Contagion_all,0,'Tail','right');
%% Statistical information about choices
% participants are seprated based on whom observed
Choice_Av = ML_results.model_free(ML_results.model_free(:,2)<0.5,:);
M_Choice_Av = mean(Choice_Av,1) ;
SD_Choice_Av = std(Choice_Av,1);
L_av = min(Choice_Av) ;
U_av= max(Choice_Av);
Choice_profile_av = [L_av' U_av' M_Choice_Av' SD_Choice_Av'];
Choice_Ta = ML_results.model_free(ML_results.model_free(:,2)>0.5,:);
M_Choice_ta = mean(Choice_Ta,1) ;
SD_Choice_ta = std(Choice_Ta,1);
L_ta = min(Choice_Ta) ;
U_ta= max(Choice_Ta);
Choice_profile_ta = [L_ta' U_ta' M_Choice_ta' SD_Choice_ta']
%% Statistical test: Choices
[StatTest_Choice_Av.h,StatTest_Choice_Av.p,StatTest_Choice_Av.ci,StatTest_Choice_Av.stat]  =...
    ttest (Choice_Av(:,1),Choice_Av(:,3),'Tail','right');
[StatTest_Choice_Ta.h,StatTest_Choice_Ta.p,StatTest_Choice_Ta.ci,StatTest_Choice_Ta.stat]  = ...
    ttest (Choice_Ta(:,1),Choice_Ta(:,3),'Tail','left');
[StatTest_Choice_AvVsTa.h,StatTest_Choice_AvVsTa.p,StatTest_Choice_AvVsTa.ci,StatTest_Choice_AvVsTa.stat]  = ...
    ttest2 (Choice_Av(:,1),Choice_Ta(:,1))

%% statistical tests: prediction session
Prediction_Performance_Time(idx_Exclusion,:) =NaN;
[StatTest_Prediction.h,StatTest_Prediction.p,StatTest_Prediction.ci,StatTest_Prediction.stat]=...
    ttest(Prediction_Performance_Time,0.633)
%% social shift vs social distance
delta_OS = ML_results.Rho(:,2)-ML_results.Rho(:,1);% distance in peers risk attitude
delta_SS = ML_results.Rho(:,3)-ML_results.Rho(:,1);% change of risk attitude
ML_Rho_O_av = ML_results.Rho(ML_results.Rho(:,2)<1,:); % data for subject who obsevre aversive agent
ML_Rho_O_ta = ML_results.Rho(ML_results.Rho(:,2)>1,:); % data for subject who obsevre risk taker agent
Contagion = delta_SS.* sign(ML_results.Rho(:,2)-1);
Contagion = delta_SS.* sign(delta_OS);
Contagion_model_free = (ML_results.model_free(:,3)-ML_results.model_free(:,1)).*...
    sign(ML_results.model_free(:,2)-0.5);
Contagion_model_free(idx_Exclusion,:) =NaN;
[r,p]=corr(Contagion, Prediction_Performance,'rows','complete')
%% statistical test: Contagion
Contagion_Mean = nanmean(Contagion);
Contagion_STD  = nanstd(Contagion);
Contagion_normalized = (Contagion-Contagion_Mean)./Contagion_STD;
Contagion_StatTest.Mean = Contagion_Mean;
Contagion_StatTest.STD  = Contagion_STD;
clear Contagion_Mean Contagion_STD
[h,p,ksstat,cv] = kstest(Contagion_normalized); % normality test
Contagion_StatTest.kstest.h = h;
Contagion_StatTest.kstest.p = p;
Contagion_StatTest.kstest.ksstat = ksstat;
Contagion_StatTest.kstest.cv = cv;
clear h p ksstat cv
[h,p,ci,stats] = ttest(Contagion',0,'Tail','right');% one-sided ttest: H1: Contagion is positive
Contagion_StatTest.ttest.h = h;
Contagion_StatTest.ttest.p = p;
Contagion_StatTest.ttest.ci = ci;
Contagion_StatTest.ttest.stats = stats;
clear h p ci stats
[StatTest_Beta.h,StatTest_Beta.p,StatTest_Beta.ci,StatTest_Beta.stats] = ...
    ttest(ML_results.Beta(:,3) , ML_results.Beta(:,1) );
[R,P_value]=corr(delta_OS,delta_SS,'rows','complete');% correlation test
SSOS_Corr_test.R = R;
SSOS_Corr_test.P_value = P_value;
clear R P_value
[R,P_value]= corr(Prediction_Performance,Contagion,'rows','complete');
COCP_Corr_test.R = R;
COCP_Corr_test.P_value = P_value;
clear R P_value
%% Figure S-3. Baseline risk attitudes
X1 = ML_results.Rho(:,1);
X1 = X1(~isnan(X1));
[Xvector1,idx2] = sort(X1); % sort participants respect to their Rho in Session 1
% Create figure
fig1 = figure('Name','Participants'' Risk-preference');
% Create axes
axes1 = axes('Parent',fig1,...
    'Position',[0.2 0.12 0.60 0.80]);
hold(axes1,'on');
% Create bar
b= bar(Xvector1,'Parent',axes1,'FaceColor',[.65 0.65 0.65],'EdgeColor',[1 1 1]);
b.BaseValue =1;
 M_O_ta =mean (ML_Rho_O_ta(:,2)) 
 S_o_ta = std (ML_Rho_O_ta(:,2)) 
 M_O_av = mean(ML_Rho_O_av(:,2)) 
 S_o_av = std (ML_Rho_O_av(:,2))
l51=line([0 40],[M_O_ta M_O_ta],'Color','b');
l51.LineWidth = 2;
l52 = line([0 40],[M_O_av...
    M_O_av],'Color','b');
l52.LineWidth = 2;
% Create ylabel
ylabel('$ Risk\  Attitude\ in\ session\ one\  (\rho_{S1})\  $','interpreter','latex','Fontsize',14);
% Create xlabel
xlabel('$Participants$','interpreter','latex','Fontsize',16);
% Set the remaining axes properties
xlim([0 size(Xvector1,1)+1])
ylim([nanmin(ML_results.Rho(:,1))-0.01 nanmax(ML_results.Rho(:,1)+0.01)])
set(axes1,'YTick',0.8:0.05:1.2);
%% Figure S-4: Performance in predict trials
fig2=figure();
s2=scatter(ones(size(Prediction_Performance,1),1).*6,Prediction_Performance,...
    'jitter','on', 'jitterAmount', 0.50);
s2.MarkerEdgeColor = [1 1 1];
s2.MarkerFaceAlpha= 0.5;
s2.MarkerFaceColor = [0.3 0.75 0.93];
hold on
X2 = [1 2 3 4];
e2=errorbar(X2, nanmean(Prediction_Performance_Time,1), nanstd (Prediction_Performance_Time,1));
e2.Color = 'b';
e2.Marker = 'o';
e2.MarkerFaceColor = 'b';
e2.CapSize = 16;
e2.LineStyle = '--';
e2.LineWidth = 1;
e2.MarkerSize = 10;
hold on
% for iii=1:4
% scatter(ones(size(Prediction_Performance_Time,1),1).*iii,Prediction_Performance_Time(:,iii),'jitter','on', 'jitterAmount', 0.20);
% hold on
% end
axis([0  7    0.3  1])
xlabel('$Trials$','interpreter','latex','Fontsize',18)
ylabel('$P(correct\  prediction)$','interpreter','latex','Fontsize',18)
xticks([1 2 3 4 6])
xticklabels({'Train','Early ','Middle ','Late ','All'})
yticks(0:0.1:1)
yticklabels({'0',' ',' ',' ',' ','0.5 ','0.6 ',' 0.7',' 0.8',' 0.9','1'})
hold on
l2=line([0 7],[Min_Accepted_Pre  Min_Accepted_Pre],'Color',[0.5 0.5 0.5]);
l2.LineStyle = '-.';
l2.LineWidth = 2;
l2.Color = [0.5 0.5 0.5];
text(2.5, 0.60,'$Chance\ Threshold $', 'FontSize', 14, 'Color', 'k','interpreter','latex');
box off
hold on 
e22=errorbar(6, nanmean(Prediction_Performance,1), nanstd (Prediction_Performance,1));
e22.Color = 'b';
e22.Marker = 'o';
e22.MarkerFaceColor = 'b';
e22.CapSize = 16;
e22.LineStyle = '--';
e22.LineWidth = 1;
e22.MarkerSize = 10;
hold on 
%% Figure 2-a: Choice behavior (model-free) 
figure()
p13 = StatTest_Choice_Av.p;p12 =StatTest_Choice_AvVsTa.p;p24 =StatTest_Choice_Ta.p; 
P_value = [NaN p12 p13 NaN;... 
          p12 NaN NaN p24;...
          p13 NaN NaN NaN;...
          NaN p24 NaN NaN];
bar_input = [M_Choice_Av(1) M_Choice_ta(1) ; M_Choice_Av(3) M_Choice_ta(3)]
errorbar_input = [SD_Choice_Av(1) SD_Choice_ta(1);SD_Choice_Av(3) SD_Choice_ta(3)]/4
superbar( bar_input','E',errorbar_input','P',P_value,'BarFaceColor','rk')
ylabel('$P_{gamble}$','interpreter','latex','Fontsize',20)
ylim([0 1])
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Group with averse peer', ' Group with seeker peer'})
xlabel('session\ 1/session\ 3','interpreter','latex','Fontsize',16)
title('$Choice\  Behaviour$','interpreter','latex','Fontsize',14)
%% Figure 2-b: Social shift (model-based)
fig3=figure();
subplot(1,2,1)
for tt=1:size(ML_Rho_O_av,1)
    plot(ML_Rho_O_av(tt,1:2:3),'color',[0.65 0.65 0.65],'marker','o',...
        'MarkerFaceColor',[0.65 0.65 0.65],'MarkerSize',8);
    hold on
end
poly_fit3=polyfit([ones(size(ML_Rho_O_av,1),1) ;2*ones(size(ML_Rho_O_av,1),1)]...
    ,[ML_Rho_O_av(:,1);ML_Rho_O_av(:,3)],1);
l31=refline(poly_fit3(1),poly_fit3(2));
l31.LineWidth = 4;
l31.Color = 'r';
xlim([1-0.05 2+0.05])
ylim([nanmin(ML_results.Rho(:,1))-0.05 nanmax(ML_results.Rho(:,1))+0.05])
xticks([1 2])
xticklabels({'One ','Three '})
xlabel('$Session$','interpreter','latex','Fontsize',14)
ylabel('$Risk\ Attitude\ (\rho)$','interpreter','latex','Fontsize',14)
title('$risk\  averse\  peer$','interpreter','latex','Fontsize',14)
l32=line([0 3],[1 1],'Color','k');
l32.LineStyle = '--';
subplot(1,2,2)
for tt=1:size(ML_Rho_O_ta,1)
    plot(ML_Rho_O_ta(tt,1:2:3),'color',[0.65 0.65 0.65],'marker','o',...
        'MarkerFaceColor',[0.65 0.65 0.65],'MarkerSize',8);
    hold on
end
poly_fit3=polyfit([ones(size(ML_Rho_O_ta,1),1) ;2*ones(size(ML_Rho_O_ta,1),1)],[ML_Rho_O_ta(:,1);ML_Rho_O_ta(:,3)],1);
l33=refline(poly_fit3(1),poly_fit3(2));
l33.LineWidth = 4;
l33.Color = 'r';
xlim([1-0.05 2+0.05])
ylim([nanmin(ML_results.Rho(:,1))-0.05 nanmax(ML_results.Rho(:,1))+0.05])
xticks([1 2])
xticklabels({'One ','Three '})
xlabel('$Session$','interpreter','latex','Fontsize',14)
l34=line([0 3],[1 1],'Color','k');
l34.LineStyle = '--';
%ylabel('$Risk\ Attitude\ (\rho)$','interpreter','latex','Fontsize',14)
title('$risk\  seeker\  peer$','interpreter','latex','Fontsize',14)
suptitle('Social Shift in risk attitude')
%% Figure 3-a: Contagion (for all participants)  
% violon plot: Degree of contagion across participants, n = 31
fig5=figure();
X5=Contagion;
X5 = X5(~isnan(X5));
[patch5,legend5,MX5,MED5]=violin(X5);
title('Violin Plot for Contagion','interpreter','latex')
xticks([])
ylabel('$Contagion\ (\Delta)$','interpreter','latex','FontSize',14)
ylim('auto')
hold on
l51=line([0.5 1.5],[0 0],'Color',[0.5 0.5 0.5]);
l51.LineStyle = '--';
l51.LineWidth = 1;
l51.Color = [0.5 0.5 0.5];
hold on
rng(7)
s5 = scatter(ones(size(X5,1),1),X5,'jitter','on', 'jitterAmount', 0.10);
s5.MarkerEdgeColor = 'k';
s5.MarkerFaceAlpha= 0.5;
s5.MarkerFaceColor = 'r';
legend5.String = {'mean' 'median' 'zero line'  'subjects contagion'};
legend5.Location = 'eastoutside';
patch5.LineStyle = 'none';
patch5.FaceColor = [0.8 0.8 0.8];
%% Figure 3-b
figure()
p11 = StatTest_Contagion_group1.p;p12 =StatTest_contagion_g1vsg2.p;p22 =StatTest_Contagion_group2.p; 

P_value = [p11 p22]
bar_input = [nanmean(Contagion_group1) nanmean(Contagion_group2) ]
errorbar_input = [nanstd(Contagion_group1) nanstd(Contagion_group2)]/4
superbar( bar_input','E',errorbar_input','P',P_value,'BarFaceColor','rk')
P_value = [p11 p12; p12 p22]; 
hold on
superbar( bar_input','E',errorbar_input','P',P_value,'BarFaceColor',[0.5 0.5 0.5])
ylabel('$Contagion$','interpreter','latex','Fontsize',20)
ylim([0 0.1])
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Aversive', 'Seeker'})
xlabel('$peer\ type$','interpreter','latex','Fontsize',20)


%% Figure 4-a: Correlation between social shift and social distance
figure()
fig6= plot(delta_OS,delta_SS,'o','color','k');
fig6.Color = 'k';
fig6.MarkerFaceColor = 'r';
fig6.MarkerSize = 9;
xlabel('$X:Distance\  in\  risk\  attitude\ (\rho_{peer}- \rho_{S1})$','interpreter','latex','Fontsize',14);
ylabel('$Y:Change\ of\  risk\  attitude\  (\rho_{S3} - \rho_{S1})$','interpreter','latex','Fontsize',14);
xlim ([-0.301 0.301]);
%least-square line
l61 =lsline;
l61.Color = 'r';
l61.LineWidth = 4;
poly_fit6 = polyfit(get(l61,'xdata'),get(l61,'ydata'),1);
%identity and horizental and vertical lines
l62= refline(1,0);
l62.Color = [0.5 0.5 0.5];
l62.LineStyle = '--';
l63=line([0 0],[-0.4 0.4],'Color',[0.5 0.5 0.5]);
l64=line([-0.4 0.4],[0 0],'Color',[0.5 0.5 0.5]);
% caption1 = sprintf('Y = %.2f * X + %.2f', poly_fit5(1), poly_fit5(2));
% text(0.01, -0.3,caption1, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold','interpreter','latex');
caption2 = sprintf('r = %.2f ; p = %.2e', SSOS_Corr_test.R, SSOS_Corr_test.P_value);
text(0.01, -0.25,caption2, 'FontSize', 12, 'Color', 'k','interpreter','latex');
ylim ([-0.301 0.301]);
%%
 mdl2=fitglm([ML_results.Rho(:,1) ML_results.Rho(:,2)],ML_results.Rho(:,3))
 mdl2=fitglm([ML_results.Rho(:,1) ML_results.Rho(:,2)],ML_results.Rho(:,3),'Intercept',false)
 mdl2=fitglm(ML_results.Rho(:,2)-ML_results.Rho(:,1),ML_results.Rho(:,3)-ML_results.Rho(:,1))
 lb = zeros(1,2);
 Aeq = ones(1,2);
 beq = 1;
 X1 = ML_results.Rho(:,1);
 X1 =X1(~isnan(X1));
 X2 = ML_results.Rho(:,2);
 X2 =X2(~isnan(X2));
 Y = ML_results.Rho(:,3);
 Y =Y(~isnan(Y));
 coef1 = lsqlin([X1 X2],Y,[],[],Aeq,beq,lb,ub)
 X1 = ML_results.model_free(:,1);
 X1 =X1(~isnan(X1));
 X2 = ML_results.model_free(:,2);
 X2 =X2(~isnan(X2));
 Y = ML_results.model_free(:,3);
 Y =Y(~isnan(Y));
 coef2 = lsqlin([X1 X2],Y,[],[],Aeq,beq,lb,ub)
 mdl3=fitglm([ML_results.model_free(:,1) ML_results.model_free(:,2)],...
     ML_results.model_free(:,3),'Intercept',false)