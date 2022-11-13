clc 
clearvars
close all
Scenario=1; %Simulation scenario
Num_Subject =100; % number of subjects
%Num_Subject =1;
Num_Trials = 35; % less than 56
% Initial value for parameters and matrix
Tom2USD =4.2;
U_sure =30/Tom2USD;
beta = 5;
P_gamble = zeros(70,1);
load('conditions.mat')
bunch_of_risky_conditions = PR;
for j=1:70
    Rho(j,1)= 0.65+(j-1)*0.01
    for i =1:Num_Subject
        rng shuffle
        % design Simulation blocks
        new_index = randperm(size(bunch_of_risky_conditions,1));
        data = bunch_of_risky_conditions(new_index,:);
        data_prob = data(:,1);
        jitter = randi([-1,1],size(data,1),1);
        data_mag= data(:,2)/4.2;
        Rho_All =Rho(j,1).*ones(Num_Trials,1);
        U_risk= data_prob.*(data_mag.^(Rho_All));
        F = U_risk - U_sure;
        S = 1./(1.+exp(-beta*F));
        data_choice= binornd(1,S);
        P_gamble(j,1) = ((i-1)*P_gamble(j,1) + mean(data_choice))/i;
       
    end
%      [Rho_hat(j,1),Beta_hat(j,1),model_free_hat(j,1),NLL_hat(j,1)]...
%         = ML_fitting([data_choice data_mag data_prob NaN(35,1) U_sure*ones(35,1) ones(35,1)])
end
%%
%Table_S1= [Rho(17:4:56) Rho_hat(17:4:56)];

%%



%%
% Figure S-2
figure()
%plot(Rho,P_gamble)
xticks(0.5:0.05:1.5)
grid on
figure()
data_plot = [Rho P_gamble];
[logitCoef,dev] = glmfit(Rho,P_gamble,'binomial','logit');
logitFit = glmval(logitCoef,Rho,'logit');
plot(Rho,logitFit,'-','LineWidth',4,'Color',[1 0.4 0.4]);
hold on
plot(Rho,P_gamble,'r*','MarkerSize',6,'MarkerEdgeColor',[1 0.4 0]);
hold on
% plot(Rho,P_gamble,'bx', Rho,logitFit,'-','LineWidth',2,'Color',[0.5 0.5 0.5]);
% hold on
xlabel('$Risk\ Attitude\ (\rho)$','interpreter','latex','Fontsize',14);
ylabel('$P_{gamble}$','interpreter','latex','Fontsize',14);
xlim([min(Rho)  max(Rho)])
xticks(min(Rho):0.05:max(Rho))
grid on
ax=gca
ax.GridLineStyle = '--' 
hold on
l1=line([min(Rho) max(Rho)],[0.9 0.9],'Color',[0.25 0.25 0.25]);
l1.LineStyle = '--';
l1.LineWidth = 1;
hold on
l2=line([min(Rho) max(Rho)],[0.1 0.1],'Color',[0.25 0.25 0.25]);
l2.LineStyle = '--';
l2.LineWidth = 1;
%     if (Num_Subject ==1)
%         plot(data_prob,data_mag,'o')
%         xlim([0 1])
%         ylim([0 50])
%         hold all
%         x=0:0.01:1;
%         y= U_sure./x;
%         plot(x,y)
%     end
