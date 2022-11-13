clc 
clearvars
close all
Num_Trials = 35; % less than 56
U_sure =30;
beta = 5;
load('conditions.mat')
data_prob = PR(:,1);
data_mag= PR(:,2);
%% Figure S-1. Set of gambles presented
    data_choice = sign(data_prob.*data_mag-U_sure);

    idxC10 = find (data_choice==-1);
    idxC11 = find (data_choice==1);
    idxC12 = find (data_choice==0);
    FigureS1=figure()
    plot(data_prob(idxC10),data_mag(idxC10),'ob','MarkerFaceColor','b','MarkerSize',7,'MarkerEdgeColor','k')
    hold on
    plot(data_prob(idxC11),data_mag(idxC11),'or','MarkerFaceColor','r' ,'MarkerSize',7,'MarkerEdgeColor','k')
    hold on
    plot(data_prob(idxC12),data_mag(idxC12),'o','MarkerSize',7,'MarkerFaceColor','none','MarkerEdgeColor','k')
    hold on
    xlim([0 1])
    ylim([0 210])
    hold all
    x=0:0.01:1;
    y= (U_sure./x);
    plot(x,y,'-k','LineWidth',2)
    grid on
    ax=gca
    ax.GridLineStyle = '--' 
    title(['']) 
    xlabel('$Reward\ probability$','interpreter','latex','Fontsize',14)
    ylabel('$Reward\ magnitude(TT)\ $','interpreter','latex','Fontsize',14)

    