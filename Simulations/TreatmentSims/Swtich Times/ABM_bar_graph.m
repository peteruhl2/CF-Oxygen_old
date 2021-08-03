%%% quick bar graph for some agent based ensembles
%%% 8/2/2021

close all

%%% lambda
x = [0.8; 0.9; 1.0; 1.1];
mean = [3.69 0.63; 8.65 2.35; 20.51 5.86; 127.01 18.1];
sd = [0.97 nan; 1.59 nan; 5.31 nan; 43.16 nan];

figure()
hold on
box on
bar(x,mean)

% er = errorbar(x,mean,sd,sd);
% er.Color = 'k';
% er.LineStyle = 'none';

xticks([0.8 0.9 1.0 1.1])
xticklabels({'80%','90%','100%','110%'})

xlabel('Oxygen flow rate (percent of normal)')
ylabel('Time to population switch (days)')
title("Days to population switch as function of \lambda")
legend('Spatially heterogenous model','Spatially homogenous model')

%%% =======================================================================
%%% dns
x2 = [1;2;3;4;5];
mean2 = [23.45 6.82; 12.08 4.08; 6.43 2.93; 2.77 3.5; 1.8 6.55];
sd2 = [3.57;2.28;1.68;0.89;0.57];

figure()
hold on
box on
bar(x2,mean2)

% er = errorbar(x2,mean2,sd2,sd2);
% er.Color = 'k';
% er.LineStyle = 'none';

xticks([1 2 3 4 5])
xticklabels({'0.5','1.0','2.0','4.0','5.0'})

xlabel('Bacterial death rate per day')
ylabel('Time to population switch (days)')
title("Days to population switch as function of d_{N}")
legend('Spatially heterogenous model','Spatially homogenous model')