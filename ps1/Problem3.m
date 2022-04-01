%% Problem 3 - Problem set 1
% Authors: Isaac Norwitch
%          Jose Quintero
%          Alex Weinberg
% Empirical Analysis 3 

% Housekeeping
clc; clear all; close all;
% Import colors
global c 
myColors();

%% Monte Carlo Simulations

% Define paramters 
N       = 1e04;
X       = [ones(N,1),randn(N,1)];
beta    = [2,3]';
sigma2  = 4;

% Simulate error termns
U = randn(N,1)*sqrt(sigma2);

% Compute outcome variable
Y = X*beta + U;

% Get estimates for beta
betaE   = (X'*X)\X'*Y;

% Calculate the estandar errors for beta
VE      = N*inv(X'*X)*var(Y-X*betaE);
betaSE  = sqrt(diag(VE)/N);

% Sanity check: compare to Matlab OLS function
lm1 = fitlm(X(:,2),Y);
latex.data = table2array(lm1.Coefficients);
latex.tableCaption = 'OLS Estimates';
latex.tableLabel = 'ea3:ps1:q3a:tab1';
latex.tablePositioning = 'htb';
latex.tableColumnAlignment = 'c';
latex.tableBorders = 0;
latex.dataFormat = {'%2.3f'};
latex.tableColLabels = lm1.Coefficients.Properties.VariableNames;
latex.tableRowLabels = {'$\beta_0$','$\beta_1$'};
latex = latexTable2(latex);
dlmcell(strcat('Tables/q3a_tab1.tex'),latex)

% Clean auxiliary variables
clear ans latex

% Montecarlo Paramters
S       = 1e04;
betaM   = NaN(2,S);

% Run simulations
for s =1:S
    % Simulate error termns
    U = randn(N,1)*sqrt(sigma2);
    % Compute outcome variable
    Y = X*beta + U;
    % Get estimates for beta
    betaM(:,s)   = (X'*X)\X'*Y;
end


% Calculate the variance and estimates
betaE2      = mean(betaM,2);
betaSE2     = sqrt(mean(betaM.^2,2)-mean(betaM,2).^2);

% Export results
latex.data = [mean(betaM,2),betaSE2];
latex.tableCaption = 'Montecarlo Estimates';
latex.tableLabel = 'ea3:ps1:q3a:tab2';
latex.tablePositioning = 'htb';
latex.tableColumnAlignment = 'c';
latex.tableBorders = 0;
latex.dataFormat = {'%2.3f'};
latex.tableColLabels = {'Estimate','SE'};
latex.tableRowLabels = {'$\beta_0$','$\beta_1$'};
latex = latexTable2(latex);
dlmcell(strcat('Tables/q3a_tab2.tex'),latex)

% Begin Figure
figure;
hold on
% Histogram
H = histogram(betaM(1,:),'Normalization','probability');
H.FaceColor = c.brBlue;
H.FaceAlpha = 1;
% Estimated mean
plot(betaE2(1)*[1,1],ylim(),'Color',c.maroon)
% Plot make-up
xlabel('$\hat{\beta}_0$')
ylabel('Probability')
xlim([1.9,2.1])
legend('Montecarlo Simulation','Estimated Mean','Location','northeast','box','off')
% Save figure
export_fig('Figures/p3qa','-pdf','-transparent')


%% Non-parametric Bootstrap

clearvars -except c

% Set up problem structure
N   = 1e04;
Y1  = 5 + randn(N,1);
Y0  = 2 + randn(N,1);
D   = (rand(N,1)>=0.5);
Y   = Y1.*D + (1-D).*Y0;

% Estimate the expected difference of Y1 and Y0
lm1 = fitlm(D,Y);
latex.data = table2array(lm1.Coefficients);
latex.tableCaption = 'OLS Estimates';
latex.tableLabel = 'ea3:ps1:q3b:tab1';
latex.tablePositioning = 'htb';
latex.tableColumnAlignment = 'c';
latex.tableBorders = 0;
latex.dataFormat = {'%2.3f'};
latex.tableColLabels = lm1.Coefficients.Properties.VariableNames;
latex.tableRowLabels = {'$\beta_0$','$\beta_1$'};
latex = latexTable2(latex);
dlmcell(strcat('Tables/q3b_tab1.tex'),latex)
clear latex

% Montecarlo Paramters
S       = 1e04;
betaM   = NaN(2,S);

% Run simulations
for s =1:S
    % Sample from observations
    data = datasample([Y,D],N);
    X    = [ones(N,1),data(:,2)];
    Ys   = data(:,1);
    % Get estimates for beta
    betaM(:,s)   = (X'*X)\X'*Ys;
end

% Calculate the variance and estimates
betaE2      = mean(betaM,2);
betaSE2     = sqrt(mean(betaM.^2,2)-mean(betaM,2).^2);

% Export results
latex.data = [mean(betaM,2),betaSE2];
latex.tableCaption = 'Bootstrap Estimates';
latex.tableLabel = 'ea3:ps1:q3b:tab2';
latex.tablePositioning = 'htb';
latex.tableColumnAlignment = 'c';
latex.tableBorders = 0;
latex.dataFormat = {'%2.3f'};
latex.tableColLabels = {'Estimate','SE'};
latex.tableRowLabels = {'$\beta_0$','$\beta_1$'};
latex = latexTable2(latex);
dlmcell(strcat('Tables/q3b_tab2.tex'),latex)


% Begin Figure
figure;
hold on
% Histogram
H = histogram(betaM(2,:),'Normalization','probability');
H.FaceColor = c.brBlue;
H.FaceAlpha = 1;
% Estimated mean
plot(betaE2(2)*[1,1],ylim(),'Color',c.maroon)
% Plot make-up
xlabel('$\hat{\beta}_1$')
ylabel('Probability')
xlim([2.9,3.1])
legend('Bootstrap Simulation','Estimated Mean','Location','northeast','box','off')
% Save figure
export_fig('Figures/p3qb','-pdf','-transparent')



