%% Housekeeping
clc; clear all; close all;
% Random seed for replication
rng(10)
% Given parameters
par.alpha = 0.67;
par.beta  = 0.2;
par.sigma = -0.9;
par.Sigma = [1,par.sigma;par.sigma,1];
par.mu    = [0,0];
par.C     = 1.5;
N         = 1e03;
% Generate model 
U   = mvnrnd(par.mu,par.Sigma,N);
Y1  = par.alpha + par.beta + U(:,1);
Y0  = par.alpha + U(:,2);
D   = (Y1-Y0-par.C>=0);
Y   = Y1.*D + (1-D).*Y0;

%% Question (a)


% MTE
ud          = normcdf(diff(U,1,2)/3.8);
idx         = floor(round(ud,2)*100);
MTEm        = accumarray(idx,Y1-Y0,[],@mean,NaN);
x           = accumarray(idx,idx/100,[],@mean,NaN);
MTEf        = @(u) par.beta - 2*(1-par.sigma)*norminv(u);

% Import colors
global c
myColors();

% Plot 
figure;
hold on 
plot(linspace(0,1,100),MTEf(linspace(0,1,100)),'Color',c.maroon)
scatter(x,MTEm,40,'MarkerEdgeColor','k','MarkerFaceColor',c.nvyBlue)
xlabel('$u_D$')
ylabel('MTE$(u_D)$')
legend('Analitic Form','From Simulation','location','northeast','box','off')
export_fig('Figures/MTEcompare','-pdf','-transparent'); 


%% Question (b-c)

% Vector of correlation
sigma = -1:0.1:1;
% Initialize treatments
ATE  = NaN(length(sigma),1);
TT   = NaN(length(sigma),1);
TUT  = NaN(length(sigma),1);
LATE = NaN(length(sigma),1);
PRTE = NaN(length(sigma),1);
% Intialize figure
mte_fig = figure;
hold on 
% Begin loop
for j=1:length(sigma)
    % Update sigma
    par.sigma = sigma(j);
    % Update matrix 
    par.Sigma = [1,par.sigma;par.sigma,1];
    % Generate model 
    U   = mvnrnd(par.mu,par.Sigma,N);
    Y1  = par.alpha + par.beta + U(:,1);
    Y0  = par.alpha + U(:,2);
    D   = (Y1-Y0-par.C>=0);
    Y   = Y1.*D + (1-D).*Y0;
    % Treatment effects
    ATE(j)  = mean(Y1-Y0);
    TT(j)   = mean(Y1(D)-Y0(D));
    TUT(j)  = mean(Y1(~D)-Y0(~D));
    % New decision for LATE
    Dn  = (Y1-Y0-(par.C-0.5)>=0);
    % Estimate LATE
    compl   = ~D & Dn;
    LATE(j) = mean(Y1(compl)-Y0(compl));
    % New decision for PRTE
    Dn  = (Y1-Y0-(par.C-1)>=0);
    Yn  = Y1.*Dn + (1-Dn).*Y0;
    % Estimate PRTE
    PRTE(j) = mean(Yn) - mean(Y);
    % Get MTE
    MTEf    = @(u) par.beta - 2*(1-par.sigma)*norminv(u);
    x       = linspace(0,1,100);
    % Add to figure (only a few to avoid over crowding)
    if(mod(j-1,4)==0)
        plot(x,MTEf(x))
    end
end
% Calculate index for plot
idx = mod((1:length(sigma))-1,4)==0;

% Finish plot 
xlabel('$u_D$')
ylabel('MTE$(u_D)$')
legend(strcat('$\rho=$',string([sigma(idx)])),'location','northeast',...
      'box','off','NumColumns',2)
export_fig('Figures/MTEgrid','-pdf','-transparent'); 

% Export treatment effects as a table 
latex.data = [ATE,TT,TUT,LATE,PRTE];
latex.tableCaption = 'Treatment effects';
latex.tableLabel = 'te';
latex.tablePositioning = 'htb';
latex.tableColumnAlignment = 'c';
latex.tableBorders = 0;
latex.dataFormat = {'%2.2f'};
latex.tableColLabels = {'ATE','TT','TUT','LATE','PRTE'};
latex.tableRowLabels = cellstr(strcat('$\rho=$',string(sigma)));
latex = latexTable2(latex);
dlmcell(strcat('Tables/treatEffect.tex'),latex)



