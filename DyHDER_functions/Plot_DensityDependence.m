function Plot_DensityDependence( Inputs,Metapop,LifeStages )
% Function to Plot Density Dependent Stock-Recruitment
% Default is to Plot Values from First Listed Site

%% Copyright & Licensing

% Dynamic Habitat Disturbance & Ecological Resilience (DyHDER)
% Copyright (C) 2020 Brendan P. Murphy

% This file is part of DyHDER.

% DyHDER is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
 
% DyHDER is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
 
% You should have received a copy of the GNU General Public License
% along with DyHDER. If not, see <https://www.gnu.org/licenses/>.

%%

%Pull Site Names
sites = fieldnames(Metapop); %site names
site_num = 1; %Pull First Site Listed

[~,n] = size(Metapop.(sites{site_num}).LM); %number of life-stages

%Inputs
f = Metapop.(sites{site_num}).OffSpring(n); %Site 1 Offspring for Stage=n
s0 = Metapop.(sites{site_num}).EggSurvival(2); %Max Offspring Survival
s1 = Metapop.(sites{site_num}).EggSurvival(1); %Offspring Survival at n/K = 1
K = Metapop.(sites{site_num}).K;
Prop_n = Metapop.(sites{site_num}).Proportion(n);
Stock = K*Prop_n;
Beta = Metapop.(sites{site_num}).Beta(1); %fitting parameter
nk = 0:0.01:1.5; % Population Densities plot range

%Set Up Figure
figure('Name','Density Dependence Parameters')
set(gcf,'units','normalized','position',[0.5 0.1 0.4 0.75])
hold on


%Ricker Model
if Inputs.DD_model == 1

    %Fecundity
    F_plot = f*s0*exp(-Beta*nk);

    %Plot
    subplot(2,1,1)
    plot(nk,F_plot,'LineWidth',1.5)
    hold on
    plot(nk,f*s1*ones(size(nk)),'k--')
    set(gca,'ylim',[0 1.25*max(F_plot)])
    yL = ylim;
    plot([1,1],yL,'k--')
    title(['Ricker Model: ',strrep(LifeStages{n},'_',' ')],'Interpreter','latex')
    xlabel('Population Density','Interpreter','latex')
    ylabel('Offspring per Female','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')

    subplot(2,1,2)
    x = linspace(0,1.5*Stock,numel(nk));
    plot(x,x.*F_plot,'LineWidth',1.5)
    hold on
    set(gca,'xlim',[0 max(x)],'ylim',[0 1.25*max(x.*F_plot)])
    yL = ylim;
    plot([Stock,Stock],yL,'k--')
    xlabel(['Number of ``',strrep(LifeStages{n},'_',' '),'"'],'Interpreter','latex')
    ylabel(['Number of Recruits'],'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')

%Beverton-Holt Model
elseif Inputs.DD_model == 2 

    F_plot = f*s0./(1+Beta*nk);

    %Plot
    subplot(2,1,1)
    plot(nk,F_plot,'LineWidth',1.5)
    hold on
    plot(nk,f*s1*ones(size(nk)),'k--')
    set(gca,'ylim',[0 1.25*max(F_plot)])
    yL = ylim;
    plot([1,1],yL,'k--')
    title(['Beverton-Holt: ',strrep(LifeStages{n},'_',' ')],'Interpreter','latex')
    xlabel('Population Density','Interpreter','latex')
    ylabel('Offspring per Female','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    grid on

    subplot(2,1,2)
    x = linspace(0,1.5*Stock,numel(nk));
    plot(x,x.*F_plot,'LineWidth',1.5)
    hold on
    set(gca,'xlim',[0 max(x)],'ylim',[0 1.25*max(x.*F_plot)])
    yL = ylim;
    plot([Stock,Stock],yL,'k--')
    xlabel(['Number of ``',strrep(LifeStages{n},'_',' '),'"'],'Interpreter','latex')
    ylabel(['Number of Recruits'],'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    grid on

end

subplot(2,1,1)
if Inputs.DD_model == 1
    text(1.15,1.05*F_plot(1),['$\beta$ = ', num2str(Beta,3)],'FontSize',14,'Interpreter','latex')
    text(0.05,(f*s1)*1.18,strcat('$f_{',strrep(LifeStages{n},'_',''),'}*\phi_{0}(1)$'),'Interpreter','latex')
elseif Inputs.DD_model == 2
    text(1.15,1.05*F_plot(1),['$\beta$ = ', num2str(Beta,3)],'FontSize',14,'Interpreter','latex')
    text(0.05,(f*s1)*1.18,strcat('$f_{',strrep(LifeStages{n},'_',''),'}*\phi_{0}(1)$'),'Interpreter','latex')
end


end

