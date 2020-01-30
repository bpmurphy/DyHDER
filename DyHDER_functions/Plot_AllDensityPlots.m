function Plot_AllDensityPlots( Inputs,Record,SimName,SiteNames )
% Plot Mean (plus 5th & 95th percentiles) of Pop. Density Time-series 
% For all Sub-populations and Metapopulation across all simulations

%Note 1: Abundance and Carrying Capacity are adjusted to not include 
%the smallest stage in the matrix model.
%Note 2: Figure layout currently only configured for 7 subpopulations 

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

if Inputs.EnvStochasticity == 1
    
    %Create Variables
    K_metapop = 0;
    Compile = [];
    
    %Number of sites
    h = numel(fieldnames(Record.(SimName)));
    
    %Set Subplot Site Locations within Figure
    spar = [1,2,3,4,8,12,16];

    figure('Name','Pop. Densities: Mean + 5th & 95th Percentiles')
    set(gcf,'units','normalized','Position',[0.2 0.07 0.6 0.8])
    hold all
    
    for i = 1:h
    
    %Adjust Carrying Capacity to Include Stage 2 and Up
    K = Record.(SimName).(SiteNames{i}).K(1)*...
        sum(Record.(SimName).(SiteNames{i}).Proportion(2:end,1));
    K_metapop = K_metapop + K;
    
    %Only Pull & Aggregate Abundance for Stage 2 and Up
    for ss = 1:Inputs.SimNum
        SimName = sprintf('Simulation%d',ss);
        Compile(ss,:) = sum(Record.(SimName).(SiteNames{i}).Stages(2:end,:))/K; %Age 1+
    end
       
    %Compute Pop. Density Stats
    for t = 1:Inputs.tmax
        MeanAbund(t) = mean(Compile(:,t));
        Abund5th(t) = prctile(Compile(:,t),5);
        Abund95th(t) = prctile(Compile(:,t),95);
    end
        
   %Subplots for each Subpopulation 
        subplot(4,4,spar(i))
        hold on
        plot([0 Inputs.tmax],[1 1],'--k')
        plot(1:Inputs.tmax,MeanAbund(1:Inputs.tmax),'-','Color',[0 0 0],'LineWidth',1.5)
        plot(1:Inputs.tmax,Abund5th(1:Inputs.tmax),'-','Color',[0.5 0.5 1],'LineWidth',0.5)
        plot(1:Inputs.tmax,Abund95th(1:Inputs.tmax),'-','Color',[0.5 0.5 1],'LineWidth',0.5)
               
        text(0.03*Inputs.tmax,0.1,['Site: ',SiteNames{i}],'FontName','Arial','FontSize',8)
    
        ax = gca;
        ax.FontSize = 8;
        ax.YLim = [0 1.5];
        ax.YTick = 0:0.5:1.5;
        ax.XLim = [0 Inputs.tmax];
        ax.XTick = 0:Inputs.tmax/4:Inputs.tmax;
        ax.YAxis.MinorTick = 'on';
        ax.YAxis.MinorTickValues = 0:0.5:1.5;
        ax.YMinorGrid = 'on';
        grid on
        box on        

    end

%Set subplot location for the total Metapopulation plot
subplot(4,4,[5 6 7 9 10 11 13 14 15])

%Compile Metapopulation Results
Compile = zeros(Inputs.SimNum,Inputs.tmax);

%Only Pull & Aggregate Abundance for Stage 2 and Up
for i = 1:h
    for ss = 1:Inputs.SimNum
        SimName = sprintf('Simulation%d',ss);
        Compile(ss,:) = Compile(ss,:)+sum(Record.(SimName).(SiteNames{i}).Stages(2:end,:)); %Age 1+
    end
end

%Normalize Abundance by Carrying Capacity
Compile = Compile/K_metapop;

%Compute Pop. Density Stats for Metapopulation
    for t = 1:Inputs.tmax
        MeanPop(t) = mean(Compile(:,t));
        StdevPop(t) = std(Compile(:,t)); 
        Abund5th(t) = prctile(Compile(:,t),5);
        Abund95th(t) = prctile(Compile(:,t),95);
    end
    
%Plot Metapopulation 
    hold on
    plot([0 Inputs.tmax],[1 1],'--k')
    plot([0 Inputs.tmax],[Inputs.Nx Inputs.Nx],'--r')
    text(0.79*Inputs.tmax,1.6*Inputs.Nx,'Quasi-Extinction','Color','r','FontName','Arial','FontSize',9)
    
    plot(1:Inputs.tmax,MeanPop(1:Inputs.tmax),'-','Color',[0 0 0],'LineWidth',1.5)
    plot(1:Inputs.tmax,Abund5th(1:Inputs.tmax),'-','Color',[0.5 0.5 1],'LineWidth',1)
    plot(1:Inputs.tmax,Abund95th(1:Inputs.tmax),'-','Color',[0.5 0.5 1],'LineWidth',1)
    
    grid on
    xlabel('Model Year','FontName','Arial','FontSize',11)
    ylabel('Population Density, N/K','FontName','Arial','FontSize',11)
    set(gca,'xtick',[0:10:Inputs.tmax])
    set(gca,'FontName','Arial','FontSize',10)
    box on
    axis([0 Inputs.tmax 0 1.2])

end

end

