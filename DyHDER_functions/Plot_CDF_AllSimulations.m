function Plot_CDF_AllSimulations( Inputs,Population,K_total )
% CDF Plot of Mean Metapopulation Densities From All Years in 
% All Simulations, as well as 5th & 95th percentiles

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

%% Compile Abundance Data

    for ss = 1:Inputs.SimNum
        SimName = sprintf('Simulation%d',ss);
        Compile(ss,:) = Population.(SimName).Abundance;
    end
    
    for t = 1:Inputs.tmax
        MeanPop(t) = mean(Compile(:,t)); 
        Pop5th(t) = prctile(Compile(:,t),5); 
        Pop95th(t) = prctile(Compile(:,t),95); 
    end
    
    %Normalize to Carrying Capacity
    MeanPop = MeanPop./K_total;
    Pop5th = Pop5th./K_total;
    Pop95th = Pop95th./K_total;

    
%% CDF of Metapopulation Density

sortedVector1 = sort(MeanPop(1:Inputs.tmax));
indexOfValueChange1 = 1:numel(sortedVector1);
relativeCounts1 = (0:length(sortedVector1)-1)/(length(sortedVector1)-1);

sortedVector2 = sort(Pop5th(1:Inputs.tmax));
indexOfValueChange2 = 1:numel(sortedVector2);
relativeCounts2 = (0:length(sortedVector2)-1)/(length(sortedVector2)-1);

sortedVector3 = sort(Pop95th(1:Inputs.tmax));
indexOfValueChange3 = 1:numel(sortedVector3);
relativeCounts3 = (0:length(sortedVector3)-1)/(length(sortedVector3)-1);

figure('Name','CDF of Metapopulation Density')
hold on
p1=plot(sortedVector1(indexOfValueChange1),relativeCounts1(indexOfValueChange1),'-','Color',[0 0 0],'LineWidth',1.5);
p2 = patch([sortedVector2(indexOfValueChange2),...
    fliplr(sortedVector3(indexOfValueChange3))],...
    [relativeCounts2(indexOfValueChange2),...
    fliplr(relativeCounts3(indexOfValueChange3))],1);
set(p2,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.2,...
    'EdgeColor','k','EdgeAlpha',0)

p3=plot([1 1],[0 1],'--b','LineWidth',0.75);
p4=plot([Inputs.Nx Inputs.Nx],[0 1],'--r','LineWidth',0.75);

leg= legend([p1,p2,p3,p4],...
    {'Mean','$5^{th}$ to $95^{th}$ Percentile',...
    'Carrying Capacity (K)','Quasi-Extinction'});
set(leg,'Location','NorthWest','interpreter','latex')

xx = xlim;
set(gca, 'xlim',[0 max(xx)],'xtick',0:0.25:2,...
    'ylim',[0 1],'ytick',0:0.25:1)
set(gca,'TickLabelInterpreter','latex')
grid on
xlabel('Population Density','Interpreter','latex')
ylabel('Fraction Less Than','Interpreter','latex')
box on
axis square


end

