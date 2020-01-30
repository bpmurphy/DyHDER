function Plot_MetapopAbundance_AllSims( Inputs,SimName,SiteNames,Population,Record,K_total,probext )
%Plot Metapopulation Abundance Time-series for All Simulations

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

h = numel(fieldnames(Record.(SimName))); %number of sites

%% PLOT: Total Abundance Through Each Simulation
figure('Name','Metapopulation Abundance: All Simulations')
set(gcf,'units','normalized','position',[0.02 0.25 0.8 0.6])

hold all
for ss = 1:Inputs.SimNum
    SimName = sprintf('Simulation%d',ss);
    plot(1:Inputs.dt:Inputs.tmax,Population.(SimName).Abundance,'-','LineWidth',0.5,...
    'Color',[0.35 0.35 0.35])
end

box on
axis([0 Inputs.tmax 0 1.25*K_total(1)])
plot([0 Inputs.tmax],K_total(1)*[Inputs.Nx Inputs.Nx],'--r')
text(1.02*Inputs.tmax,K_total(1)*Inputs.Nx,'$N_{x}$','Interpreter','latex')
hold on
plot([0 Inputs.tmax],[K_total(1) K_total(1)],'--b')
text(1.02*Inputs.tmax,K_total(1),'$N_{K}$','Interpreter','latex')
box on
grid on
xlabel('Model Year','Interpreter','latex')
ylabel('Metapopulation Abundance, N','Interpreter','latex')
set(gca,'xtick',[0:10:Inputs.tmax])
set(gca,'TickLabelInterpreter','latex','FontSize',11)

if Inputs.terminate == 1
    text(0.8*Inputs.tmax,1.28*K_total(1),...
    ['Prob. of Extinction = ', num2str(100*probext), '${\%}$'],...
    'Interpreter','latex')
end


%% ALTERNATIVE PLOT: Subpopulation Abundance Time-series (if Deterministic)

if Inputs.EnvStochasticity == 0
       
    figure('Name','Subpopulation Abundances: Deterministic')
    set(gcf,'units','normalized','position',[0.02 0.25 0.8 0.6])
    hold all
    
    rnd = round(linspace(0.65,0.95,h)*Inputs.tmax);
    for i = 1:h
        plot(1:Inputs.dt:Population.(SimName).TotalYears,Record.(SimName).(SiteNames{i}).Abundance,'-')
        text(rnd(i),Record.(SimName).(SiteNames{i}).Abundance(rnd(i)),SiteNames{i})
    end

    set(gcf,'units','normalized','position',[0.02 0.25 0.8 0.6])
    hold on
    box on
    grid on
    xlabel('Model Year','Interpreter','latex')
    ylabel('Site Abundance, n','Interpreter','latex')
    set(gca,'xtick',[0:10:Inputs.tmax])
    set(gca,'TickLabelInterpreter','latex','FontSize',11)
    
end



end

