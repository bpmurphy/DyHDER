function Plot_PreAndPostHistograms( Inputs,Population,K_total,binsize )
% Plot histograms of metapopulation abundances from all simulations
% between pre- and post-disturbance years prescribed by user

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
    
    binedges = 0:binsize:K_total(1)+2*binsize;
    post = Inputs.PostDisturbance;
    pre = Inputs.PreDisturbance;

    for ss = 1:Inputs.SimNum
        SimName = sprintf('Simulation%d',ss);
        PreCat(ss) = Population.(SimName).Abundance(pre);
        PostCat(ss) = Population.(SimName).Abundance(post);
    end

    figure('Name','Pre- And Post-Disturbance Histograms')
    histogram(PreCat,binedges,'Normalization','probability')
    hold on
    histogram(PostCat,binedges,'Normalization','probability')
    legend('Pre-Disturbance','Post-Disturbance')
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('Metapopulation Abundance','interpreter','latex')
    ylabel('Frequency','interpreter','latex')
    set(gca,'Ygrid','on')
    leg = legend('Location','NorthWest');
    set(leg,'Interpreter','latex','FontSize',11)
    
end


end

