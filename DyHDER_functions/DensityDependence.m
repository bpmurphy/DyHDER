function [ LM_adj,Record ] = DensityDependence( t,Inputs,Record,SimName,Metapop,SiteNames,LM_adj )
% Further adjust habitat-adjusted reproduction rates in each matrix based  
% on local density using the selected density-dependent model

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
[~,n] = size(Record.(SimName).(SiteNames{1}).LM); %number of life-stages

%Calculate for each site
for i = 1:h
    
    %Find reproductive life-stages
    Index = find(Record.(SimName).(SiteNames{i}).OffSpring>0);
    
    %Adjust Offspring Survival: Habitat-dependent Fecundity Adjustment
    Adj_Factor = LM_adj.(SiteNames{i})(1,n)./Metapop.(SiteNames{i}).LM(1,n);
    s1 = Adj_Factor*Metapop.(SiteNames{i}).EggSurvival(1);
    s0 = Adj_Factor*Metapop.(SiteNames{i}).EggSurvival(2);
    
    %Call-in Variables
    K = Record.(SimName).(SiteNames{i}).K(t-1);
    Abund = Record.(SimName).(SiteNames{i}).Abundance(t-1);
    Density = Record.(SimName).(SiteNames{i}).Density(t-1);
    OffSpring = Metapop.(SiteNames{i}).OffSpring;
    
    for nn = min(Index):max(Index)
        
        % Ricker Model
        if Inputs.DD_model == 1 
            
            %Recalculate Best-fit Beta & Adjust Reproduction Rate
            BetaFec = log(s0/s1);
            LM_adj.(SiteNames{i})(1,nn) = s0*OffSpring(nn)*exp(-BetaFec*Density);
                        
        % Beverton-Holt Model
        elseif Inputs.DD_model == 2
            
            %Recalculate Best-fit Beta & Adjust Reproduction Rate
            BetaFec = (s0/s1)-1;
            LM_adj.(SiteNames{i})(1,nn) = (s0*OffSpring(nn))/(1+(BetaFec*Density));
                                    
        end
    
    end
    
    Record.(SimName).(SiteNames{i}).Beta(t)=BetaFec;

end

