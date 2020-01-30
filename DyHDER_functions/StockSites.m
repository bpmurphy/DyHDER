function [ Record ] = StockSites( t,SimName,Record,Stocking,LifeStages )
%Stock/Harvest Sites based on Prescribed Stocking/Harvesting Scenarios

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

if ismember(t,Stocking.StockYears)    
    
    stock_ind = find(Stocking.StockYears==t);
    
    for kk = 1:numel(stock_ind)
        
        %Identify Stage to Add/Remove
        stage_ind = strmatch(Stocking.StockAge(kk),LifeStages,'exact');
        
        %Add/Remove Individuals from Total Abundance Calculation
        Record.(SimName).(Stocking.StockSites{kk}).Abundance(t) = ...
            Record.(SimName).(Stocking.StockSites{kk}).Abundance(t) + Stocking.StockAbund(kk);
        
        %In case asked to harvest more than present, Abundance = 0
        if Record.(SimName).(Stocking.StockSites{kk}).Abundance(t) < 1
            Record.(SimName).(Stocking.StockSites{kk}).Abundance(t) = 0;
        end
        
        %Add/Remove Individuals from Stage Abundance Calculation
        Record.(SimName).(Stocking.StockSites{kk}).Stages(stage_ind,t) = ...
            Record.(SimName).(Stocking.StockSites{kk}).Stages(stage_ind,t) + Stocking.StockAbund(kk);
        
        %In case asked to harvest more than present, Abundance = 0
        if Record.(SimName).(Stocking.StockSites{kk}).Stages(stage_ind,t) < 1
            Record.(SimName).(Stocking.StockSites{kk}).Stages(stage_ind,t) = 0;
        end
        
        %Recalculate Stage Proportions
        Record.(SimName).(Stocking.StockSites{kk}).Proportion(:,2) = ...
            Record.(SimName).(Stocking.StockSites{kk}).Stages(:,2)/Record.(SimName).(Stocking.StockSites{kk}).Abundance(2);
        
    end
    
end

end

