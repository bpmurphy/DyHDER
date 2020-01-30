function [ Record ] = AnnualizedGrowthRate( Record,SimName,SiteNames )
%Calculate Annualized Growth Rate for the Metapopulation

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

for i = 1:h
    
    GrowthRate = Record.(SimName).(SiteNames{i}).Abundance(2:end)./ ...
        Record.(SimName).(SiteNames{i}).Abundance(1:end-1);
    
    Record.(SimName).(SiteNames{i}).GrowthRate = GrowthRate;

end

