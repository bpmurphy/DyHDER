function [ Record ] = MatrixProjection( t,Record,SimName,SiteNames,LM_adj )
%Project Adjusted Subpopulation Matrices (LM_adj) and Record Conditions
 
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
    
    %Project Matrix
    Record.(SimName).(SiteNames{i}).Stages(:,t) = ...
        LM_adj.(SiteNames{i}) * Record.(SimName).(SiteNames{i}).Stages(:,t-1);
    
    %Record Stage-based Abundances
    Record.(SimName).(SiteNames{i}).Abundance(t) = ...
        sum(Record.(SimName).(SiteNames{i}).Stages(:,t));
    %Record Stage Proportions
    Record.(SimName).(SiteNames{i}).Proportion(:,t) = ...
        Record.(SimName).(SiteNames{i}).Stages(:,t)/Record.(SimName).(SiteNames{i}).Abundance(t);
    %Record Site Pop. Density
    Record.(SimName).(SiteNames{i}).Density(t) = ...
        Record.(SimName).(SiteNames{i}).Abundance(t)/Record.(SimName).(SiteNames{i}).K(t);
    
end

end

