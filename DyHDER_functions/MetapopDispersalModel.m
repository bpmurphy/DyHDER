function [ Record ] = MetapopDispersalModel( t, DispersalMetrics, HSI, Record, SimName, SiteNames, LifeStages )
% Metapopulation Dispersal Model from Murphy et al., (20xx) 
% Individuals are moved between subpopulations using a set of stepwise
% probabilistic functions. Individuals emigrate (leave & do not return) 
% based on a movement-independent function that considers habitat quality
% and population density. Individuals are then distributed (immigrate) to 
% other sites using a movement function that considers habitat suitability,
% population density, site-to-site distance, and directionally-
% dependent connectivity. 

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

% For-loop for each life stage
for j = 1:n

    % Construct Stage-Dependent Markov Transition Matrix
    [ M ] = CreateMovementMatrix( j, t, DispersalMetrics, HSI, Record, SimName, SiteNames, LifeStages );
    
    % Pull Stage Abundance from Each Site
    for i = 1:h
        Stage_All(i) = Record.(SimName).(SiteNames{i}).Stages(j,t);
    end

    % Redistribute given Life-Stage 
    Stage_All = Stage_All * M; 
    
    % Record Post-Dispersal Stage Abundance
    for i = 1:h
        Record.(SimName).(SiteNames{i}).Stages(j,t) = Stage_All(i);
    end

end

% Record Post-Dispersal Site Abundances
Record.(SimName).(SiteNames{i}).Abundance(t) = ...
        sum(Record.(SimName).(SiteNames{i}).Stages(:,t));

end

