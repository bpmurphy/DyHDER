function [ HSI ] = HabitatSuitabilityIndex( SuitabilityValues,Metapop,Inputs,SiteNames,LifeStages )
% Calculate Time-Series for Habitat Suitability Indices based on User
% Prescribed Habitat Parameters (Inputs.hsi_params)

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

for t = 1:Inputs.dt:Inputs.tmax

h = numel(fieldnames(Metapop));

for i = 1:h
    
    if isempty(Inputs.hsi_params)==0
    
    for j = 1:numel(LifeStages)
        
        for k = 1:numel(Inputs.hsi_params)
            HSI_Opts(k) = SuitabilityValues.(SiteNames{i}).(Inputs.hsi_params{k}).(LifeStages{j})(t);
        end
            
        % Adjust Based on Fuzzy Logic Option
        if strcmp(Inputs.fuzzytype,'product')==1
            prod_logic = prod(HSI_Opts);
            HSI.(SiteNames{i}).(LifeStages{j})(t) = prod_logic;
        elseif strcmp(Inputs.fuzzytype,'minimum')==1
            min_logic = min(HSI_Opts);
            HSI.(SiteNames{i}).(LifeStages{j})(t) = min_logic;
        elseif strcmp(Inputs.fuzzytype,'geomean')==1
            pwr_logic = prod(HSI_Opts)^(1/numel(HSI_Opts));
            HSI.(SiteNames{i}).(LifeStages{j})(t) = pwr_logic;
        else
            error('Error: Must Define Fuzzy Logic Method')
        end
    
    end
    
    %If not habitat metrics are assigned to create HSI, then all = 1
    else
        
        for j = 1:numel(LifeStages)
            HSI.(SiteNames{i}).(LifeStages{j})(t) = 1;
        end
        
    end
 
end
    
end

end