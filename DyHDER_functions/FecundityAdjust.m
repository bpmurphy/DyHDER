function [ LM_adj ] = FecundityAdjust( t,LM_adj,HabitatOptimality,Metapops,Inputs,SiteNames,LifeStages )
%Adjust Vital Rates based on Prescribed Habitat Adjustment Parameters
%(i.e. Inputs.surv_adj, Inputs.adv_adj)

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

h = numel(fieldnames(HabitatOptimality));

for i = 1:h
                  
    if isempty(Inputs.fec_adj)==0
    
    for j = 2:numel(LifeStages)
           
        for k = 1:numel(Inputs.fec_adj)
            Fec_Opts(k) = HabitatOptimality.(SiteNames{i}).(Inputs.fec_adj{k}).(LifeStages{j})(t);
        end
    
        fec_rate = Metapops.(SiteNames{i}).LM(1,j);
        
        % Adjust Based on Fuzzy Logic Option
        if strcmp(Inputs.fuzzytype,'product')==1
            prod_logic = prod(Fec_Opts);
            LM_adj.(SiteNames{i})(1,j) = fec_rate * prod_logic;
        elseif strcmp(Inputs.fuzzytype,'minimum')==1
            min_logic = min(Fec_Opts);
            LM_adj.(SiteNames{i})(1,j) = fec_rate * min_logic;
        elseif strcmp(Inputs.fuzzytype,'geomean')==1
            pwr_logic = prod(Fec_Opts)^(1/numel(Fec_Opts));
            LM_adj.(SiteNames{i})(1,j) = fec_rate * pwr_logic;
        else
            error('Error: Must Define Fuzzy Logic Method')
        end

    end
    
    else %No Adjustment
        
        for j = 2:numel(LifeStages)
            LM_adj.(SiteNames{i})(1,j) = Metapops.(SiteNames{i}).LM(1,j);
        end
        
    end
    
end
    

end