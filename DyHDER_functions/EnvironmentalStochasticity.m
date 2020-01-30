function [ LM_adj ] = EnvironmentalStochasticity( LM_adj,Metapop,SiteNames )
% Environmental Stochasticity: Adjust All Rates Based on Assigned Variance
% using random normal distributions

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

h = numel(fieldnames(Metapop)); %number of sites
[~,n] = size(Metapop.(SiteNames{1}).LM); %number of life-stages

for i = 1:h
    LM_adj.(SiteNames{i}) = LM_adj.(SiteNames{i}) +  ... %mean +
        Metapop.(SiteNames{1}).Variance.*randn(n,n); %stdev*randn
end

end
