function [ HabitatMetrics ] = ExtractHabitat( InputsFolder,Filenames,SiteNames )
% Extract Time-series of Habitat Metrics
% First column from input folder is unused within this script, but it 
% should = t and be the same number of and represent the same timseteps as  
% prescribed in model input parameters i.e., 1:dt:tmax

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

ParentFolder = cd(InputsFolder); %Set parent folder & open subfolder

    for j = 1:numel(SiteNames)

    [hab_metrics,hab_params,~] = xlsread(Filenames.Habitat,SiteNames{j});

        for i = 2:numel(hab_params)
            param = hab_params{i};
            hab_params{i}=param(find(isspace(param)==0));
            HabitatMetrics.(SiteNames{j}).(hab_params{i}) = hab_metrics(:,i);
        end

    end

cd(ParentFolder)
    
end

