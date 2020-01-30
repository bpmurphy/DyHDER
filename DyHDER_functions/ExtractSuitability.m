function [ SuitabilityRelations ] = ExtractSuitability( InputsFolder,Filenames,LifeStages )
% Extract Species-specific Habitat Suitability Relations for each of the 
% Habitat Metrics prescribed in the model.  

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

[numfiles,~] = size(Filenames.Suitability);


for i = 1:numfiles
    
    %Pull Filename & Remove Spaces
    Filename = Filenames.Suitability(i,:);
    Filename = Filename(~isspace(Filename));
    
    %Pull Data
    [num,txt,raw] = xlsread(Filename);
    
    %ID Metric
    metric = txt{2,1};
    metric = metric(find(isspace(metric)==0));
    
    %Remove NaN columns if exist
    temp = any(isnan(num), 1);
    num(:, temp) = [];
    
    %Create Structure
    SuitabilityRelations.(metric).Value = num(:,1);
    
    [~,opt_col] = size(num);
    
    %Extract data either for all or individual life-stages
    if opt_col - 1 == 1
        SuitabilityRelations.(metric).All = num(:,2);
    elseif opt_col - 1 == numel(LifeStages)
        
        for jj = 1:numel(LifeStages)
            SuitabilityRelations.(metric).(LifeStages{jj}) = num(:,jj+1);
        end
        
    else
       error('Error: Suitability Input File Incorrectly Formatted')
    end
    
end


cd(ParentFolder)

end

