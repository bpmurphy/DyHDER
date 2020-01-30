function [ Filenames,ParentFolder ] = FileReadIn( InputsFolder )
%Read-in all .xlsx files of required input data

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

Filenames.('Demographics') = ls('Demographic*.xlsx');
Filenames.('Connectivity') = ls('Connectivity*.xlsx');
Filenames.('Dispersal') = ls('Dispersal*.xlsx');
Filenames.('Distance') = ls('Distance*.xlsx');
Filenames.('ScalarDistance') = ls('ScalarDistance*.xlsx');
Filenames.('Habitat') = ls('Habitat*.xlsx');
Filenames.('Stocking') = ls('Stocking*.xlsx');
Filenames.('Suitability') = ls('Suitability*.xlsx');

cd(ParentFolder) %Pass back to parent folder

end

