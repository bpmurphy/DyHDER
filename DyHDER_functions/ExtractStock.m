function [ Stocking ] = ExtractStock( InputsFolder,Filenames )
%Extract Stocking Parameters

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

ParentFolder = cd(InputsFolder);

if isempty(Filenames.('Stocking'))==0

    [stock_num,stock_txt,~] = xlsread(Filenames.('Stocking'));
    Stocking.StockSites = stock_txt(2:end,1);
    Stocking.StockYears = stock_num(:,1);
    Stocking.StockAge = stock_txt(2:end,3);
    Stocking.StockAbund = stock_num(:,3);

end

cd(ParentFolder);

end

