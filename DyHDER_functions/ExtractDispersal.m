function [ DispersalMetrics ] = ExtractDispersal_Flexible2( InputsFolder,Filenames,SiteNames )
%Read in Migration Parameters

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

%% Connectivity Matrix 
% Elements of c = 0-1
% Matrix dimensions = h x h

% Read in Connectivity Matrix Data
[conn_values,conn_names,~] = xlsread(Filenames.('Connectivity'));
ConnNameOrder = conn_names(2,3:end);
[m,n]=size(conn_values);

% Determine Site Order Relative to 'SiteNames' Order from Demographic Data
% and Rearrange Indexing
for i = 1:numel(SiteNames)
    SiteFind = strfind(ConnNameOrder,SiteNames{i});
    IndexIn(i) = find(not(cellfun('isempty',SiteFind)));
end
IndexMatConn = sortrows([(1:length(IndexIn))',IndexIn'],2);
IndexConn = IndexMatConn(:,1);

% Write Reordered Connectivity Matrix
Connect = [];
for i = 1:n
    for j = 1:m
        row=IndexConn(j);
        col=IndexConn(i);
        Connect(row,col)=conn_values(j,i);
    end
end

DispersalMetrics.Connectivity = Connect;


%% Site-to-site Distance Matrix 
% Elements of d in km
% Matrix dimensions = h x h

% Read in Distance Matrix Data
[dist_values,dist_names,~] = xlsread(Filenames.('Distance'));
DistNameOrder = dist_names(2,3:end);
[m,n]=size(dist_values);

% Determine Site Order Relative to 'SiteNames' Order from Demographic Data
% and Rearrange Indexing
for i = 1:numel(SiteNames)
    SiteFind = strfind(DistNameOrder,SiteNames{i});
    IndexIn(i) = find(not(cellfun('isempty',SiteFind)));
end
IndexMatDist = sortrows([(1:length(IndexIn))',IndexIn'],2);
IndexDist = IndexMatDist(:,1);

% Write Reordered Distance Matrix
Distance = [];
for i = 1:n
    for j = 1:m
        row=IndexDist(j);
        col=IndexDist(i);
        Distance(row,col)=dist_values(j,i);
    end
end

DispersalMetrics.Distance = Distance;


%% Life-stage Dependent Dispersal Probabilities 
% Elements of q = 0-1
% Array length = number of life-stages (ordered from youngest to oldest)
DispersalMetrics.Dispersal = xlsread(Filenames.('Dispersal'));

%% Scaling Lengths for Migration Distance for Each Life-stage
% Elements in km
% Array length = number of life-stages (ordered from youngest to oldest)
DispersalMetrics.DistanceScalar = xlsread(Filenames.('ScalarDistance'));


cd(ParentFolder)

end

