function [ Record ] = Eigenvalue( t,Record,SimName,LM_adj,SiteNames )
%Compute Strictly Dominant Eigenvalue of Lefkovitch Matrix

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
    
    % Calculate Eigenvalues of Leslie Matrix
    [~,Diag] = eig(LM_adj.(SiteNames{i}));    
    
    %Find Strictly dominant eigenvalue
    Eigenvalue = max(Diag(Diag>0)); 
    
    %Record Strictly Dominant Eigenvalue
    if isempty(Eigenvalue)==1
        Record.(SimName).(SiteNames{i}).Eigenvalue(t) = 0;
    else
        Record.(SimName).(SiteNames{i}).Eigenvalue(t) = Eigenvalue; 
    end

end

