function [varargout] = CompileResults( Inputs,Record,SimName,SiteNames)
% Compile Mean (and 5th & 95th percentiles) of Pop. Density 
% For all Sub-populations and Metapopulation across all simulations
% From the last n = Inputs.CompileYears of the model

% Note: Abundance and Carrying Capacity are adjusted to not include 
% the smallest stage in the matrix model.

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

   
%Create Variables
K_metapop = 0;
Compile = [];

%Number of sites
h = numel(fieldnames(Record.(SimName)));

%Number of final model years to compile
cmpyr = Inputs.CompileYears - 1;

%In case CompileYears > Model Run Years, Pull All model years
if cmpyr > Inputs.tmax
    cmpyr = Inputs.tmax - 1;
    warning('Inputs.CompileYears > Inputs.tmax: Model Results Compiled from All Years')
end

for i = 1:h

%Adjust Carrying Capacity to Include Stage 2 and Up
K = Record.(SimName).(SiteNames{i}).K(1)*...
    sum(Record.(SimName).(SiteNames{i}).Proportion(2:end,1));
K_metapop = K_metapop + K;

%Only Pull & Aggregate Abundance for Stage 2 and Up
for ss = 1:Inputs.SimNum
    SimName = sprintf('Simulation%d',ss);
    Compile(ss,:) = sum(Record.(SimName).(SiteNames{i}).Stages(2:end,:))/K; %Age 1+
end

%Compute Pop. Density Stats
for t = 1:Inputs.tmax
    MeanAbund(t) = mean(Compile(:,t));
    Abund5th(t) = prctile(Compile(:,t),5);
    Abund95th(t) = prctile(Compile(:,t),95);
end

%Record Subpopulation Results
ModelResults(i,:) = ...
    [mean(mean(Compile(:,Inputs.tmax-cmpyr:Inputs.tmax))),...
    mean(prctile(Compile(:,Inputs.tmax-cmpyr:Inputs.tmax),5)),...
    mean(prctile(Compile(:,Inputs.tmax-cmpyr:Inputs.tmax),95))];

end

%Compile Metapopulation Results
Compile = zeros(Inputs.SimNum,Inputs.tmax);

%Only Pull & Aggregate Abundance for Stage 2 and Up
for i = 1:h
    for ss = 1:Inputs.SimNum
        SimName = sprintf('Simulation%d',ss);
        Compile(ss,:) = Compile(ss,:)+...
            sum(Record.(SimName).(SiteNames{i}).Stages(2:end,:)); %Age 1+
    end
end

%Normalize Abundance by Carrying Capacity
Compile = Compile/K_metapop;

%Compute Pop. Density Stats for Metapopulation
    for t = 1:Inputs.tmax
        MeanPop(t) = mean(Compile(:,t));
        StdevPop(t) = std(Compile(:,t)); 
        Abund5th(t) = prctile(Compile(:,t),5);
        Abund95th(t) = prctile(Compile(:,t),95);
    end
    
%Record Metapopulation Results
ModelResults = [ModelResults;...
        mean(mean(Compile(:,Inputs.tmax-cmpyr:Inputs.tmax))),...
        mean(prctile(Compile(:,Inputs.tmax-cmpyr:Inputs.tmax),5)),...
        mean(prctile(Compile(:,Inputs.tmax-cmpyr:Inputs.tmax),95))];
  
%Write Output Table of Pop. Density Results
T = table(ModelResults(:,1),ModelResults(:,2),ModelResults(:,3),...
    'RowNames',[SiteNames,{'Total'}],'VariableNames',{'Mean','P5th','P95th'});
varargout{1} = T;
   

end