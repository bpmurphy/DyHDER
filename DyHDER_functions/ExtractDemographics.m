function [ Metapop,Demographics,LifeStages,SiteNames ] = ...
    ExtractDemographics( InputsFolder, Filenames, LifeStages )
%Extract Demographics Parameters from Input .xlsx files


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

%Check File - Must Have 6 Sheets
[~,sheetname] = xlsfinfo(Filenames.('Demographics')); 
if size(sheetname,2) ~= 6
    error('Incorrect Number of Sheets from Input File')
end

%Number of Life-Stages
n = numel(LifeStages);

%Read File
[surv,surv_txt,~] = xlsread(Filenames.('Demographics'),'Survival');

%Site & Stage Survival Estimates
SurvSites = surv_txt(1,2:end);
for i = 1:numel(SurvSites)
    for j = 1:numel(LifeStages)
    Demographics.Survival.(SurvSites{i}).(LifeStages{j})=surv(j,i);
    end
end

%Site & Stage Transition Rates
[trans,trans_txt,~] = xlsread(Filenames.('Demographics'),'Transition');
TransSites = trans_txt(1,2:end);
for i = 1:numel(TransSites)
    for j = 1:numel(LifeStages)-1 %Oldest stage can't have transition
    Demographics.Transition.(TransSites{i}).(LifeStages{j})=trans(j,i);
    end
end

%Site & Stage Reproduction Rates
[fec,fec_txt,~] = xlsread(Filenames.('Demographics'),'Reproduction');
FecSites = fec_txt(1,2:end);
for i = 1:numel(FecSites)
    for j = 1:numel(LifeStages)
    Demographics.Reproduction.(FecSites{i}).(LifeStages{j})=fec(j,i);
    end
end

%Temporal Variance of Vital Rates
[var,var_txt,~] = xlsread(Filenames.('Demographics'),'Variance');
VarSites = var_txt(1,2:end);
for i = 1:numel(VarSites)
    for j = 1:numel(LifeStages)
    Demographics.Variance.(VarSites{i}).(LifeStages{j})=var(j,i);
    Demographics.FecundityVariance.(VarSites{i}).(LifeStages{j})=var(j+n,i);
    end
end

%Site Carrying Capacity
[K,K_txt,~] = xlsread(Filenames.('Demographics'),'CarryingCapacity');
KSites = K_txt(1,2:end);
for i = 1:numel(KSites)
    Demographics.CarryingCapacity.(KSites{i})=K(:,i);
end

%Site Offspring Survival Rates
[eggsurv,eggsurv_txt,~] = xlsread(Filenames.('Demographics'),'Offspring');
eggsurvSites = eggsurv_txt(1,2:end);
for i = 1:numel(eggsurvSites)
    Demographics.EggSurvival.(eggsurvSites{i})=eggsurv(:,i);
end

SiteNames = surv_txt(1,2:end);
h = numel(SiteNames);

%% Partition Survival between Survive & Stay and Survive & Transition

for i = 1:h
    for j = 1:n
        survival(j) = Demographics.Survival.(SiteNames{i}).(LifeStages{j}); 
        if j < n
        transition(j) = Demographics.Transition.(SiteNames{i}).(LifeStages{j}); 
        end
        fecundity(j) = Demographics.Reproduction.(SiteNames{i}).(LifeStages{j}); 
    end
    
    survstay = [survival(1:end-1).*(1-transition) survival(end)];
    survtrans = survival(1:end-1).*(transition);

%% Build Lefkovitch Matrix

    LM = diag(survstay);
    index = diag(survtrans,-1)~=0;
    LM(index) = survtrans;
    
    LM(1,:) = fecundity;
    
    Demographics.LM.(SiteNames{i}) = LM;
    
    %Check Matrix Subpopulation Stability
    [~,Diag] = eig(LM); % Calculate eigenvectors and eigenvalues of Leslie Matrix
    Lambda = max(Diag(Diag>0)); %Strictly dominant eigenvalue, i.e. growth rate
    props(:,1) = ((LM^100)*(10*ones(n,1)))/sum((LM^100)*(10*ones(n,1))); % Stable Life Stage Proportions
    
    %Generate Initial Population (N=K) with Stable Life Stage Proportions
    K = Demographics.CarryingCapacity.(SiteNames{i});
    pop_i(:,1) = props*K;

    Metapop.(SiteNames{i}).LM = LM;
    Metapop.(SiteNames{i}).Stages = pop_i;
    Metapop.(SiteNames{i}).Abundance = sum(pop_i);
    Metapop.(SiteNames{i}).Proportion = props;
    Metapop.(SiteNames{i}).Lambda = Lambda;
    Metapop.(SiteNames{i}).K = K;
    Metapop.(SiteNames{i}).Density = sum(pop_i)/K;

end


%% Create Matrix of Temporal Variance of Respective Rates

for i = 1:h
    for j = 1:n
        variance(j) = Demographics.Variance.(SiteNames{i}).(LifeStages{j});
        fecvar(j) = Demographics.FecundityVariance.(SiteNames{i}).(LifeStages{j});
        if j < n
        transition(j) = Demographics.Transition.(SiteNames{i}).(LifeStages{j}); 
        end
    end
    
    varstay = [variance(1:end-1).*(1-transition) variance(end)];
    vartrans = variance(1:end-1).*(transition);
    
    VarMat = diag(varstay);
    index = diag(vartrans,-1)~=0;
    VarMat(index) = vartrans;
    
    VarMat(1,:) = fecvar;
    
    Metapop.(SiteNames{i}).Variance =  VarMat;
    
end

%% Record Density-Dependent Offspring Survival Rates
for i = 1:h
    Metapop.(SiteNames{i}).EggSurvival = Demographics.EggSurvival.(SiteNames{i});
end

cd(ParentFolder)

end

