function [ LM_adj ] = MatrixAdjustment( t,SuitabilityValues,Demographics,Inputs,SiteNames,LifeStages )
% Adjust Matrix Demographic Rates based on User Prescribed Habitat 
% Parameters (Inputs.surv_adj, Inputs.adv_adj, Inputs.fec_adj) and build
% new adjusted matrices for each subpopulation 

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

h = numel(fieldnames(SuitabilityValues)); %number of sites
n = numel(LifeStages); %number of life-stages

for i = 1:h
    
    survadjust = [];
    transadjust = [];
    fecadjust = [];
    
    %Pull Baseline Values
    for j = 1:numel(LifeStages)
        survival(j) = Demographics.Survival.(SiteNames{i}).(LifeStages{j}); 
        if j < n
        transition(j) = Demographics.Transition.(SiteNames{i}).(LifeStages{j}); 
        end
        fecundity(j) = Demographics.Reproduction.(SiteNames{i}).(LifeStages{j}); 
    end
    
    
    %% Survival Adjustment
    
    if isempty(Inputs.surv_adj)==0
        
    for j = 1:numel(LifeStages)
           
        for k = 1:numel(Inputs.surv_adj)
            Surv_Opts(k) = SuitabilityValues.(SiteNames{i}).(Inputs.surv_adj{k}).(LifeStages{j})(t);
        end
            
        % Adjust Based on Selected Fuzzy Aggregation Method
        if strcmp(Inputs.fuzzytype,'product')==1
            adjust = prod(Surv_Opts);
        elseif strcmp(Inputs.fuzzytype,'minimum')==1
            adjust = min(Surv_Opts);
        elseif strcmp(Inputs.fuzzytype,'geomean')==1
            adjust = prod(Surv_Opts)^(1/numel(Surv_Opts));
        else
            error('Error: Must Define Fuzzy Logic Method')
        end
        
        survadjust(j) = adjust;

    end
    
    else 
        
        survadjust = ones(1,numel(LifeStages));
        
    end
    
    %% Transition Adjustment (Does Not Include Last Stage)
    
    if isempty(Inputs.adv_adj)==0
    
    for j = 1:numel(LifeStages)-1
           
        for k = 1:numel(Inputs.adv_adj)
            Adv_Opts(k) = SuitabilityValues.(SiteNames{i}).(Inputs.adv_adj{k}).(LifeStages{j})(t);
        end
        
        % Adjust Based on Selected Fuzzy Aggregation Method
        if strcmp(Inputs.fuzzytype,'product')==1
            adjust = prod(Adv_Opts);
        elseif strcmp(Inputs.fuzzytype,'minimum')==1
            adjust = min(Adv_Opts);
        elseif strcmp(Inputs.fuzzytype,'geomean')==1
            adjust = prod(Adv_Opts)^(1/numel(Adv_Opts));
        else
            error('Error: Must Define Fuzzy Logic Method')
        end
        
        transadjust(j) = adjust;

    end
    
    else 
        
        transadjust = ones(1,numel(LifeStages)-1);
        
    end
    
    
    %% Fecundity Adjustment
                  
    if isempty(Inputs.fec_adj)==0
    
    for j = 1:numel(LifeStages)
           
        for k = 1:numel(Inputs.fec_adj)
            Fec_Opts(k) = SuitabilityValues.(SiteNames{i}).(Inputs.fec_adj{k}).(LifeStages{j})(t);
        end
        
        % Adjust Based on Selected Fuzzy Aggregation Method
        if strcmp(Inputs.fuzzytype,'product')==1
            adjust = prod(Fec_Opts);
        elseif strcmp(Inputs.fuzzytype,'minimum')==1
            adjust = min(Fec_Opts);
        elseif strcmp(Inputs.fuzzytype,'geomean')==1
            adjust = prod(Fec_Opts)^(1/numel(Fec_Opts));
        else
            error('Error: Must Define Fuzzy Logic Method')
        end
        
        fecadjust(j) = adjust;

    end
       
    else 
        
        fecadjust = ones(1,numel(LifeStages));
       
    end
    
    %Construct Habitat-Adjusted Matrix
    surv_t = survival.*survadjust;
    trans_t = transition.*transadjust;
    fec_t = fecundity.*fecadjust;
    
    survstay = [surv_t(1:end-1).*(1-trans_t) surv_t(end)];
    survtrans = surv_t(1:end-1).*(trans_t);
    
    LM = diag(survstay);
    index = diag(survtrans,-1)~=0;
    LM(index) = survtrans;
    
    LM(1,:) = fec_t;
    LM_adj.(SiteNames{i}) = LM; %Habitat Adjusted Lefkovitch Matrix
    
end
    
end