%% DyHDER: Dynamic Habitat Disturbance & Ecological Resilience

% Description: This script runs the DyHDER matrix population model 
% and all associated functions herein. DyHDER was designed to  
% incorporate subpopulation habitat condition and connectivity 
% into a population viability analysis (PVA) framework.

% Written by: Dr. Brendan P. Murphy (Utah State University)
% Version: 1.0
% Published/Released: Jan. 30, 2020

% For more information on DyHDER, or if you use DyHDER in your research, 
% cite the following open-source, peer-reviewed publication:

% Murphy, B.P., T.E. Walsworth, P. Belmont, M.M. Conner, and P. Budy.
% 2020. Dynamic Habitat Disturbance and Ecological Resilience (DyHDER): 
% modeling population responses to habitat condition. Ecosphere 11:1-26. 
% DOI: 10.1002/ecs2.3023.


%% Copyright & Licensing

% DyHDER: Dynamic Habitat Disturbance & Ecological Resilience
% Copyright (C) 2020 Brendan P. Murphy 

% Developer can be contacted at bpmurphy@aggiemail.usu.edu

% This program is free software; you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation; either version 3 of the License, or 
% any later version.

% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License 
% along with this program. If not, see <https://www.gnu.org/licenses/>.


%% Directions: 
% Before Running DyHDER_Master.m: 

% 1. Directory Set Up: 
% This DyHDER_Master.m script must be located in a folder along with a 
% subfolder that contains all associated functions: 'DyHDER_scripts' 
% Another subfolder of any name (assigned below) must also be included 
% containing all required data input files (Excel spreadsheets). 

% 2. Model Parameters:
% Before running this master code, users must prescribe all model 
% parameters and options below (Sections A - E). 

% For more detailed description of input variables and how to setup input
% files, read the document: DyHDER_README_v1.0.pdf


%% Initiliaze Workspace 
clear
close all


%% A. Prescribe Model Parameters

% Name of Subfolder containing all Input Spreadsheets
% Date Type: String, ex: 'foldername'
InputsFolder = 'Inputs_LoganRiverExample';

%Name of Subfolder Containing All DyHDER Functions
% Date Type: String
DyHDERFunctions = 'DyHDER_functions';

% Names of each animal life-stage (ordered from youngest to oldest)
% Date Type: Cell Array of Comma-Delimited Strings (no spaces)
LifeStages = {'Age0','Juvenile','SmallAdult','LargeAdult'};

% Number of model simulations
% Data Type: Integer
Inputs.SimNum = 100; 

% Model maximum run time, in years
% Data Type: Integer
Inputs.tmax = 100; 

% Model time step, in years
% Data Type: Double
% Note: DyHDER to-date has only been run using annual timesteps, 
% (Inputs.dt=1) and bugs could arise with differing timestep
Inputs.dt = 1; 

% Select Years to Compare Distributions of Metapopulation Abundance
% Data Type: Integer
Inputs.PreDisturbance = round(0.30*Inputs.tmax); %years
Inputs.PostDisturbance = round(0.90*Inputs.tmax); %years

% Select Number of Final Model Years to Compile Sub- & Metapopulation 
% Pop. Density Data for Outputs
% Data Type: Integer
Inputs.CompileYears = 25;


%% B. Quasi-Extinction Thresholds:
% Whichever of the two inputs below equals the larger number of individuals
% for a given population in a given year will be the threshold applied:

% Fraction of Carrying Capacity (for both sub- & metapopulation)
% Data Type: Double (must be < 1)
Inputs.Nx = 0.05; 

% Absolute number of individuals
% Data Type: Integer
Inputs.QEx = 10;


%% C. Turn On/Off Optional Model Components

% Inputs must be either 0 or 1, where
% 0 = Off
% 1 = On 

% 1. Terminate model simulations if metapopulation abundance falls below
% quasi-extinction threshold
Inputs.terminate = 1; 

% 2. Suspend subpopulation matrix projection if subpopulation abundance 
% falls below quasi-extinction threshold
Inputs.Halt = 1; 

% 3. Environmental Stochasticity
Inputs.EnvStochasticity = 1; 

% 4. Dispersal
Inputs.Disperse = 1; 

% 5. Stocking/Harvesting
Inputs.Stock = 0; 

% 6. Save Inputs and Outputs?
Inputs.SaveData = 1;


%% D. Select Density Dependence Model and Fuzzy Aggregegation Method 

% Select Density Dependence Model:
% Ricker = 1 
% Beverton-Holt = 2
Inputs.DD_model = 2; 

% Assign Fuzzy Aggregation Method 
% Date Type: String = 'product', 'minimum', or 'geomean'
Inputs.fuzzytype = 'geomean';


%% E. Assign Habitat Metrics for Vital-Rate Adjustment
% Date Type: Cell Array of Comma-Delimited Strings (no spaces)
% Note: If no habitat metric is applied, leave empty cell brackets: {}
% Note: Input string must exactly match that from Excel Spreadsheets

Inputs.surv_adj = {'DO'};
Inputs.adv_adj = {'Temperature'};
Inputs.fec_adj = {};
Inputs.hsi_params = {'DO','Temperature'};




%% Inputs Complete %%






%% F. Run DyHDER model


%% Check Parameter Conflicts & Add DyHDER Functions to Path

% Parameter Check/Override
if Inputs.EnvStochasticity == 0
    Inputs.SimNum = 1; % Only 1 Simulation if Deterministic
end

% Add DyHDER Functions to Path
addpath(strcat(cd,'\',DyHDERFunctions))


%% Read-in Filenames of Input Files

[ Filenames,ParentFolder ] = FileReadIn( InputsFolder );

%% Read-in Demographic Inputs

[ Metapop,Demographics,LifeStages,SiteNames ] = ExtractDemographics( InputsFolder, Filenames, LifeStages );

%% Fit Density-Dependent Parameters based on Selected Model

[ Metapop ] = BestFitDensityDependence( Inputs,Metapop,SiteNames );

%% Read-in Dispersal Parameters

[ DispersalMetrics ] = ExtractDispersal( InputsFolder,Filenames,SiteNames );

%% Read-in Habitat Metric Time-series 

[ HabitatMetrics ] = ExtractHabitat( InputsFolder,Filenames,SiteNames );

%% Read-in Habitat Suitability Relations 

[ SuitabilityRelations ] = ExtractSuitability( InputsFolder,Filenames,LifeStages );

%% Convert Habitat Metrics to Suitability Values

[ SuitabilityValues ] = ConvertHabitatToSuitability( HabitatMetrics,SuitabilityRelations,SiteNames,LifeStages );

%% (Optional) Stocking/Harvesting 

if Inputs.Stock == 1
    [ Stocking ] = ExtractStock( InputsFolder,Filenames );
end

%% Compute Time-Series of Site Habitat Suitability Indices 

HSI = HabitatSuitabilityIndex(SuitabilityValues,Metapop,Inputs,SiteNames,LifeStages ); 

%% Miscellaneous Parameters/Set-up Procedures

h = numel(SiteNames); % Number of Subpopulation sites
n = numel(LifeStages); % Number of Life-Stages

% All Input Parameters Now Read-in from Input Folder: 
% Set Directory Back to Parent Folder
cd(ParentFolder) 


%% For-loop of Simulations (if multiple stochastic simulations)

for ss = 1:Inputs.SimNum
    
%% Reset Initial Conditions for Each Simulation

% Stochastic Behavior: Randomize MATLAB Random Generator 
rand('state',sum(1e5.*clock));

% Set-up Simulation Records
SimName = sprintf('Simulation%d',ss);
Record.(SimName) = Metapop;

% Calculate Initial Metapopulation Abundance
N(1) = 0; %Default Value 
K_total(1) = 0; %Default Value 
for i = 1:h
    N(1) = N(1) + Record.(SimName).(SiteNames{i}).Abundance(1); 
    K_total(1) = K_total(1) + Record.(SimName).(SiteNames{i}).K(1);
end

%% For-loop for Temporal Iteration for t = Timestep 2 : Max Model Years

for t = 2:Inputs.dt:Inputs.tmax

%% Carrying Capacity Assignment per Subpopulation
    
    % For number of Subpopulation sites
    for i = 1:h
        Record.(SimName).(SiteNames{i}).K(t) = Record.(SimName).(SiteNames{i}).K(t-1);
    end
    
%% Habitat-Based Vital Rate Adjustment: Survival, Transition & Reproduction

    LM_adj = MatrixAdjustment( t,SuitabilityValues,Demographics,Inputs,SiteNames,LifeStages );

%% Density-Dependent (Reproduction) Adjustment

    [ LM_adj,Record ] = DensityDependence( t,Inputs,Record,SimName,Metapop,SiteNames,LM_adj );

%% (Optional) Environmental Stochasticity 

    if Inputs.EnvStochasticity == 1
        LM_adj = EnvironmentalStochasticity( LM_adj,Metapop,SiteNames );  
    end
        
%% Truncate Normal Distributions: 

% If Survival or Transition Probabilities > 1, then = 1
    for i = 1:h
        [row,col]=ind2sub([numel(LifeStages),numel(LifeStages)],find(LM_adj.(SiteNames{i})>1));
        index1 = find(row>1);
        if isempty(index1)==0
            for j = 1:index1
                LM_adj.(SiteNames{i})(row(index1(j)),col(index1(j)))=1;
            end
        end
    end
    
% If any vital rate < 0, then = 0
    for i = 1:h
        LM_adj.(SiteNames{i})(LM_adj.(SiteNames{i})<0)=0;
        Record.(SimName).(SiteNames{i}).AdjustedMatrix{t} = LM_adj.(SiteNames{i});
    end
        
    
%% Matrix Projection
    
% (Optional) Suspend Matrix Projection If < max(Inputs.QEx,Frac_Thresh)
    if Inputs.Halt == 1
        
        for i = 1:h
            Frac_Thresh = Inputs.Nx * Record.(SimName).(SiteNames{i}).K(t);
            Halt_Threshold = max(Inputs.QEx,Frac_Thresh); 
        
            %If Above Threshold: Project Population Using Adjusted Matrices
            if Record.(SimName).(SiteNames{i}).Abundance(t-1) > Halt_Threshold
                
                Record = MatrixProjection( t,Record,SimName,SiteNames,LM_adj );
  
            %If Below Threshold: Do Not Project Population
            else
                
                Record.(SimName).(SiteNames{i}).Stages(:,t) = ...
                    Record.(SimName).(SiteNames{i}).Stages(:,t-1);

                Record.(SimName).(SiteNames{i}).Abundance(t) = ...
                    Record.(SimName).(SiteNames{i}).Abundance(t-1);

                Record.(SimName).(SiteNames{i}).Proportion(:,t) = ...
                    Record.(SimName).(SiteNames{i}).Proportion(:,t-1);

                Record.(SimName).(SiteNames{i}).Density(t) = ...
                    Record.(SimName).(SiteNames{i}).Density(t-1);
                
            end
        end
    
    %Project Populations Regardless of Quasi-extinction Threshold
    elseif Inputs.Halt == 0
        
        for i = 1:h
            Record = MatrixProjection( t,Record,SimName,SiteNames,LM_adj );
        end
        
    end

    
%% Compute Strictly Dominant Eigenvalue of Fully Adjusted Matrices

    Record = Eigenvalue( t,Record,SimName,LM_adj,SiteNames );
    
%% (Optional) Metapopulation Dispersal
    
    if Inputs.Disperse == 1 
        
        Record = MetapopDispersalModel( t, DispersalMetrics, HSI, Record, SimName, SiteNames, LifeStages );
        
    end

%% (Optional) Stocking/Harvesting

    if Inputs.Stock == 1

        Record = StockSites( t,SimName,Record,Stocking,LifeStages );

    end

%% Compute Final Metapopulation Abundance & Compare to Quasi-Extinction

    % Write new cell values (default = 0)
    N(t) = 0; %Metapopuluation Abundance in time t
    K_total(t) = 0; %Metapopuluation Carrying Capacity in time t
    
    for i = 1:h
        N(t) = N(t) + Record.(SimName).(SiteNames{i}).Abundance(t);
        K_total(t) = K_total(t) + Record.(SimName).(SiteNames{i}).K(t);
    end
    
%% (Optional) Check for Quasi-extinction (Terminate Simulation)

    if Inputs.terminate == 1
        if N(t) < Inputs.Nx * K_total(t)
           
            %Write out NaN to max time
            N(t:Inputs.dt:Inputs.tmax) = NaN; 
            K_total(t:Inputs.dt:Inputs.tmax) = NaN;
            
            for i = 1:h
                Record.(SimName).(SiteNames{i}).Abundance(t:Inputs.dt:Inputs.tmax) = NaN;
                Record.(SimName).(SiteNames{i}).Stages(:,t:Inputs.dt:Inputs.tmax) = NaN;
            end
            
            break %If Below Quasi-Extinction, Then End Simulation 
            
        end
    end
          
end


%% Record Final Conditions for Each Simulation

Population.(SimName).Abundance = N;
Population.(SimName).FinalAbundance = N(end);
Population.(SimName).TotalYears = t;
Population.(SimName).CarryingCapacity = K_total;

%Compute Annualized Growth Rate for Metapopulation
[ Record ] = AnnualizedGrowthRate( Record,SimName,SiteNames );


end


%% Calculate Percent of Extinctions From All Run Simulations

for i = 1:Inputs.SimNum
    SimName = sprintf('Simulation%d',i);
    Compile(i) = Population.(SimName).FinalAbundance;
end

prob_ext = sum(Compile<=(Inputs.Nx*K_total(1)))/Inputs.SimNum; 

%% G. PLOTTING
% Any or all of these plot functions can be commented out if not desired

% Metapopulation Abundance (If Deterministic, Script Plots All Subpops)
Plot_MetapopAbundance_AllSims( Inputs,SimName,SiteNames,Population,Record,K_total,prob_ext )

% Density Dependence Model 
Plot_DensityDependence( Inputs,Metapop,LifeStages )

% Histograms of Pre & Post Disturbance Population Abundance Distributions
binsize = K_total(1)/30; %set bin size (default: 1/30th of initial K)
Plot_PreAndPostHistograms( Inputs,Population,K_total,binsize )

% CDF of All Years & All Simulation Metapopulation Density
if Inputs.SimNum > 1
    Plot_CDF_AllSimulations( Inputs,Population,K_total )
end

% Plot Density Time-Series For Metapopulation and Subpopulations
% Currently layout is only properly configured for seven sites
if h == 7
    Plot_AllDensityPlots( Inputs,Record,SimName,SiteNames )
end

% Compile Normalized Metapopulation and Subpopulation Results
ModelResults = CompileResults( Inputs,Record,SimName,SiteNames);

%% H. Save Inputs & Outputs
% Written to a date and time-stamped output folder within the parent folder

if Inputs.SaveData == 1
t = char(datetime('now','TimeZone','local','Format','MMM_dd_y_HH_mm'));
OutputFolder = strcat('Outputs_',t);
mkdir(OutputFolder)
cd(OutputFolder)
save('Inputs.mat','Inputs','Filenames','HSI',...
    'DispersalMetrics','HabitatMetrics','Metapop','Demographics',...
    'SuitabilityRelations','SuitabilityValues')
save('Outputs.mat','Population','Record','prob_ext')
writetable(ModelResults,'Results.xls')
cd(ParentFolder)
end


%% I. Clear Variables Except for Key Input & Output Structures & Tables

clearvars -except Inputs Filenames HSI DispersalMetrics ...
    Demographics HabitatMetrics Metapops Population Record prob_ext ...
    SuitabilityRelations SuitabilityValues ModelResults 
