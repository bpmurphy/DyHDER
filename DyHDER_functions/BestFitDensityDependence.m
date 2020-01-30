function [ Metapop ] = BestFitDensityDependence( Inputs,Metapop,SiteNames )
% Compute the best-fit beta parameter from density-dependent offspring 
% survival rates based on the selected density-dependent model

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

%Calculated for each site
for i = 1:numel(SiteNames) 
    
    % Survival to age-1 at carrying capacity
    s1 = Metapop.(SiteNames{i}).EggSurvival(1) ;
    
    % Maximum possible survival to age-1 (as site n/K approaches 0)
    s0 = Metapop.(SiteNames{i}).EggSurvival(2) ;
    
    % Ricker Model
    if Inputs.DD_model == 1 
        BetaFec = log(s0/s1); % %Fitting Parameter

    % Beverton-Holt Model
    elseif Inputs.DD_model == 2 
        BetaFec = (s0/s1)-1; %Fitting Parameter

    end

    % Find maternity rate from Leslie Matrix and Prescribed Survival Rate
    Fecundity = Metapop.(SiteNames{i}).LM(1,:);
    Metapop.(SiteNames{i}).OffSpring = Fecundity/s1;
    Metapop.(SiteNames{i}).Beta = BetaFec;
    
end

end

