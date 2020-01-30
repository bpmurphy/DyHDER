function [ M ] = CreateMovementMatrix( j, t, DispersalMetrics, HSI, Record, SimName, SiteNames, LifeStages )
% Step-wise functions to create life-stage dependent transition matrix.
% For more detail on model see Murphy et al., (20xx)

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

%% Inputs:

%   C = h x h sized matrix of site to site connectivity
%       with rows, nu, representing dispersal origin
%       and columns, zeta, representing dispersal destination
%       and elements, 0 <= c(nu,zeta) <= 1 representing relative ease of
%       dispersal between locations nu and zeta independent of life-stage
%   D = h x h size matrix of site-to-site distance, independent of site
%       connectivity
%   D_s = Exponential scalar for distance dependent sedentary behavior
%   alfa = vector of individual site attractivity, length h
%       with elements, 0 <= alf(zeta) <= 1 representing attractivity of
%       locations species can disperse to; alf can be life-stage dependent
%       and can be factor of habitat suitabilility and/or population
%       density of destination site
%   q = propensity for dispersal of given species life-stage independent of
%       site attractivity, distance, connectivity, etc.


%%

[a,b]=size(DispersalMetrics.Connectivity);

%Rewrite variables to simpler names
C = DispersalMetrics.Connectivity;
D = DispersalMetrics.Distance;
D_s = DispersalMetrics.DistanceScalar(j);
q = DispersalMetrics.Dispersal(j);

%Extract Site Densities
for i = 1:b
    density(i) = Record.(SimName).(SiteNames{i}).Density(t);
    alfa(i) = HSI.(SiteNames{i}).(LifeStages{j})(t);
end

S = exp(-D./D_s); % Distance-dependent Sedentary Propensity

for nu = 1:b
    
    %Compute Emigration Probability
    P_emigration = q/((1/density(nu))*alfa(nu)^2*(1-q)+q);
    
    % Create temporary variable summing Immigration Probabilities for
    % Denominator in Site Immigration Probability
    temp = S(nu,:).*C(nu,:).*(alfa.^2).*(1./density);
    temp(:,nu) = [];
    summation = sum(temp);
    
    %If no probability of dispersal to other sites 
    if summation == 0
        
        for zeta = 1:a
            if zeta == nu
                M(nu,zeta)=1;
            else
                M(nu,zeta)=0;
            end
        end
        
    %If probability of dispersal to other sites
    else
    
        for zeta = 1:a
            if zeta == nu
                M(nu,zeta)=1-P_emigration; %Probability of staying put
            
            %Probability of immigrating to other sites
            else
                M(nu,zeta)=P_emigration*...
                    ((C(nu,zeta)*(alfa(zeta)^2)*S(nu,zeta)*...
                    (1/density(zeta))))./(summation);
            end
        end
        
    end
    
end


end

