function [ SuitabilityValues ] = ConvertHabitatToSuitability( HabitatMetrics,SuitabilityRelations,SiteNames,LifeStages )
% Convert Time-Series of Habitat Metrics to Time-Series of Suitability 
% Values based on species-specific habitat suitability relations

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

%Check Input Names for Conversion
a=fieldnames(HabitatMetrics.(SiteNames{1}));
b=fieldnames(SuitabilityRelations);
if isempty(setdiff(a,b))==0
    error('Suitability and Habitat Metric Names Do Not Match')
end

tmax = numel(HabitatMetrics.(SiteNames{1}).(a{1}));

%Site Iteration
for i = 1:numel(SiteNames)
    
    %Metric Iteration
    for j = 1:numel(a)
        
        check_fields = fieldnames(SuitabilityRelations.(a{j}));
        
        %If One Habitat Suitability Relation for All Life-Stages
        if strcmp(check_fields(2,:),'All')==1
            
            Opt = [];
            for t = 1:tmax

                %Find Value if Exact Match
                habmet = HabitatMetrics.(SiteNames{i}).(a{j})(t);
                exact = find(SuitabilityRelations.(a{j}).Value==habmet);

                if isempty(exact) == 0
                    Opt(t,1) = SuitabilityRelations.(a{j}).All(exact);
                
                %Linear Interpolation if not Exact Match
                elseif isempty(exact) == 1
                    up = find(SuitabilityRelations.(a{j}).Value>habmet,1,'first');
                    down = find(SuitabilityRelations.(a{j}).Value<habmet,1,'last');
                    x = [SuitabilityRelations.(a{j}).Value(down),SuitabilityRelations.(a{j}).Value(up)];
                    y = [SuitabilityRelations.(a{j}).All(down),SuitabilityRelations.(a{j}).All(up)];
                    exact = interp1(x,y,habmet);
                    Opt(t,1) = exact;
                end

            end

            for kk = 1:numel(LifeStages)

                SuitabilityValues.(SiteNames{i}).(a{j}).(LifeStages{kk})=Opt;

            end
        
        %If Separate Habitat Suitability Relations for Each Life-Stage
        else
        
            for kk = 1:numel(LifeStages)
            
            Opt = [];
            for t = 1:tmax

                %Find Value if Exact Match
                habmet = HabitatMetrics.(SiteNames{i}).(a{j})(t);
                exact = find(SuitabilityRelations.(a{j}).Value==habmet);
    
                %Linear Interpolation if not Exact Match
                if isempty(exact) == 1
                    up = find(SuitabilityRelations.(a{j}).Value>habmet,1,'first');
                    down = find(SuitabilityRelations.(a{j}).Value<habmet,1,'last');
                    x = [SuitabilityRelations.(a{j}).Value(down),SuitabilityRelations.(a{j}).Value(up)];
                    y = [SuitabilityRelations.(a{j}).(LifeStages{kk})(down),SuitabilityRelations.(a{j}).(LifeStages{kk})(up)];
                    exact = interp1(x,y,habmet);
                end

                Opt(t,1) = exact; %Set Suitability Value

            end
            
            SuitabilityValues.(SiteNames{i}).(a{j}).(LifeStages{kk})=Opt;
            
            end
            
        end
        
    end

end

end

