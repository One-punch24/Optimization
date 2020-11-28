% % 20 cities
clear,clc
close all
% load city_distance.mat
% load cityXY.mat
load x50.mat
city_distance = preprocess(cityXY);
cityXY = cityXY';

City_Number=50;         %city numbers
Race_Number=50;        %population
Iteration=10000;          %iterations
P_Cross=0.3;            %cross prob.
P_Mutation=0.7;         %variatino prob.
race=zeros(Race_Number,City_Number+2);
%initialize the race
for i=1:Race_Number                         
    temp=randperm(City_Number);
    route=[City_Number+1,temp,City_Number+1];
    route=ga_hamilton(route);        %optimal loop path algorithm--->initial value optimization
    race(i,:)=route;
end
for t=1:Iteration
    adaptation=ga_adaptation(race);         %calculate the adaptation function
    race=ga_choose(race,adaptation);        %selection function
    race=ga_cross(race,P_Cross);            %cross
    race=ga_mutation(race,P_Mutation);      %variation
    [path,val]=ga_plot(race);
    pause(0.1);
    fprintf('%dth generation, optimal: %d\n',t,val);
end
