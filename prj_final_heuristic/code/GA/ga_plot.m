function [ best_path,val ] = ga_plot( race )
%     load city_location.mat
    load x50.mat
    city_distance = preprocess(cityXY);
    cityXY = cityXY';
    [m,n]=size(race);
    point=zeros(n,2);
    adaptation=ga_adaptation(race);
    [val,index]=min(adaptation);
    best_path=race(index,:);
    for i=1:n
        point(i,:)=cityXY(best_path(i),:);
    end
    plot(point(:,1),point(:,2),'-o');
end

