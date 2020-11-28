function solution_new = ga_hamilton(path)
load x50.mat
city_distance = preprocess(cityXY);
len=length(path);
solution_new=path;
for i=3:len-3
    for j=i+1:len-2
        d1=city_distance(solution_new(i-1),solution_new(i))+city_distance(solution_new(i),solution_new(i+1));
        d2=city_distance(solution_new(j-1),solution_new(j))+city_distance(solution_new(j),solution_new(j+1));
        d3=city_distance(solution_new(i-1),solution_new(j))+city_distance(solution_new(j),solution_new(i+1));
        d4=city_distance(solution_new(j-1),solution_new(i))+city_distance(solution_new(i),solution_new(j+1));
        if((d1+d2)>(d3+d4))
            temp=solution_new(i);
            solution_new(i)=solution_new(j);
            solution_new(j)=temp;
        end
    end
end
end

