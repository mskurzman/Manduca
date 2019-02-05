function [legs,muscles] = manducaGenerateInitialSolution()
legs = zeros(10,5);
muscles = zeros(10,4);
%legs matrix
%min 2 legs locked. no 5 legs locked
for i=1:10
num_legs = randi([2,4]);
locations = randsample(5,num_legs);
legs(i,locations)=1;
end
%muscles matrix
%max 3 muscles contracting. no contraction if adjacent legs are locked
for j=1:10
    count=0;
    for k=1:4
    if (legs(j,k)~=1 || legs(j,k+1)~=1) && count<3
        if rand()<0.5
            muscles(j,k)=100;
            count=count+1;
        end
    end
    end
end

