
function [neigh_legs,neigh_muscles] = permuteManduca(legs,muscles,curr_temp)
neigh_legs=legs;
neigh_muscles=muscles;
curr_temp=10;
m=10*exp(-1/curr_temp);
n=5*exp(-1/curr_temp);
%m=10*curr_temp/initial_temp;
%n=5*curr_temp/initial_temp;
rows=randsample(10,int16(m));
for i=1:n+1
    loc=randi(5);
    if neigh_legs(rows,loc)==0
        neigh_legs(rows,loc)=1;
    else
        neigh_legs(rows,loc)=0;
    end
end
%check for constraints and change matrix as needed
for j=1:10
    summ=sum(neigh_legs(j,:));
    if summ==5
        neigh_legs(j,randi(5))=0;
    elseif summ==0
        neigh_legs(j,randsample(5,2))=1;
    else
        while summ==1
            neigh_legs(j,randi(5))=1;
            summ=sum(neigh_legs(j,:));
        end
    end
end

%muscles
k=4*exp(-1/curr_temp);
%uncomment this to select new random rows for muscles
%rowsM=randsample(10,int16(m));
for i=1:k+1
    loc=randi(4);
    if neigh_muscles(rows,loc)==0
        neigh_muscles(rows,loc)=100;
    else
        neigh_muscles(rows,loc)=0;
    end
end

for j=1:10
    for k=1:4
        if neigh_legs(j,k)==1 && neigh_legs(j,k+1)==1 && neigh_muscles(j,k)==100
        neigh_muscles(j,k)=0;
        end
    end
    summ=sum(neigh_muscles(j,:));
    if summ==400
        neigh_legs(j,randi(4))=0;
    end
end




end

















