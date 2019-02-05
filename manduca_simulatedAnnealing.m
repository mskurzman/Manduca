function [cur_legs,cur_muscles,cur_fit] = manduca_simulatedAnnealing(initial_temp)
global report
[cur_legs,cur_muscles]=manducaGenerateInitialSolution();
cur_fit=manducaFitness(cur_legs);
done=false;
currentIt=0;
curr_temp=initial_temp;
while (~done)
    [neigh_legs,neigh_muscles]=permuteManduca(cur_legs,cur_muscles,curr_temp);
    neighScore=manducaFitness(neigh_legs,neigh_muscles,report);       
    if accept(cur_fit,neighScore,curr_temp)         
        cur_fit=neighScore;
        cur_legs=neigh_legs;
        cur_muscles=neigh_muscles;
    end
    curr_temp=manduca_updateTemp(curr_temp);
    currentIt=currentIt+1;
    done=currentIt==5000;
end
end


function [accept] = accept(curFit,neighborFit,T)
deltaS=curFit-neighborFit;
if(deltaS<0)
    accept=true;
else
    prob=exp(-deltaS/T);
    accept=rand()<prob;
end
end 

function [new_temp] = manduca_updateTemp(current_temp)
new_temp=0.999*current_temp;
end
