function EqT = calcEqT(Matrix)

% Calculates the equilibrium concentration based on the average of the last
% 10th of the run. 

EqT = mean(Matrix(round(0.9*end,0):end)); 

end