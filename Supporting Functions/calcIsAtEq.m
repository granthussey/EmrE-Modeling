function EqRatio = calcIsAtEq(Matrix)

% Returns a value to aid in determining if a matrix is at equilibrium. 
% Divides the last fifth of the run by the previous last fifth of the run.

EqRatio = mean(Matrix(round(0.95*end:end,0)))/mean(Matrix(round(0.9*end,0):(round(0.95*end,0))));

end