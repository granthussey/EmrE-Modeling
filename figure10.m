
% pH Conditions
pHs = [6.5,7.5];

% Run simulations for first figure 10.A then figure 10.B
[x1,Z1,R1,T1] = figure10scan(pHs,8.2,7);
[x2,Z2,R2,T2] = figure10scan(pHs,7,5);

% Save the run data within a cell array
WT = {x1,Z1};
Downshift = {x2,Z2};

% Pass Fig 10.A data to plotting program
ContourEmrE(WT);
title('pKa_{1} = 8.2, pKa_{2} = 7')
ylim([0.1,10]);

% Pass Fig 10.B data to plotting program
ContourEmrE(Downshift);
title('pKa_{1} = 7, pKa_{2} = 5')
ylim([0.1,10]);

