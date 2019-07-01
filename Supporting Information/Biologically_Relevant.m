

% This script analyzes the equilibration of drug across the membrane using
% physically-relevant parameters.

% Conditions
Times = [0, 100];   % in seconds
pHs = [7.5,6.5];

% Physical parameters
EmrE = 20e-9;
Vint = 2.64e-7;
Vtot = 0.001;
Vext = Vtot - Vint;
Volumes = [Vint Vext];
Drugs = [50e-9 0];

% Calculate rate constants
EmrE_pK = [8.2, 6.8, 7, 7.2, 6];
EmrE_AltK = [100;...           % k17 EH
    100;...                    % k18 EH
    40;...                     % k19 E
    40;...                     % k20 E
    7.3;...                    % k21 ED
    7.3;...                    % k22 ED
    8.9;...                    % k23 EHD
    8.9;...                    % k24 EHD
    220;...                    % k25 EHH
    220];                      % k26 EHH
k = calcTenStatesRates(EmrE_pK, EmrE_AltK);

% Run simulation
[t,Sol,Dint,Dext, MolSol, PerSol, MassBal] = runTenState(k, Times, pHs, Drugs, EmrE, Volumes);


