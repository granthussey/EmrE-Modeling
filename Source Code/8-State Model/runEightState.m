function [t,Sol,Tint,Text, MolecularSol, PercentSol, MassBalance] = runEightState(k, Times, pHs, Drugs, EmrE, Volumes)

%% Notes

% pH is [pHint pHext]
% Drugs is [intDint iniDext]
% Volumes is [Vint Vext]

% Since mass action systems has a globally attracting fixed point, 
% initial conditions don't matter.

% I seed each non-drug state with an equal fraction of EmrE to minimize
% numerical error.

%% Set Initial Conditions

Vint = Volumes(1);
Vext = Volumes(2);
Vtot = Volumes(1) + Volumes(2);
Smem = calcSmem(Vint);            % Program outputs units of dm^2 
                                    % for an inputted internal volume. 
                                    
Etot = (EmrE)*Vtot*6.02e23/Smem;  % Etot in units of mc/dm^2


% Initialize drug parameters
Dtot = ((Drugs(1)*Vint) + (Drugs(2)*Vext))/Vtot;   % Dtot in units of M
Dint = Drugs(1);                                  

Hext = 10^(-pHs(2));                 % Hext in units of M
Hint = 10^(-pHs(1));                 % Hint in units of M

% Set initial conditions by variables above
iniCond = [Etot/4,Etot/4,Etot/4,0,0,0,0,Dint];

% VolFactor
VolFactInt  = Smem/(Vint*6.02e23);       % volFactorInt -> To M in terms of internal volume.
VolFactTot  = Smem/(Vtot*6.02e23);       % VolFactTot   -> To M in terms of total volume.

%% ODE System

% This is the newly minted dxdt matrix
dxdt = @(t,x) [k(3)*Hint*x(3) + k(16)*x(5) + k(18)*(Etot-x(1)-x(2)-x(3)-x(4)-x(5)-x(6)-x(7)) - k(4)*x(1) - k(15)*x(8)*x(1) - k(17)*x(1);...
    k(2)*(Etot-x(1)-x(2)-x(3)-x(4)-x(5)-x(6)-x(7)) + k(10)*x(6) + k(19)*x(3) - k(1)*Hext*x(2) - k(9)*((Vtot/Vext)*(Dtot-(x(8)*Vint/Vtot)-(x(4)+x(5)+x(6)+x(7))*VolFactTot))*x(2) - k(20)*x(2);...
    k(4)*x(1) + k(12)*x(7) + k(20)*x(2) - k(3)*Hint*x(3) - k(11)*x(8)*x(3) - k(19)*x(3);...
    k(5)*Hext*x(6) + k(13)*((Vtot/Vext)*(Dtot-(x(8)*Vint/Vtot)-(x(4)+x(5)+x(6)+x(7))*VolFactTot))*(Etot-x(1)-x(2)-x(3)-x(4)-x(5)-x(6)-x(7)) + k(23)*x(5) - k(6)*x(4) - k(14)*x(4) - k(24)*x(4);...
    k(7)*Hint*x(7) + k(15)*x(8)*x(1) + k(24)*x(4) - k(8)*x(5) - k(16)*x(5) - k(23)*x(5);...
    k(6)*x(4) + k(9)*((Vtot/Vext)*(Dtot-(x(8)*Vint/Vtot)-(x(4)+x(5)+x(6)+x(7))*VolFactTot))*x(2) + k(21)*x(7) - k(5)*Hext*x(6) - k(10)*x(6) - k(22)*x(6);...
    k(8)*x(5) + k(11)*x(8)*x(3) + k(22)*x(6) - k(7)*Hint*x(7) - k(12)*x(7) - k(21)*x(7);...
    VolFactInt*(k(12)*x(7)+k(16)*x(5)-k(11)*x(3)*x(8)-k(15)*x(1)*x(8))];

tic
[t,Sol] = ode15s(dxdt,[Times(1) Times(2)],iniCond);
toc

Tint = Sol(:,8);
Text = (Dtot - ((Sol(:,8).*Vint./Vtot)) - (sum(Sol(:,4:7),2).*Smem./Vtot./6.02e23)).*Vtot./Vext; % Back-calculated via mass balance

% An easy-to-use matrix for getting molecular solution data
MolecularSol = zeros(length(t),10);                     % Initialize the matrix
MolecularSol(:,1) = (Etot - sum(Sol(:,1:7),2)).*Smem;   % Back-calculate the mass-balanced EHext species
MolecularSol(:,2:8) = (Sol(:,1:7)).*Smem;               
MolecularSol(:,9) =  Text*Vext*6.02e23;
MolecularSol(:,10) = Tint*Vint*6.02e23;

% An easy-to-use matrix for comparing mass balance.
MassBalance = sum(MolecularSol(:,:),2);  % This should stay constant.

PercentSol = MolecularSol;
PercentSol(:,1:8) = PercentSol(:,1:8)./(Etot*Smem);
PercentSol(:,9:10) = PercentSol(:,9:10)./(Dtot*Vtot*6.02e23);
PercentSol = PercentSol.*100;

end


