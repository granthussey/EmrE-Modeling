function [t,Sol,Tint,Text, MolecularSol, PercentSol, MassBalance] = runTenState(k, Times, pHs, Drugs, EmrE, Volumes)

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
Smem =  calcSmem(Vint);            % Program outputs units of dm^2
                                    % for an inputted internal volume. 

Etot = (EmrE)*Vtot*6.02e23/Smem;  % Etot in units of mc/dm^2


% Initialize drug parameters
Dtot = ((Drugs(1)*Vint) + (Drugs(2)*Vext))/Vtot;   % Dtot in units of mc
Dint = Drugs(1);                                

Hext = 10^(-pHs(2));                 % Hext in units of M
Hint = 10^(-pHs(1));                 % Hint in units of M

% Set initial conditions by variables above
iniCond = [Etot/6,Etot/6,Etot/6,0,0,0,0,Etot/6,Etot/6,Dint];

% VolFactor
% Changes units of mc/dm^2 to M.
VolFactInt  = Smem/(Vint*6.02e23);       % volFactorInt -> To M in terms of internal volume.
VolFactTot  = Smem/(Vtot*6.02e23);       % VolFactTot   -> To M in terms of total volume.

%% ODE System

% This is the newly minted dxdt matrix
dxdt = @(t,x) [k(3)*Hint*x(3) + k(16)*x(5) + k(18)*(Etot-x(1)-x(2)-x(3)-x(4)-x(5)-x(6)-x(7)-x(8)-x(9)) + k(30)*x(9) - k(4)*x(1) - k(15)*x(10)*x(1) - k(17)*x(1) - k(29)*x(1)*Hint;...
    k(2)*(Etot-x(1)-x(2)-x(3)-x(4)-x(5)-x(6)-x(7)-x(8)-x(9)) + k(10)*x(6) + k(19)*x(3) - k(1)*Hext*x(2) - k(9)*((Vtot/Vext)*(Dtot-(x(10)*Vint/Vtot)-(x(4)+x(5)+x(6)+x(7))*VolFactTot))*x(2) - k(20)*x(2);...
    k(4)*x(1) + k(12)*x(7) + k(20)*x(2) - k(3)*Hint*x(3) - k(11)*x(10)*x(3) - k(19)*x(3);...
    k(5)*Hext*x(6) + k(13)*((Vtot/Vext)*(Dtot-(x(10)*Vint/Vtot)-(x(4)+x(5)+x(6)+x(7))*VolFactTot))*(Etot-x(1)-x(2)-x(3)-x(4)-x(5)-x(6)-x(7)-x(8)-x(9)) + k(23)*x(5) - k(6)*x(4) - k(14)*x(4) - k(24)*x(4);...
    k(7)*Hint*x(7) + k(15)*x(10)*x(1) + k(24)*x(4) - k(8)*x(5) - k(16)*x(5) - k(23)*x(5);...
    k(6)*x(4) + k(9)*((Vtot/Vext)*(Dtot-(x(10)*Vint/Vtot)-(x(4)+x(5)+x(6)+x(7))*VolFactTot))*x(2) + k(21)*x(7) - k(5)*Hext*x(6) - k(10)*x(6) - k(22)*x(6);...
    k(8)*x(5) + k(11)*x(10)*x(3) + k(22)*x(6) - k(7)*Hint*x(7) - k(12)*x(7) - k(21)*x(7);...
    k(27)*Hext*(Etot-x(1)-x(2)-x(3)-x(4)-x(5)-x(6)-x(7)-x(8)-x(9))+k(25)*x(9)-k(28)*x(8)-k(26)*x(8);...
    k(29)*Hint*x(1)+k(26)*x(8)-k(30)*x(9)-k(25)*x(9);...
    VolFactInt*(k(12)*x(7)+k(16)*x(5)-k(11)*x(3)*x(10)-k(15)*x(1)*x(10))];


tic
[t,Sol] = ode15s(dxdt,[Times(1) Times(2)],iniCond);
toc

Tint = Sol(:,10);
Text = (Dtot - ((Sol(:,10).*Vint./Vtot)) - (sum(Sol(:,4:7),2).*Smem./Vtot./6.02e23)).*Vtot./Vext;

% An easy-to-use matrix for getting molecular solution data
MolecularSol = zeros(length(t),12);                     % Initialize the matrix
MolecularSol(:,1) = (Etot - sum(Sol(:,1:9),2)).*Smem;   % Back-calculate the mass-balanced EHext species
MolecularSol(:,2:10) = (Sol(:,1:9)).*Smem;
MolecularSol(:,11) =  Text*Vext*6.02e23;
MolecularSol(:,12) = Tint*Vint*6.02e23;

% An easy-to-use matrix for comparing mass balance.
MassBalance = sum(MolecularSol(:,:),2);  % This should stay constant.

PercentSol = MolecularSol;
PercentSol(:,1:10) = PercentSol(:,1:10)./(Etot*Smem);
PercentSol(:,11:12) = PercentSol(:,11:12)./(Dtot*Vtot*6.02e23);
PercentSol = PercentSol.*100;


end


