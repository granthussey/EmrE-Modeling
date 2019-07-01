function k = calcTenStatesRates(pKs, AlternatingK)

% Accepts pKa and pKd values in a 1x5 vector
% Accepts AlternatingK, a 10x1 vector of AA values

% This is a modified code from the 8-state. Similarly, the extra
% states/transitions/rates added on with the extra two states are in their
% own separate, annotated-module. 

pKs = 10.^(-1.*pKs);    % Make it a K value

protonOn = 1e10;        % per M per second
drugOn = 1e7;           % per M per second

Ka1 = pKs(1);          % for apo to singly-bound
Ka2 = pKs(2);          % for drug bound to simultaenously-bound
Ka3 = pKs(3);          % for singly-bound to doubly-bound

Kd1 = pKs(4);          % for apo to drug bound
Kd2 = pKs(5);          % for proton-bound to simultaenously-bound

ProtonBindingK = [protonOn;...          % k1             
                  protonOn*Ka1;...      % k2
                  protonOn;...          % k3  
                  protonOn*Ka1;...      % k4
                  protonOn;...          % k5
                  protonOn*Ka2;...      % k6
                  protonOn;...          % k7
                  protonOn*Ka2];        % k8
             
DrugBindingK = [drugOn;...              % k9
                drugOn*Kd1;...          % k10
                drugOn;...              % k11
                drugOn*Kd1;...          % k12
                drugOn;...              % k13
                drugOn*Kd2;...          % k14
                drugOn;...              % k15
                drugOn*Kd2];            % k16

DoubleProtK = [protonOn;...                  % k27
              protonOn*Ka3;...               % k28
              protonOn;...                   % k29
              protonOn*Ka3];...              % k30
           
         
k = [ProtonBindingK; DrugBindingK; AlternatingK; DoubleProtK];


end