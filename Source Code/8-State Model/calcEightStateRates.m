function k = calcEightStateRates(pKs, AlternatingK)

% Accepts pKa and pKd values in a 1x4 vector
% Accepts AlternatingK, a 8x1 vector of AA values

pKs = 10.^(-1.*pKs);    % Make it a K value

protonOn = 1e10;        % per M per second
drugOn = 1e7;           % per M per second

Ka1 = pKs(1);          % for apo to singly-bound
Ka2 = pKs(2);          % for drug bound to simultaenously-bound

Kd1 = pKs(3);          % for apo to drug bound
Kd2 = pKs(4);          % for proton-bound to simultaenously-bound



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
          
k = [ProtonBindingK; DrugBindingK; AlternatingK];

end