function [x,Z,T1,T2] = figure10scan(pHs,pKa1,pKa3)

%% Setup

% This is an EmrE-like AA matrix but with adjustable alternating access
% rates for the ED and EHD alternating access steps.

AltK = @(ED,EHD) [100;...      % k17 PROTON
    100;...                    % k18 PROTON
    40;...                     % k19 APO
    40;...                     % k20 APO
    ED;...                     % k21 DRUG
    ED;...                     % k22 DRUG
    EHD;...                    % k23 DUO
    EHD;...                    % k24 DUO
    220;...                    % k25 DoublyProt
    220];                      % k26 DoublyProt

% This is an anonymous function to easilly calculate K from pK.
pKtoRate = @(pK) (10.^(-1.*pK));

%% Initialize substrate-off parameters.

% Initialize ranges for pKa and pKd values
pKd1_Range = (3:.125:9)';
pKd2_Range = flipud((3:.125:9)');
npKd = length(pKd1_Range);

A = [pKa1.*ones(npKd,1), pKa3.*ones(npKd,1), pKd1_Range, pKd2_Range]; 

% Solve for pKa2 based on thermodynamic cycle
pKa2 = A(:,1) + A(:,4) - A(:,3);

pK_Matrix = [A(:,1),pKa2,A(:,2),A(:,3),A(:,4)];

% Double-check thermo cycles
Check = pK_Matrix(:,1) + pK_Matrix(:,5) == pK_Matrix(:,4) + pK_Matrix(:,2);
ThrownOut = nnz(~Check);  % This is 0 as it is coded correctly. 

if ThrownOut > 0
    warning('Something is wrong with the current thermocycles.');
end

% Initialize a parallel list recording each Kd ratio.
pK_RatioList = pKtoRate(pK_Matrix(:,4))./pKtoRate(pK_Matrix(:,5));
npK_Ratios = length(pK_RatioList);

%% Initialize alternating access parameters.

kAAED = linspace(1,100,npKd);
kAAEHD = fliplr(kAAED);

AA_Matrix = [kAAED',kAAEHD'];
AA_RatioList = AA_Matrix(:,1)./AA_Matrix(:,2);

nAA_Ratios = length(AA_RatioList);

%% Initialize run

nRuns = nAA_Ratios*npK_Ratios;

% Initialize save data.
EqTint = zeros(nRuns,1);
EqText = zeros(nRuns,1);
Uptake = zeros(nRuns,1);
Efflux = zeros(nRuns,1);
Overwrite = zeros(nRuns,1);
RunList = zeros(nRuns,7);
Ratio_Runlist = zeros(nRuns,2);

% Initialize ContourW63G variables
x = {AA_RatioList, pK_RatioList};
Z = zeros(nAA_Ratios,npK_Ratios);

% Extra variables for bug testing
M = cell(nRuns,3);

%% Run the model

iter = 1;

for iAA = 1:nAA_Ratios
    
    curAA = AltK(AA_Matrix(iAA,1),AA_Matrix(iAA,2));
    
    for ipK = 1:npK_Ratios
        
        curpK = pK_Matrix(ipK,:);
        curK = calcTenStatesRates(curpK, curAA);
        [curT,~,curTint,curText, ~, ~, ~] = runTenState(curK, [0 1e8], pHs, [25e-9 25e-9], 20e-7, [1e-7 1e-7]);
        
        EqTint(iter) = calcEqT(curTint);
        EqText(iter) = calcEqT(curText);
        Uptake(iter) = EqTint(iter)/EqText(iter);
        Efflux(iter) = EqText(iter)/EqTint(iter);
        
        M{iter,1} = curT;
        M{iter,2} = curTint;
        M{iter,3} = curText;
        
        % Save the uptake data in those indices
        if ~(Z(iAA,ipK) == 0)
            Overwrite(iter) = 1;
        end
        
        Z(iAA,ipK) = Uptake(iter);
        RunList(iter,:) = [curpK, AA_Matrix(iAA,1), AA_Matrix(iAA,2)];
        Ratio_Runlist(iter,:) = [curpK(4)/curpK(5),AA_Matrix(iAA,1)/AA_Matrix(iAA,2)];
        
        iter = iter + 1;
        
    end
    
end

T1 = table((1:nRuns)',RunList(:,1),RunList(:,2),RunList(:,3),RunList(:,4),RunList(:,5),RunList(:,6),RunList(:,7));
T1.Properties.VariableNames = {'iter','pKa1','pKa2','pKa3','pKd1','pKd2','kAAED','kAAEHD'};

T2 = table((1:nRuns)',Ratio_Runlist(:,1),Ratio_Runlist(:,2),Uptake,Efflux,M(:,1),M(:,2));
T2.Properties.VariableNames = {'iter','pK_Ratio','AA_Ratio','Uptake','Efflux','t','Tint'};

end

