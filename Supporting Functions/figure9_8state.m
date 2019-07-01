function [T] = figure9_8state

%% Setup

AltK_Anti = [100;...           % k17 PROTON
    100;...                    % k18 PROTON
    10;...                     % k19 APO
    10;...                     % k20 APO
    100;...                    % k21 DRUG
    100;...                    % k22 DRUG
    10;...                     % k23 DUO
    10];                       % k24 DUO

pKs = [8, 7, 8, 7];

k = calcEightStateRates(pKs,AltK_Anti);

Times = [0 1e8];
Drugs = [25e-9 25e-9];
EmrE = 20e-7;
Volumes = [1e-7 1e-7];

%% Initialize run parameters

pHint = 7.4;
pHext = (5.4:.05:9.4)';

nRuns = length(pHext);

pHs = [pHint.*ones(nRuns,1), pHext];

dpHs = pHs(:,1) - pHs(:,2);

nRuns = length(dpHs);

EqTint = zeros(nRuns,1);
EqText = zeros(nRuns,1);
Uptake = zeros(nRuns,1);
Efflux = zeros(nRuns,1);

%% Begin runs

for iter = 1:nRuns
    
    curpH = pHs(iter,:);
    [~,~,curTint,curText, ~, ~, ~] = runEightState(k, Times, curpH, Drugs, EmrE, Volumes);
    
    EqTint(iter) = calcEqT(curTint);
    EqText(iter) = calcEqT(curText);
    Uptake(iter) = EqTint(iter)/EqText(iter);
    Efflux(iter) = EqText(iter)/EqTint(iter);
    
end

%% Summarize data

data = Efflux;  
edges = fliplr([inf, 1, -inf]);
binData = discretize(data,edges);
PercentMax = zeros(nRuns,1);

% Bin 1 is for those with Efflux < 1, must use Uptake instead.
idx = binData == 1;
PercentMax(idx) = 100.*(Uptake(idx)./(10.^(-dpHs(idx))));

% Bin 2 is for those with Efflux > 1, use Efflux.
idx = binData == 2;
PercentMax(idx) = 100.*(Efflux(idx)./(10.^(dpHs(idx))));

T = table(pHs,dpHs,PercentMax,Uptake,Efflux);
T.Properties.VariableNames = {'pHs','dpH','PercMax','Uptake','Efflux'};

end
