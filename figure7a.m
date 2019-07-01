function [T,fig] = figure7a(pHs)

Poss_AA = round(logspace(0,2,100),1);  % Logspace vector from 1 to 100 rounded to the nearest 10ths place.

Conds = [Poss_AA',flipud(Poss_AA'),Poss_AA',flipud(Poss_AA')];   % Create conditions as described in figure 6

pK = [7 7 7 7];

[nRuns, ~] = size(Conds);

EqTint = zeros(nRuns,1);
EqText = zeros(nRuns,1);
Uptake = zeros(nRuns,1);
Efflux = zeros(nRuns,1);

for iter = 1:nRuns
        
    curAA = [Conds(iter,1);...                 % k17 PROTON
            Conds(iter,1);...                  % k18 PROTON
            Conds(iter,2);...                  % k19 APO
            Conds(iter,2);...                  % k20 APO
            Conds(iter,3);...                  % k21 DRUG
            Conds(iter,3);...                  % k22 DRUG
            Conds(iter,4);...                  % k23 DUO
            Conds(iter,4)];
    
    curK = calcEightStateRates(pK, curAA);
    [~,~,curTint,curText, ~, ~, ~] = runEightState(curK, [1 1e8], pHs, [25e-9 25e-9], 20e-7, [1e-7 1e-7]);
    
    EqTint(iter) = calcEqT(curTint);
    EqText(iter) = calcEqT(curText);
    Uptake(iter) = EqTint(iter)/EqText(iter);
    Efflux(iter) = EqText(iter)/EqTint(iter);
        
end

T = table(Conds(:,1),Conds(:,2),Conds(:,3),Conds(:,4),Conds(:,1)./Conds(:,2),Uptake,Efflux);
T.Properties.VariableNames = {'EH','E','ED','EHD','AA_Ratio','Uptake','Efflux'};

fig = figure;
plot(T.AA_Ratio,T.Uptake);
set(gca,'XScale', 'log','YScale', 'log');
xlabel('R_{AA}');
ylabel('T_{r}');

end




