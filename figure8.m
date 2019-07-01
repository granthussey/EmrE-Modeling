function [T,fig] = figure8(pHs)

% Initialize alternating access matrix
AltK = @(x) [x;x;x;x;x;x;x;x;x;x];
AlternatingK = AltK(1);

% Initialize ranges for pKa and pKd values
spacer = .1;
pKa_Range = (4:spacer:10)';
pKd_Range = (3:spacer:9)';
npK = length(pKa_Range);

% Arrange pKa and pKd values for running simulations
A = [pKa_Range,flipud(pKa_Range),(zeros(1,npK))',pKd_Range,flipud(pKd_Range)];
Check = (A(:,1) + A(:,5) == A(:,2) + A(:,4));
pK_Matrix = A(Check,:);         % Makes sure that all pKa/pKd values are thermodynamically feasible

% Initialize the pKa3 values to be used
pKa3_Range = [4, 5, 7, 7.5, 8, 9];
npKa3 = length(pKa3_Range);

nRuns = npK * npKa3;            % Total number of simulations are a product of npK
                                % values for Roff scan and npKa3 values to be tested.

% Initialize save data indexed by simulation number
EqTint  = zeros(nRuns,1);
EqText  = zeros(nRuns,1);
Uptake  = zeros(nRuns,1);
Efflux  = zeros(nRuns,1);
Runlist = zeros(nRuns,5);

% Initialize save data indexed by Roff and pKa3 value
Z_Uptake = zeros(npK,npKa3);
Z_Roff   = zeros(npK,npKa3);
Z_pKa1   = zeros(npK,npKa3);

% nested for loop to run simulations
iter = 1;
for ipKa3 = 1:npKa3
    
    curpKa3 = pKa3_Range(ipKa3);
    
    for ipK = 1:npK
        
        curpK = pK_Matrix(ipK,:);
        curpK(3) = curpKa3;
        
        curK = calcTenStatesRates(curpK, AlternatingK);
        
        [~,~,curTint,curText, ~, ~, ~] = runTenState(curK, [0 1e8], pHs, [25e-9 25e-9], 20e-7, [1e-7 1e-7]);
        
        EqTint(iter) = calcEqT(curTint);
        EqText(iter) = calcEqT(curText);
        Uptake(iter) = EqTint(iter)/EqText(iter);
        Efflux(iter) = EqText(iter)/EqTint(iter);
        
        curkOff_pK1 = 1e10*10^(-1*curpK(1));
        curkOff_pK2 = 1e10*10^(-1*curpK(2));
        curRoff = curkOff_pK1/curkOff_pK2;
        
        Z_Uptake(ipK,ipKa3) = Uptake(iter);
        Z_Roff(ipK,ipKa3)   = curRoff;
        Z_pKa1(ipK,ipKa3)   = curpK(1);
        
        Runlist(iter,:) = curpK;
        
        iter = iter + 1;
        
    end
    
    
end

curkOff_pK1 = 1e10.*10.^(-1.*Runlist(:,1));
curkOff_pK2 = 1e10.*10.^(-1.*Runlist(:,2));
Ratio_kOff = curkOff_pK1./curkOff_pK2;

T = table(Runlist(:,1),Runlist(:,2),Runlist(:,3),Runlist(:,4),Runlist(:,5),Ratio_kOff,Uptake,Efflux);
T.Properties.VariableNames = {'pKa1','pKa2','pKa3','pKd1','pKd2','Ratio_kOff','Uptake','Efflux'};


% Initialize a matrix to remove simulations where pKa^{EH} < pKa^{EHH}
Logic = zeros(npK,npKa3);

for iter = 1:npKa3
    Logic(:,iter) = (Z_pKa1(:,iter) >= pKa3_Range(iter)); % Fill in that matrix
end

Z_Roff   = Z_Roff.*Logic;
Z_Uptake = Z_Uptake.*Logic;

Z_Roff(Z_Roff == 0) = NaN;
Z_Uptake(Z_Uptake == 0) = NaN;

% Plot each trace for its respective pKa3 value

fig = figure;
for iter = 1:npKa3
    loglog(Z_Roff(:,iter),Z_Uptake(:,iter));
    hold on
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('kOff');
ylabel('Uptake');
legend('pKa = 4','pKa = 5','pKa = 7','pKa = 7.5','pKa = 8','pKa = 9');


end