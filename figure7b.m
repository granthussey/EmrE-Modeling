function [T,fig] = figure7b(pHs)

% Initialize variables to test three uniformly-constant AA rates 
AltK = @(x) [x;x;x;x;x;x;x;x];
Poss_AA = [1, 10, 100];
[~,nAA] = size(Poss_AA);

% Initialize ranges for pKa and pKd values
Poss_pKa = (4:.05:10)';
Poss_pKd = (3:.05:9)';
pKConds = [Poss_pKa,flipud(Poss_pKa),Poss_pKd,flipud(Poss_pKd)];
[npK,~] = size(pKConds);

nSim = npK*nAA;                 % Total number of simulations are a product of unique Roff
                                % values and alternating access values to be tested.
                                
% Initialize save data indexed by simulation number
EqTint = zeros(nSim,1);
EqText = zeros(nSim,1);
Uptake = zeros(nSim,1);
Efflux = zeros(nSim,1);
pKs    = zeros(nSim,4);
AAs    = zeros(nSim,4);

% Initialize save data indexed by Roff and AA value
graphEqTint = zeros(npK,nAA);
graphEqText = zeros(npK,nAA);
graphUptake = zeros(npK,nAA);
graphEfflux = zeros(npK,nAA);
graph_Ratio_kOff = zeros(npK,nAA);

iter = 1;
for iAA = 1:nAA
    
    curAA = AltK(Poss_AA(iAA));
    
    for ipK = 1:npK
        
        curpK = pKConds(ipK,:);
        curK = calcEightStateRates(curpK, curAA);
        
        [~,~,curTint,curText, ~, ~, ~] = runEightState(curK, [1 1e8], pHs, [25e-9 25e-9], 20e-7, [1e-7 1e-7]);
        
        EqTint(iter) = calcEqT(curTint);
        EqText(iter) = calcEqT(curText);
        Uptake(iter) = EqTint(iter)/EqText(iter);
        Efflux(iter) = EqText(iter)/EqTint(iter);
        
        pKs(iter,:)  = pKConds(ipK,:);
        AAs(iter,:)  = [curAA(1),curAA(3),curAA(5),curAA(6)];
                
        graphEqTint(ipK,iAA) = EqTint(iter);
        graphEqText(ipK,iAA) = EqText(iter);
        
        graphUptake(ipK,iAA) = Uptake(iter);
        graphEfflux(ipK,iAA) = Efflux(iter);
        
        GRAPH_kOff_pK1 = 1e10.*10.^(-1.*pKConds(ipK,1));
        GRAPH_kOff_pK2 = 1e10.*10.^(-1.*pKConds(ipK,2));
        
        graph_Ratio_kOff(ipK,iAA) = GRAPH_kOff_pK1./GRAPH_kOff_pK2;
        
        iter = iter + 1;
        
    end
end

kOff_pK1 = 1e10.*10.^(-1.*pKs(:,1));
kOff_pK2 = 1e10.*10.^(-1.*pKs(:,2));
Ratio_kOff = kOff_pK1./kOff_pK2;

T = table(AAs(:,1),AAs(:,2),AAs(:,3),AAs(:,4),pKs(:,1),pKs(:,2),pKs(:,3),pKs(:,4),Ratio_kOff,AAs(:,1)./AAs(:,2),pKs(:,1)-pKs(:,2),Uptake,Efflux);
T.Properties.VariableNames = {'EH','E','ED','EHD','pKa1','pKa2','pKd1','pKd2','Ratio_kOff','AA_Ratio','pKDiff','Uptake','Efflux'};

fig = figure;
loglog(graph_Ratio_kOff,graphUptake)
xlabel('R_{off}')
ylabel('T_{r}')
legend('k_{AA} = 1','k_{AA} = 10','k_{AA} = 100');

end



