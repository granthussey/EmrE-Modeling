function [T,fig] = figure7c(pHs)

% Initialize alternating access possibilities
Poss_AA = [1, 2, 3, 4, 5];
AAconds = [Poss_AA',flipud(Poss_AA'),Poss_AA',flipud(Poss_AA')];
[nAA,~] = size(AAconds);

% Initialize ranges for pKa and pKd values
Poss_pKa = (4:.01:10)';
Poss_pKd = (3:.01:9)';
pKConds = [Poss_pKa,flipud(Poss_pKa),Poss_pKd,flipud(Poss_pKd)];
[npK,~] = size(pKConds);

% Total number of simulations are a product of both varibles to iterate
% over.
nSim = npK*nAA;

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
graph_AARatios   = zeros(npK,nAA);

iter = 1;

for iAA = 1:nAA
    
    curAA = [AAconds(iAA,1);...             % k17 PROTON
        AAconds(iAA,1);...                  % k18 PROTON
        AAconds(iAA,2);...                  % k19 APO
        AAconds(iAA,2);...                  % k20 APO
        AAconds(iAA,3);...                  % k21 DRUG
        AAconds(iAA,3);...                  % k22 DRUG
        AAconds(iAA,4);...                  % k23 DUO
        AAconds(iAA,4)];
    
    for ipK = 1:npK
        
        curpK = pKConds(ipK,:);
        
        curK = calcEightStateRates(curpK, curAA);
        [~,~,curTint,curText, ~, ~, ~] = runEightState(curK, [1 1e8], pHs, [25e-9 25e-9], 20e-8, [5e-8 5e-8]);
        
        EqTint(iter) = calcEqT(curTint);
        EqText(iter) = calcEqT(curText);
        Uptake(iter) = EqTint(iter)/EqText(iter);
        Efflux(iter) = EqText(iter)/EqTint(iter);
        
        pKs(iter,:)  = pKConds(ipK,:);
        AAs(iter,:)  = [AAconds(iAA,1),AAconds(iAA,2),AAconds(iAA,3),AAconds(iAA,4)];
                
        graphEqTint(ipK,iAA) = EqTint(iter);
        graphEqText(ipK,iAA) = EqText(iter);
        
        graphUptake(ipK,iAA) = Uptake(iter);
        graphEfflux(ipK,iAA) = Efflux(iter);
        
        GRAPH_kOff_pK1 = 1e10.*10.^(-1.*pKConds(ipK,1));
        GRAPH_kOff_pK2 = 1e10.*10.^(-1.*pKConds(ipK,2));
        
        graph_Ratio_kOff(ipK,iAA) = GRAPH_kOff_pK1./GRAPH_kOff_pK2;
        graph_AARatios(ipK,iAA) = AAconds(iAA,1)./AAconds(iAA,2);
        
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
legend('R_{AA} = .2','R_{AA} = .5','R_{AA} = 1','R_{AA} = 2','R_{AA} = 5')

end



