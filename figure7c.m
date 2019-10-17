

DrugRange = [1e-9, 400e-9, 135e-6, 50e-3];
nDrug = length(DrugRange);

Data = cell(1,nDrug);

% Run
for iter = 1:nDrug
    
    [~,graph] = figure7b_VaryDrug([6.5,7.5], DrugRange(iter), 1);
    Data{iter} = graph;
    
end

% Graph
figure;
for iter = 1:nDrug
    
    curGraph = Data{iter};
    loglog(curGraph.graph_Ratio_kOff,curGraph.graphUptake);
    hold on
    
end


% Format graph
xlabel('R_{off}');
ylabel('T_{r}');
ylim([0.1, 10]);

%rang = num2cell(DrugRange);
%leg = cellfun(@(x) num2str(x),rang,'UniformOutput',false);
%legend(leg);

legend('1 nM','400 nM','135 uM','50 mM');

title('with k_{AA} = 1');



DrugRange = [1e-9, 400e-9, 135e-6, 50e-3];
nDrug = length(DrugRange);

Data = cell(1,nDrug);

% Run
for iter = 1:nDrug
    
    [~,graph] = figure7b_VaryDrug([6.5,7.5], DrugRange(iter), 1000);
    Data{iter} = graph;
    
end

% Graph
figure;
for iter = 1:nDrug
    
    curGraph = Data{iter};
    loglog(curGraph.graph_Ratio_kOff,curGraph.graphUptake);
    hold on
    
end


% Format graph
xlabel('R_{off}');
ylabel('T_{r}');
ylim([0.1, 10]);

%rang = num2cell(DrugRange);
%leg = cellfun(@(x) num2str(x),rang,'UniformOutput',false);
%legend(leg);

legend('1 nM','400 nM','135 uM','50 mM');
title('with k_{AA} = 1000');

