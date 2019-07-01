

% Initialize figure
figure;

% Get data for 10-state trace
[T] = figure9_10state;

% Plot
semilogy(T.pHs(:,2),T.Efflux);
ylabel('[drug]_{ext} /[drug]_{int}');
xlabel('External pH');
hold on

% Get data for 8-state trace
[T] = figure9_8state;

% Plot
semilogy(T.pHs(:,2),T.Efflux);
ylabel('[drug]_{ext} /[drug]_{int}');
xlabel('External pH');

% Fix the x-axis
xlim([5.4,9.4]);