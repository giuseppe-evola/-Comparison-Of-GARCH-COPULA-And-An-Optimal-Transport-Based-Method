data = readtable('standardized_residuals.csv');
data.Properties.VariableNames = {'ITA', 'GER', 'FRA', 'SPA'};

bins_list = [30, 40, 50, 75];        
eps_list = [0.001, 0.002, 0.005, 0.02];  

summary_table = algo_parameters_test(data, bins_list, eps_list);
writetable(summary_table, 'parameter_tuning_summary.csv');