function summary_table = algo_parameters_test(resids, bins_list, eps_list)
    % resids: table with standardized uniform residuals (columns: ITA, GER, FRA, SPA)

    real_data = table2array(resids(1:1500, :));    
    summary = [];
    keys = {};
    simulations = {}; % Store only top 3 simulations
    idx = 1;

    for b = 1:length(bins_list)
        bins = bins_list(b);
        for e = 1:length(eps_list)
            epsilon = eps_list(e);
            fprintf('Running combination: bins = %d, epsilon = %.4f\n', bins, epsilon);

            % Step 1: compute empirical marginals
            bins1 = linspace(min(resids.ITA), max(resids.ITA), bins);
            bins2 = linspace(min(resids.GER), max(resids.GER), bins);
            bins3 = linspace(min(resids.FRA), max(resids.FRA), bins);
            bins4 = linspace(min(resids.SPA), max(resids.SPA), bins);

            mu1 = histcounts(resids.ITA, [bins1, inf], 'Normalization', 'probability')';
            mu2 = histcounts(resids.GER, [bins2, inf], 'Normalization', 'probability')';
            mu3 = histcounts(resids.FRA, [bins3, inf], 'Normalization', 'probability')';
            mu4 = histcounts(resids.SPA, [bins4, inf], 'Normalization', 'probability')';

            % Step 2: build cost matrix
            [X1, X2, X3, X4] = ndgrid(bins1, bins2, bins3, bins4);
            C = (X1 - X2).^2 + (X1 - X3).^2 + (X1 - X4).^2 + ...
                (X2 - X3).^2 + (X2 - X4).^2 + (X3 - X4).^2;

            % Step 3: run Sinkhorn algorithm
            [~, pi_matrix] = sinkhorn_mot_4(mu1, mu2, mu3, mu4, C, epsilon, 20, 1e-5);

            % Step 4: simulate 100000 samples (used for statistics only)
            coarse_inds = randsample(numel(pi_matrix), 100000, true, pi_matrix(:));
            [ci, cj, ck, cl] = ind2sub(size(pi_matrix), coarse_inds);

            bin_width1 = (max(resids.ITA) - min(resids.ITA)) / bins;
            bin_width2 = (max(resids.GER) - min(resids.GER)) / bins;
            bin_width3 = (max(resids.FRA) - min(resids.FRA)) / bins;
            bin_width4 = (max(resids.SPA) - min(resids.SPA)) / bins;

            samples = [
                bins1(ci)' + rand(100000, 1) * bin_width1, ...
                bins2(cj)' + rand(100000, 1) * bin_width2, ...
                bins3(ck)' + rand(100000, 1) * bin_width3, ...
                bins4(cl)' + rand(100000, 1) * bin_width4];

            % Step 5: compute tau score on full simulated data
            tau_real = corr(real_data, 'type', 'Kendall');
            tau_sim  = corr(samples, 'type', 'Kendall');
            tau_diff = abs(tau_real - tau_sim);
            tau_score = sum(tau_diff(:));

            % Step 6: compute moment score (on full simulation)
            moment_score = 0;
            for v = 1:4
                diff_mean = abs(mean(real_data(:,v)) - mean(samples(:,v)));
                diff_var  = abs(var(real_data(:,v))  - var(samples(:,v)));
                diff_skew = abs(skewness(real_data(:,v)) - skewness(samples(:,v)));
                diff_kurt = abs(kurtosis(real_data(:,v)) - kurtosis(samples(:,v)));
                moment_score = moment_score + diff_mean + diff_var + diff_skew + diff_kurt;
            end

            % Save key, raw scores, and simulation
            keys{idx,1} = sprintf('%d bins - %.4f eps', bins, epsilon);
            summary(idx, :) = [tau_score, moment_score];
            simulations{idx} = samples;
            idx = idx + 1;
        end
    end

    % Normalize scores
    tau_raw = summary(:,1);
    moment_raw = summary(:,2);

    tau_norm = (tau_raw - min(tau_raw)) ./ (max(tau_raw) - min(tau_raw));
    moment_norm = (moment_raw - min(moment_raw)) ./ (max(moment_raw) - min(moment_raw));

    % Compute final z_score
    z_score = tau_norm + moment_norm;

    % Build summary table
    summary_table = table(keys, tau_norm, moment_norm, z_score, ...
        'VariableNames', {'Parameters', 'TAU_SCORE_NORM', 'MOMENT_SCORE_NORM', 'Z_SCORE'});

    % Select and export best 3
    [~, best_idx] = mink(z_score, 3);
    fprintf('\nTop 3 parameter combinations (models saved):\n');
    fileID = fopen('top3_models.txt', 'w');
    fprintf(fileID, 'Top 3 parameter combinations (saved as CSV):\n');

    for i = 1:3
        fprintf('%d) %s\n', i, summary_table.Parameters{best_idx(i)});
        fprintf(fileID, '%d) %s\n', i, summary_table.Parameters{best_idx(i)});

        % Export top simulation
        samples_top = simulations{best_idx(i)};
        T = array2table(samples_top, 'VariableNames', {'ITA', 'GER', 'FRA', 'SPA'});
        filename_export = sprintf('sinkhorn_simulations_TOP%d.csv', i);
        writetable(T, filename_export);
    end
    fclose(fileID);
end