function [iNaN, iNaN_bed, iNaN_int] = fimf_filter_bw(dhdtmat_bed, dhdtmat_int, t_bed, t_int, Nbad, med_filt_param, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, plot_results)
%FIMF_FILTER_BW Apply filtering and jump detection to dh/dt timeseries.

    % Initialize filtered outputs
    y_bed = [];
    y_int = [];

    % Process bed data
    if ~isempty(dhdtmat_bed)
        N = size(dhdtmat_bed, 2);
        dhmean_bed = y_bed_raw;

        % Jump detection for bed data - ordered BW criterion
        diff_dhdt_bed = diff(dhdtmat_bed, [], 2); % difference all realizations in order of BW 
        sumdiff_bed = abs(sum(sign(diff_dhdt_bed), 2)); % convert to sign and sum it up

        %i_jumpN_bed = sumdiff_bed >= (N - Nbad);
        %max_dhdt_bed = max(abs(dhdtmat_bed), [], 2);

        sumdiffdh_bed = abs(sum(sign(dhdtmat_bed), 2)); %not used - just plotted

        i_jumpNmax_bed = sumdiff_bed >= (N - Nbad) & (abs(dhmean_bed) > med_filt_param * median(abs(dhmean_bed), 'omitnan')); %if most (N-Nbad) BW differences same sign and the BW mean large enough, throw away

        % Apply masks
        nan_indices_bed = find(i_jumpNmax_bed);
        all_indices_bed = unique([nan_indices_bed; nan_indices_bed - 1; nan_indices_bed + 1]);
        all_indices_bed = all_indices_bed(all_indices_bed >= 1 & all_indices_bed <= length(dhmean_bed));
        dhmean_ijumpNnext_bed = dhmean_bed;
        dhmean_ijumpNnext_bed(all_indices_bed) = NaN;

        y_bed = dhmean_ijumpNnext_bed;
    end

    % Process internal data
    if ~isempty(dhdtmat_int)
        dhmean_int = y_int_raw;

        % Jump detection for internal data - large jump on one side of axis
        % for all BW
        diff_dhdt_int = diff(dhdtmat_int, [], 2); % difference all realizations in order of BW
        sumdiff_int = abs(sum(sign(diff_dhdt_int), 2)); % convert to sign and sum it up

        % This is actually used later on
        sumdiffdh_int = abs(sum(sign(dhdtmat_int), 2)); % Look at the sign of the BW realizations and sum it up
        all_pos_sign = all(dhdtmat_int - median(dhmean_int) > 0, 2);
        all_neg_sign = all(dhdtmat_int - median(dhmean_int) < 0, 2);
        i_jumpsignN = (all_pos_sign | all_neg_sign) & (abs(dhmean_int) > med_filt_param * median(abs(dhmean_int), 'omitnan')); % if all BW same sign and their mean large enough, throw away

        % Apply masks
        dhmean_ijumpsignN = dhmean_int;
        dhmean_ijumpsignN(i_jumpsignN) = NaN;

        y_int = dhmean_ijumpsignN;
    end

    % Compute individual and combined NaN masks
    iNaN_bed = isnan(y_bed);
    iNaN_int = isnan(y_int);
    iNaN = iNaN_bed | iNaN_int;

    % Plot results if enabled
    if plot_results
        Lwide = 1;
        nplot = 3;

        % Reconstruct combo for plotting
        y_bed_combo_plot = y_bed_raw;
        y_bed_combo_plot(iNaN) = NaN;
        y_int_combo_plot = y_int_raw;
        y_int_combo_plot(iNaN) = NaN;

        % --- FIGURE 1: Bed data (realizations, differences, sumdiff) ---
        if ~isempty(dhdtmat_bed)
            figure;
            h = [];
            colors = pink(size(dhdtmat_bed, 2) + 1);

            subplot(nplot, 1, 1); h = gca;
            for k = 1:size(dhdtmat_bed, 2)
                plot(t_bed(1:end-1), dhdtmat_bed(:, k), 'Color', colors(k, :)); hold on;
            end
            plot(t_bed(1:end-1), y_bed_raw, 'k', 'LineWidth', 2);
            plot(t_bed(1:end-1), y_bed, 'r', 'LineWidth', 2);
            ylabel('Realizations');
            datetick('x', 'mm/yy', 'keeplimits');
            title(['Bandwidth Filter, Nbad = ' num2str(Nbad)], 'fontweight', 'normal');

            subplot(nplot, 1, 2); h = [h; gca];
            for k = 1:size(diff_dhdt_bed, 2)
                plot(t_bed(1:end-1), diff_dhdt_bed(:, k), 'Color', colors(k, :)); hold on;
            end
            ylabel('Differences');
            datetick('x', 'mm/yy', 'keeplimits');

            subplot(nplot, 1, 3); h = [h; gca];
            plot(t_bed(1:end-1), sumdiff_bed, 'b-', 'LineWidth', 1); hold on;
            plot(t_bed(1:end-1), sumdiffdh_bed, 'r-', 'LineWidth', 1);
            ylabel('Sum of Differences');
            datetick('x', 'mm/yy', 'keeplimits');

            linkaxes(h, 'x');
            legend('sumdiff', 'sumdiffdh');
        end

        % --- FIGURE 2: Internal data (realizations, differences, sumdiff) ---
        if ~isempty(dhdtmat_int)
            figure;
            h = [];
            colors = pink(size(dhdtmat_int, 2) + 1);

            subplot(nplot, 1, 1); h = gca;
            for k = 1:size(dhdtmat_int, 2)
                plot(t_int(1:end-1), dhdtmat_int(:, k), 'Color', colors(k, :)); hold on;
            end
            plot(t_int(1:end-1), y_int_raw, 'k', 'LineWidth', 2);
            plot(t_int(1:end-1), y_int, 'r', 'LineWidth', 2);
            ylabel('Realizations');
            datetick('x', 'mm/yy', 'keeplimits');
            title(['Bandwidth Filter, Nbad = ' num2str(Nbad)], 'fontweight', 'normal');

            subplot(nplot, 1, 2); h = [h; gca];
            for k = 1:size(diff_dhdt_int, 2)
                plot(t_int(1:end-1), diff_dhdt_int(:, k), 'Color', colors(k, :)); hold on;
            end
            ylabel('Differences');
            datetick('x', 'mm/yy', 'keeplimits');

            subplot(nplot, 1, 3); h = [h; gca];
            plot(t_int(1:end-1), sumdiff_int, 'b-', 'LineWidth', 1); hold on;
            plot(t_int(1:end-1), sumdiffdh_int, 'r-', 'LineWidth', 1);
            ylabel('Sum of Differences');
            datetick('x', 'mm/yy', 'keeplimits');

            linkaxes(h, 'x');
            legend('sumdiff', 'sumdiffdh');
        end

        % --- FIGURE 3: Raw, filtered, and combined timeseries ---
        figure;
        subplot(2, 1, 1); h = gca;
        plot(t_bed(1:end-1), y_bed_raw, 'k-', 'LineWidth', Lwide); hold on;
        plot(t_bed(1:end-1), y_bed, 'r-', 'LineWidth', Lwide);
        plot(t_bed(1:end-1), y_bed_combo_plot, 'c-', 'LineWidth', Lwide);
        grid on;
        ylabel('Bed');
        datetick('x', 'mm/yy', 'keeplimits');
        title(['Bandwidth Filter, Nbad = ' num2str(Nbad)], 'fontweight', 'normal');

        subplot(2, 1, 2); h = [h; gca];
        plot(t_int(1:end-1), y_int_raw, 'k-', 'LineWidth', Lwide); hold on;
        plot(t_int(1:end-1), y_int, 'r-', 'LineWidth', Lwide);
        plot(t_int(1:end-1), y_int_combo_plot, 'c-', 'LineWidth', Lwide);
        linkaxes(h, 'x');
        datetick('x', 'mm/yy', 'keeplimits');
        grid on;
        ylabel('Internal');

        % --- FIGURE 4: Uncertainty bounds ---
        figure;
        subplot(2, 1, 1); h = gca;
        plot(t_bed(1:end-1), y_bed_raw - y_bed_std, 'c-', 'LineWidth', Lwide); hold on;
        plot(t_bed(1:end-1), y_bed_raw + y_bed_std, 'c-', 'LineWidth', Lwide);
        plot(t_bed(1:end-1), y_bed, 'r-', 'LineWidth', Lwide);
        plot(t_bed(1:end-1), y_bed_combo_plot, 'k-', 'LineWidth', Lwide);
        grid on;
        ylabel('Bed');
        datetick('x', 'mm/yy', 'keeplimits');
        ylim([-ymax ymax]);
        title(['Bandwidth Filter, Nbad = ' num2str(Nbad)], 'fontweight', 'normal');

        subplot(2, 1, 2); h = [h; gca];
        plot(t_int(1:end-1), y_int_raw - y_int_std, 'c-', 'LineWidth', Lwide); hold on;
        plot(t_int(1:end-1), y_int_raw + y_int_std, 'c-', 'LineWidth', Lwide);
        plot(t_int(1:end-1), y_int, 'r-', 'LineWidth', Lwide);
        plot(t_int(1:end-1), y_int_combo_plot, 'k-', 'LineWidth', Lwide);
        linkaxes(h, 'x');
        datetick('x', 'mm/yy', 'keeplimits');
        grid on;
        ylabel('Internal');
        ylim([-ymax ymax]);
    end
end