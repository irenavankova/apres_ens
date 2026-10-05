function [vsr_lin, vsr_lin_se, vsr_quad, vsr_quad_se, v_bed, v_bed_se, z_bed_comb] = fimf_vsr_multifreq_estimation(...
    sitename, bed_data, int_data, t_bed, t_int, f1, d2y, varargin)
% FIMF_VSR_ESTIMATION Multi-frequency VSR mean estimation.

    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'iNaN_array', [], @isnumeric);
    addParameter(p, 'opt_combine_method', 'median', @ischar);
    addParameter(p, 'opt_remove_outliers_vel', true, @islogical);
    addParameter(p, 'opt_vel_method', 'line_fit', @ischar);
    addParameter(p, 'opt_vel_method_deriv_bed', 'mean', @ischar);
    addParameter(p, 'opt_vel_method_deriv_int', 'median', @ischar);
    addParameter(p, 'opt_bed_source', 'tp', @ischar);
    addParameter(p, 'ylim_val', [0 170], @isnumeric);
    addParameter(p, 'ytick_dist', 20, @isnumeric);
    addParameter(p, 'drnf_vsr', [0 40; 124 4000], @isnumeric);
    addParameter(p, 'drnf_artefacts', [], @isnumeric);
    addParameter(p, 'plot_results', true, @islogical); 
    addParameter(p, 't_start', [], @isnumeric); % Optional explicit start time
    addParameter(p, 't_end', [], @isnumeric);   % Optional explicit end time
    addParameter(p, 'hax', [], @(x) true);      % Optional axes handle for overwriting
    addParameter(p, 'max_nan_threshold', 0.20, @isnumeric); % % Default 20% threshold ---
    parse(p, varargin{:});

    % Extract parameters
    iNaN_array = p.Results.iNaN_array;
    opt_combine_method = p.Results.opt_combine_method;
    opt_remove_outliers_vel = p.Results.opt_remove_outliers_vel;
    opt_vel_method = p.Results.opt_vel_method;
    opt_vel_method_deriv_bed = p.Results.opt_vel_method_deriv_bed;
    opt_vel_method_deriv_int = p.Results.opt_vel_method_deriv_int;
    opt_bed_source = p.Results.opt_bed_source;
    ylim_val = p.Results.ylim_val;
    ytick_dist = p.Results.ytick_dist;
    drnf_vsr = p.Results.drnf_vsr;
    drnf_artefacts = p.Results.drnf_artefacts;
    plot_results = p.Results.plot_results;
    hax = p.Results.hax;
    max_nan_threshold = p.Results.max_nan_threshold;

    % Determine bed data source
    is_bed_tp = strcmp(opt_bed_source, 'tp');

    % Number of frequency bands
    N = length(f1);

    % Initialize matrices to accumulate frequency instance data
    all_int_vel = [];
    all_bed_vel = [];
    all_bed_z = [];
    all_ampCor_mean = [];
    all_ampCor_std = [];

    % Time bounds (using explicitly passed or falling back to pre-loaded defaults)
    t1_use = p.Results.t_start;
    if isempty(t1_use), t1_use = t_int(1); end
    t2_use = p.Results.t_end;
    if isempty(t2_use), t2_use = t_int(end); end

    for k = 1:N
        % Get internal data
        ct = int_data{k};

        % Setup Time Bounds for this slice
        time = ct.dts_xcor.time;
        tind = find(time >= t1_use & time <= t2_use);
        
        if isempty(tind)
            continue; % Safety skip if window is entirely outside data bounds
        end

        % Get bed data
        bed_struct = bed_data{k};

        % Apply provided filtering (NaN exclusion)
        if ~isempty(iNaN_array)
            invalid_idx = intersect(tind, iNaN_array);
            
            % Abort if timeseries points exceed NaN threshold ---
            if (length(invalid_idx) / length(tind)) > max_nan_threshold
                continue; % Skips this window, leaving all_int_vel empty to trigger the NaN fallback
            end
            
            ct.dts_xcor.dh(invalid_idx, :) = NaN;
            if is_bed_tp
                bed_struct.thickness(invalid_idx) = NaN;
            else
                bed_struct.dh(invalid_idx) = NaN;
            end
        end

        % Store depth arrays
        z_all = ct.dts_xcor.dhRange;
        if is_bed_tp
            all_bed_z = cat(1, all_bed_z, bed_struct.x_input);
        else
            all_bed_z = cat(1, all_bed_z, bed_struct.dhRange);
        end

        % Calculate Mean Velocities for this frequency instance
        [v_int, v_bed_k] = fimf_calc_dhdt(ct.dts_xcor, bed_struct, tind, d2y, opt_vel_method, opt_vel_method_deriv_bed, opt_vel_method_deriv_int, opt_bed_source);

        % Accumulate results
        all_int_vel = cat(1, all_int_vel, v_int);
        all_bed_vel = cat(1, all_bed_vel, v_bed_k);

        % Ensure ampCor indices do not exceed its size
        tind_amp = tind(tind <= size(ct.dts_xcor.ampCor, 1));
        if ~isempty(tind_amp)
            all_ampCor_mean = cat(1, all_ampCor_mean, mean(ct.dts_xcor.ampCor(tind_amp, :), 1, 'omitnan'));
            all_ampCor_std = cat(1, all_ampCor_std, std(ct.dts_xcor.ampCor(tind_amp, :), 1, 'omitnan'));
        end
    end

    % Safe fallback if window resulted in empty extraction
    if isempty(all_int_vel)
        vsr_lin = NaN; vsr_lin_se = NaN; vsr_quad = NaN; vsr_quad_se = NaN; 
        v_bed = NaN; v_bed_se = NaN; z_bed_comb = NaN;
        return;
    end

    %% Statistics across central frequencies
    if opt_remove_outliers_vel
        ind_out_int = isoutlier(all_int_vel, 'median', 1);
        all_int_vel(ind_out_int) = NaN;

        ind_out_bed = isoutlier(all_bed_vel, 'median', 1);
        all_bed_vel(ind_out_bed) = NaN;

        ind_out_bed_z = isoutlier(all_bed_z, 'median', 1);
        all_bed_z(ind_out_bed_z) = NaN;
    end

    % Apply user-specified combine method for VSR estimation and bed positioning
    if strcmp(opt_combine_method, 'median')
        v_all = median(all_int_vel, 1, 'omitnan');
        v_bed = median(all_bed_vel, 1, 'omitnan');
        z_bed_comb = median(all_bed_z, 1, 'omitnan');
    else
        v_all = mean(all_int_vel, 1, 'omitnan');
        v_bed = mean(all_bed_vel, 1, 'omitnan');
        z_bed_comb = mean(all_bed_z, 1, 'omitnan');
    end

    % Bed rate of change standard error across frequency estimates
    v_bed_se = std(all_bed_vel, 0, 1, 'omitnan') ./ sqrt(N);

    % Use standard deviation divided by sqrt(N) for ensemble error weights
    v_all_se = std(all_int_vel, 0, 1, 'omitnan') ./ sqrt(N);
    
    if ~isempty(all_ampCor_mean)
        mean_ampCor = mean(all_ampCor_mean, 1, 'omitnan');
        std_ampCor = mean(all_ampCor_std, 1, 'omitnan');
    else
        mean_ampCor = []; std_ampCor = [];
    end

    %% VSR Calculation and Fitting
    [ind_fit, q_best, q_se, p_best, p_best_se] = fimf_calc_vsr(z_all', v_all', v_all_se', drnf_vsr);

    % --- NEW: Catch NaN fits, assign NaN outputs, and early-exit to skip plotting ---
    if any(isnan(q_best)) || any(isnan(p_best))
        vsr_lin = NaN; vsr_lin_se = NaN; 
        vsr_quad = NaN; vsr_quad_se = NaN;
        return;
    end

    %% Final Plot Generation and VSR calculation at bed depth
    [vsr_lin, vsr_lin_se, vsr_quad, vsr_quad_se] = fimf_plot_vsr(sitename, z_all, all_bed_z, v_all, all_int_vel, ...
        v_bed, z_bed_comb, all_bed_vel, v_all_se, ind_fit, ...
        q_best, q_se, p_best, p_best_se, all_ampCor_mean, ...
        mean_ampCor, std_ampCor, drnf_artefacts, ylim_val, ytick_dist, plot_results, hax);
end