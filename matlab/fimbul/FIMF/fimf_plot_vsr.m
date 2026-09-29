function fimf_plot_vsr(sitename, z_all, all_bed_z, v_all, all_int_vel, v_bed_comb, z_bed_comb, all_bed_vel, v_all_se, ind_fit, q_best, q_se, p_best, p_best_se, all_ampCor_mean, mean_ampCor, std_ampCor, drnf_artefacts, ylim_val, ytick_dist)
    % Dedicated plotting function (Full extent line fits & purple cross for bed mean)
    
    [hax,plotwidth] = ap_sub_sub(2,1,3,10,0.25,0);
    
    % Weights for Ftest
    w = 1./v_all_se;
    w = w/max(w);

    %% hax(1): Velocity Fits (Top Plot)
    set(gcf,'CurrentAxes',hax(1))
    
    % Separate indices into used vs not used
    ind_not_fit = setdiff(1:length(z_all), ind_fit);
    
    % 1. Plot ALL individual internal velocity estimates across frequencies
    for k = 1:size(all_int_vel, 1)
        plot(all_int_vel(k, ind_not_fit), z_all(ind_not_fit), '.', 'color', [0.9 0.8 0]); hold on
        plot(all_int_vel(k, ind_fit), z_all(ind_fit), '.', 'color', [0.5 0.8 1]); hold on
    end
    
    % Plot ALL individual basal velocity estimates at their respective ranges
    if ~isempty(all_bed_vel) && ~isempty(all_bed_z)
        for k = 1:size(all_bed_vel, 1)
            plot(all_bed_vel(k), all_bed_z(k), '.', 'color', [0.8 0.6 1]); hold on % Light purple
        end
    end
    
    % 2. Plot MEAN/MEDIAN velocity estimates 
    plot(v_all(ind_not_fit), z_all(ind_not_fit), '.', 'color', [0.9 0.4 0], 'MarkerSize', 10); hold on
    plot(v_all(ind_fit), z_all(ind_fit), '.', 'color', [0 0 0.5], 'MarkerSize', 10); hold on
    
    % Basal mean/median plotted as a purple cross (+)
    if ~isempty(v_bed_comb) && ~isempty(z_bed_comb)
        plot(v_bed_comb, z_bed_comb, '+', 'color', [0.5 0 0.5], 'LineWidth', 1.5, 'MarkerSize', 8); hold on 
    end

    % --- Extended depth vector spanning the full plotting extent for fits ---
    z_full = linspace(ylim_val(1), ylim_val(2), 500)';

    % Apply Quadratic Function over full extent
    q = q_best + q_se; yQL = ct_tide_compute_quadratic(z_full, q); 
    q = q_best - q_se; yQU = ct_tide_compute_quadratic(z_full, q); 
    yQ = ct_tide_compute_quadratic(z_full, q_best); 
    patch([yQL' fliplr(yQU')], [z_full' fliplr(z_full')], 'k', 'Facealpha', 0.2, 'Edgecolor', 'none'); hold on
    plot(yQ, z_full, 'k', 'LineWidth', 1); hold on

    % Apply Linear Function over full extent
    yL = p_best(1)*z_full + p_best(2); 
    pU = p_best - p_best_se; yLU = pU(1)*z_full + pU(2); 
    pL = p_best + p_best_se; yLL = pL(1)*z_full + pL(2); 
    patch([yLL' fliplr(yLU')], [z_full' fliplr(z_full')], 'r', 'Facealpha', 0.2, 'Edgecolor', 'none'); hold on
    plot(yL, z_full, 'r', 'LineWidth', 1); hold on
    xlim([-2.5 6]); xlabel('vel (m a^{-1})');
    
    % Calculate and display Ftest values
    try
        y_data = v_all(ind_fit); 
        y_lin_fit = (p_best(1)*z_all(ind_fit) + p_best(2)); 
        y_quad_fit = ct_tide_compute_quadratic(z_all(ind_fit), q_best);
        w_fit = w(ind_fit);
        
        [dof, faj, ftab] = tg_rev_Ftest_weighted(y_data, y_lin_fit, y_quad_fit, w_fit);
        
        text(0.55, 0.94, sprintf('F = %.1f (Crit = %.1f)', faj, ftab), ...
             'Units', 'normalized', 'fontsize', 8, 'FontWeight', 'normal');
    catch
    end

    %% hax(2): Correlation Coefficient (Bottom Plot)
    set(gcf,'CurrentAxes',hax(2))
    
    plot(all_ampCor_mean', z_all, 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.5); hold on
    patch([mean_ampCor-std_ampCor fliplr(mean_ampCor+std_ampCor)], ...
          [z_all fliplr(z_all)], 'k', 'Facealpha', 0.2, 'Edgecolor', 'none'); hold on
    plot(mean_ampCor, z_all, 'k', 'LineWidth', 1.5); hold on
    
    xlim([0.75 1.01]); xlabel('Correlation coef');

    %% Apply final aesthetic formatting and custom limits/ticks
    abcd = 'ab';
    for j = 1:length(hax)
        set(gcf,'CurrentAxes',hax(j))
        set(gca, 'Fontsize', 8, 'Linewidth', 0.5, ...
                 'ytick', ylim_val(1):ytick_dist:ylim_val(2), ...
                 'ylim', ylim_val, 'ydir', 'reverse')
        if j > 1
            set(gca,'YTickLabel', [])
            ylabel('')
        else
            ylabel('Range (m)')
        end
        text(0.05, 0.96, [abcd(j)], 'Units', 'normalized', 'fontsize', 8, 'FontWeight', 'normal')
        axis tight
        set(gca, 'ylim', ylim_val, 'ytick', ylim_val(1):ytick_dist:ylim_val(2));
        try, ct_tide_plot_shade_drnf(xlim, drnf_artefacts, 'b'); catch, end
    end
    linkaxes(hax,'y')
    xlim(hax(1), [-2 1.5])
end

function [y] = ct_tide_compute_quadratic(x,b)
    y = b(1)+x.*b(2)+x.^2.*b(3); 
end

function [dof,faj,ftab] = tg_rev_Ftest_weighted(y_data,y_lin,y_quad,w)
    if ~isempty(find(w<0, 1))
        error('Cannot have negative weights')
    end
    w = w/max(w);
    y_data = reshape(y_data,length(y_data),1);
    y_lin = reshape(y_lin,length(y_lin),1);
    y_quad = reshape(y_quad,length(y_quad),1);
    w = reshape(w,length(w),1);
    SE_Q = sum(w.*(y_quad-mean(y_data)).^2); 
    SE_L = sum(w.*(y_lin-mean(y_data)).^2); 
    SU_Q = sum(w.*(y_quad-y_data).^2); 
    N = length(y_data);
    dof = sum(w);
    m = 1; 
    var_est = SU_Q/((N-m)*sum(w)/N);
    faj = (SE_Q-SE_L) / var_est;
    xx = -1:0.01:20;
    F1 = fcdf(xx,m,dof);
    [~, b_idx] = min(abs(F1-0.99));
    ftab = xx(b_idx);
end