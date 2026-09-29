function [v_int, v_bed] = fimf_calc_dhdt(dts, dts_bed, tind, opt_vel_method, opt_deriv_combine, opt_bed_source)
d2y = 365.25; 

time = dts.time(tind);
depth = dts.dhRange;
v_int = nan(1, length(depth));

for j = 1:length(depth)
    y = dts.dh(tind, j);
    valid = ~isnan(y) & ~isnan(time');

    if sum(valid) > 3
        v_int(j) = calc_single_vel(time(valid), y(valid), opt_vel_method, opt_deriv_combine, d2y);
    end
end

v_bed = nan;
if ~isempty(dts_bed)
    if strcmp(opt_bed_source, 'tp')
        time_bed = dts_bed.tim(tind);
        y_bed = dts_bed.thickness(tind);
    else
        time_bed = dts_bed.time(tind);
        y_bed = dts_bed.dh(tind);
    end

    valid = ~isnan(y_bed) & ~isnan(time_bed');

    if sum(valid) > 3
        v_bed = calc_single_vel(time_bed(valid), y_bed(valid), opt_vel_method, opt_deriv_combine, d2y);
    end
end
end

function v = calc_single_vel(t, y, method, combine_opt, d2y)
t = t(:); y = y(:);

switch method
    case 'deriv_outlier'
        dy = gradient(y) ./ gradient(t); 
        ind = isoutlier(dy, 'median'); 
        dy(ind) = [];
        v = mean(dy, 'omitnan') * d2y;

    case 'deriv_simple'
        dy = gradient(y) ./ gradient(t);
        if strcmp(combine_opt, 'median')
            v = median(dy, 'omitnan') * d2y;
        else
            v = mean(dy, 'omitnan') * d2y;
        end

    case 'line_fit'
        warning('off', 'stats:statrobustfit:IterationLimit');
        p = robustfit(t, y, 'bisquare');
        warning('on', 'stats:statrobustfit:IterationLimit');
        v = p(2) * d2y; 

    otherwise
        error('Unknown velocity calculation method specified.');
end
end