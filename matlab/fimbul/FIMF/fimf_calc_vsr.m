function [ind_fit, q_best, q_se, p_best, p_best_se] = fimf_calc_vsr(z_all, v_all, v_all_se, drnf_vsr)
w = 1 ./ v_all_se;
w = w / max(w); 

[ind_fit, q_best, q_se] = ct_tide_quadratic_fit_local(z_all, v_all, w, drnf_vsr, 'lsq_weighted');
[~, p_best, ~, p_best_se] = ct_tide_line_fit_local(z_all, v_all, w, drnf_vsr, 'lsq_weighted');
end

function [ind,b,b_se] = ct_tide_quadratic_fit_local(x_in,y_in,w_in,drnf,opt_fit,startingVals)
if nargin < 6
    startingVals = [10e-3 10e-5 10e-2];  
end
x = x_in; y = y_in;
ind = find(isnan(x) == 0);
ind = intersect(ind,find(isnan(y) == 0));
for j = 1:size(drnf,1)
    if drnf(j,1) > drnf(j,2)
        disp(['Badly defined depth range, skipping ' num2str(drnf(j,1)) '-' num2str(drnf(j,2))])
    else
        ind = intersect(ind,find(x <= drnf(j,1) | x >= drnf(j,2)));
    end
end
ind = unique(ind);

% --- NEW: Return NaN if insufficient points for a 3-parameter fit ---
if length(ind) < 3
    b = [NaN; NaN; NaN];
    b_se = [NaN; NaN; NaN];
    return;
end

x = x(ind); x = reshape(x,length(x),1);
y = y(ind); y = reshape(y,length(y),1);
if ~isempty(w_in)
    w = w_in(ind); w = reshape(w,length(w),1);
end
modelFun = @(p,xx) p(1)+xx.*p(2)+xx.^2.*p(3); 
if strcmp(opt_fit,'robust')
    opts.RobustWgtFun = 'bisquare';
    [b,~,~,CovB,~,ErrorModelInfo] = nlinfit(x,y, modelFun, startingVals, opts);
elseif strcmp(opt_fit,'lsq_weighted')
    [b,~,~,CovB,~,ErrorModelInfo] = nlinfit(x,y, modelFun, startingVals, 'Weights',w);
end
b_se = sqrt(diag(CovB))'; 
end

function [ind,p,rmse,se] = ct_tide_line_fit_local(x_in,y_in,w_in,drnf,opt_fit)
x = x_in; y = y_in;
ind = find(isnan(x) == 0);
ind = intersect(ind,find(isnan(y) == 0));
for j = 1:size(drnf,1)
    if drnf(j,1) > drnf(j,2)
        disp(['Badly defined depth range, skipping ' num2str(drnf(j,1)) '-' num2str(drnf(j,2))])
    else
        ind = intersect(ind,find(x <= drnf(j,1) | x >= drnf(j,2)));
    end
end
ind = unique(ind);

% --- NEW: Return NaN if insufficient points for a 2-parameter fit ---
if length(ind) < 2
    p = [NaN; NaN];
    se = [NaN; NaN];
    rmse = NaN;
    return;
end

x = x(ind); x = reshape(x,length(x),1);
y = y(ind); y = reshape(y,length(y),1);
if ~isempty(w_in)
    w = w_in(ind); 
end
if strcmp(opt_fit,'robust')
    [p,stats] = robustfit(x,y,'bisquare');
    p = flipud(p);
    rmse = stats.robust_s;
    se = flipud(stats.se);
elseif strcmp(opt_fit,'lsq_weighted')
    x = reshape(x,1,length(x)); 
    y = reshape(y,1,length(y));
    A = [x' ones(length(x),1)];
    if ~isempty(w_in)
        w = reshape(w,1,length(w));
        [p,se,mse] = lscov(A,y',w');
    else
        [p,se,mse] = lscov(A,y');
    end
    rmse = sqrt(mse);
end
end