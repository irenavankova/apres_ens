function [t, dhdtmat, f1] = fimf_load_and_create_timeseries(sitename, f0, fend, db, bw, i1t, opt_bed, dhRange_low, dhRange_high)
%FIMF_LOAD_AND_CREATE_TIMESERIES Load .mat files and compute dh/dt timeseries.
%   [t, dhdtmat, f1] = FIMF_LOAD_AND_CREATE_TIMESERIES(...)

f1 = f0:db:fend;
N = length(f1);
t = [];
dhdtmat = [];

for k = 1:N
    sname = [sitename num2str(f1(k)) num2str(f1(k)+bw)];

    if opt_bed
        data = load([sname '_bed_tpq.mat']);
        if isempty(t)
            t = data.tp.tim(i1t:end);
            t = t(:); % Force column vector
        end
        dh = data.tp.thickness(i1t:end);
        dh = dh(:); % Force column vector
    else
        data = load([sname '_ts_fine_TT.mat']);
        [~, ind] = min(abs(dhRange_low - data.ct.dts_uwrp.dhRange));
        [~, ind2] = min(abs(dhRange_high - data.ct.dts_uwrp.dhRange));
        dh = mean(data.ct.dts_xcor.dh(i1t:end, ind:ind2), 2);
        dh = dh(:); % Force column vector
        if isempty(t)
            t = data.ct.dts_xcor.time(i1t:end);
            t = t(:); % Force column vector (remove original transpose)
        end
    end

    % Compute dh/dt (now guaranteed to be a column vector)
    dhdt = diff(dh) ./ diff(t);

    % Initialize dhdtmat on first iteration
    if isempty(dhdtmat)
        dhdtmat = zeros(length(dhdt), N);
    end

    % Assign to matrix
    dhdtmat(:, k) = dhdt;
end
end