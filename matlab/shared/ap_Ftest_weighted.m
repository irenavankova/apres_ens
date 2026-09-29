function [dof,faj,ftab] = ap_Ftest_weighted(y_data,y_lin,y_quad,w)
% faj is the F test value as calculated by Jenkins 2006
% ftab is the table value of F given the model parameters (m) and degrees of
% greedom (dof)

%y_data = vused;
%y_lin = p(1)*dused + p(2);
%y_quad = yquad;

if ~isempty(find(w<0))
    error('Cannot have negative weights')
end

w = w/max(w);

y_data = reshape(y_data,length(y_data),1);
y_lin = reshape(y_lin,length(y_lin),1);
y_quad = reshape(y_quad,length(y_quad),1);
w = reshape(w,length(w),1);

SE_Q = sum(w.*(y_quad-mean(y_data)).^2); % explained variance of quadratic
SE_L = sum(w.*(y_lin-mean(y_data)).^2); % explained variance of linear
SU_Q = sum(w.*(y_quad-y_data).^2); % unexplained variance of quadratic
%dof = length(y_data);
N = length(y_data);
dof = sum(w);
m = 1; %model paramters
%var_est = SU_Q/(dof-m);
var_est = SU_Q/((N-m)*sum(w)/N);

faj = (SE_Q-SE_L) / var_est;

xx = -1:0.01:20;
F1 = fcdf(xx,m,dof);
[a, b] = min(abs(F1-0.99));
%xx(b)

% F12 = fcdf(xx,m,N);
% [a2, b2] = min(abs(F12-0.99));
% xx(b2)

fconf = F1(b);

ftab = xx(b);
