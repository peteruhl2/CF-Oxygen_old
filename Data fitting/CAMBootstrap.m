%%% Program to compute 95% confidence intervals for the CF ode model
%%% by bootstrapping the erros
%%% 3/10/21

close all;

data = xlsread('C:\Users\peter\OneDrive\Desktop\cyst fib\julia stuff\ODEs\Data fitting\cf data','Rescaled');
tdata = data(:,1);
cdata = data(:,2);
fdata = data(:,3);

%%% =======================================================================

% fixed parameters
global beta r d k lambda mu eta t_treat N0

N0 = 6.7e8;
lambda = 1.32;
t_treat = 28.;

%%% set up stuff for the bootstrap loop ===================================

runs = 100;
tspan = 0 : 0.5 : 40;
steps = length(tspan);
C_results = zeros(steps,runs);
F_results = zeros(steps,runs);

% set results for number of parameters
n_params = 5;
p_results = zeros(n_params,runs);

%%% best fitting parameters
r = 30.5;

beta = 0.95*r; % try < 16
b = 2.9440e-05;
n = 1.7652;

d = 6.0;

ep = 0.5642;
mu = 1.32; % 1/5 min

k = 10^12;
eta = 2.7e-9;
q = 31.0727;

frac = 0.8459;
c0 = frac*N0;
f0 = (1 - frac)*N0;
x0 = 1.32;


% tspan = [0 : 0.5 : 40];
% [t, y] = ode15s(@(t,y) cf_eqs(t,y,p), tspan, y0);
% 
% 
% 
% [sol,C_err,F_err] = err_vec(p,tdata,cdata,fdata)


tic
for i = 1:runs
    i
    p = [b,n,q,ep,frac];

    % get ode solution and error vectors
    [sol,C_err,F_err] = err_vec(p,tdata,cdata,fdata);

    % attach time points to errors
    all_errs = [tdata C_err F_err];

    % call bootstrap function
    new_errs = Bootstrap(all_errs);

    % add bootstraped errors to sol and get new data
    new_sol = add_error(sol, new_errs);
    new_cdata = new_sol(:,2);
    new_fdata = new_sol(:,3);

    %%% =======================================================================

    try
        tic
        [p,fval,flag,output] = fminsearch(@cf_err,p,options,tdata,new_cdata,new_fdata);
        toc
    catch
        % do nothing if there's an error
    end

    % solve ode's
    frac = p(5);
    c0 = frac*N0;
    f0 = (1 - frac)*N0;
    x0 = 1.32;
    y0 = [c0; f0; x0];
    
    % tspan = [0 : 0.5 : 40];
    [t, y] = ode15s(@(t,y) cf_eqs(t,y,p), tspan, y0);
    % J = cf_err(p,tdata,cdata,fdata)
    p

    %%% relative abundances
    Ct = y(:,1)./(y(:,1) + y(:,2));
    Ft = y(:,2)./(y(:,1) + y(:,2));

    %%% store results
    C_results(:,i) = Ct;
    F_results(:,i) = Ft;
    p_results(:,i) = p;

end % bootstrap loop
toc

% get percentiles =========================================================
C_mean = zeros(length(tspan),1);
C_25 = zeros(length(tspan),1);
C_975 = zeros(length(tspan),1);

F_mean = zeros(length(tspan),1);
F_25 = zeros(length(tspan),1);
F_975 = zeros(length(tspan),1);


%%% fix the first dimension 1/19/21
p_mean = zeros(n_params,1);
p_25 = zeros(n_params,1);
p_975 = zeros(n_params,1);

%%% compute percentiles
for i = 1:length(tspan)
    C_mean(i) = mean(C_results(i,:));
    C_25(i) = prctile(C_results(i,:),2.5);
    C_975(i) = prctile(C_results(i,:),97.5);
    
    F_mean(i) = mean(F_results(i,:));
    F_25(i) = prctile(F_results(i,:),2.5);
    F_975(i) = prctile(F_results(i,:),97.5);
    
end

%%% get percentiles for p
for i = 1:n_params
    p_mean(i) = mean(p_results(i,:));
    p_25(i) = prctile(p_results(i,:),2.5);
    p_975(i) = prctile(p_results(i,:),97.5);
end


%%% =======================================================================
%%% Plots =================================================================
%%% relative abundances
Ct = y(:,1)./(y(:,1) + y(:,2));
Ft = y(:,2)./(y(:,1) + y(:,2));


hold on; box on;
plot(t,Ct,'Linewidth',2)
plot(t,Ft,'Linewidth',2)
plot(tdata,cdata,'bx', 'LineWidth',2)
plot(tdata,fdata,'rx', 'LineWidth',2)
xlabel('Time (days)')
ylabel('Population density')
title('Climax and Attack Populations')
legend('C model','F model','C data','F data','Location','w')


%%% =======================================================================

%%% Functions =============================================================

%%% cf ode function
function yp = cf_eqs(t,y,p)
global beta r d k lambda mu eta t_treat

b = p(1); 
n = p(2);
q = p(3);
ep = 0; 


if t >= t_treat
    ep = p(4);
end

c = y(1);
f = y(2);
x = y(3);

yp = zeros(3,1);

yp(1) = (beta*x^n/(b^n + x^n))*c*(1 - (c + f)/k) - d*c;
yp(2) = r*f*(1 - (f + c)/k) - d*f - ep*f - q*f*x;
yp(3) = lambda - mu*x - eta*c*x;

end

%%% objective function for cf_fitter
%%% this one returns SSE
%%% objective function for cf_fitter
function J = cf_err(p,tdata,cdata,fdata)
global N0

frac = p(5);

c0 = frac*N0;
f0 = (1 - frac)*N0;
x0 = 1.32;

y0 = [c0; f0; x0];
[t,y] = ode15s(@cf_eqs,tdata,y0,[],p);

Ct = y(:,1)./(y(:,1) + y(:,2));
Ft = y(:,2)./(y(:,1) + y(:,2));

errx = Ct - cdata;
erry = Ft - fdata;

J = errx'*errx + erry'*erry;
end

%%% Function to return two error vectors (and ode solution in rel. abund.)
function [sol,C_err,F_err] = err_vec(p,tdata,cdata,fdata)
global N0
% solve ode
frac = p(5);

c0 = frac*N0;
f0 = (1 - frac)*N0;
x0 = 1.32;
y0 = [c0; f0; x0];
[t,y] = ode15s(@cf_eqs,tdata,y0,[],p);

% return error vectors
Ct = y(:,1)./(y(:,1) + y(:,2));
Ft = y(:,2)./(y(:,1) + y(:,2));

C_err = Ct - cdata;
F_err = Ft - fdata;

% return solution too
sol = [t y];
sol(:,2) = Ct;
sol(:,3) = Ft;

end

%%% function to bootstrap the error vector
%%% errs is a three column matrix with the time data
%%% and two vectors of errors
function new_err = Bootstrap(errs)
N = length(errs);

% sample with replacement and sort
new_err = datasample(errs,N);
new_err = sortrows(new_err);

end


%%% function to add bootstraped error to ode solution
function new_sol = add_error(sol, err)
N = length(err);

for i = 1:N
    % find index of time point in sol
    sol_ind = find(sol(:,1) == err(i,1));
    % climax error
    sol(sol_ind,2) = sol(sol_ind,2) + err(i,2);
    % attack error
    sol(sol_ind,3) = sol(sol_ind,3) + err(i,3);    
end

new_sol = sol;

end