%%% using this to look at transient properites of the best fitting
%%% parameters and look at the levels for the coexistence equilibrium
%%% model is

% C' = (beta*x^n/(b^n + x^n))*c*(1 - (c + f)/k) - dc*c;
% F' = (r + beta*(1 - x^n/(b^n + x^n)))*f*(1 - (f + c)/k) - df*f - ep*f - q*f*x;
% X' = lambda - mu*x - eta*(c)*x;


%%% 4/26/2020

% close all;

data = xlsread('C:\Users\peter\OneDrive\Desktop\cyst fib\julia stuff\ODEs\Data fitting\cf data','Rescaled');
tdata = data(:,1);
% cdata = data(:,2)/100;
% fdata = data(:,3)/100;
cdata = data(:,2);
fdata = data(:,3);

%%% =======================================================================

% fixed parameters
global k lambda t_b t_c N0 mu

% global parameters for treatment
global t_start t_end treat_true
 
N0 = 6.7e8;
t_b = 19*0;
t_c = 33*Inf;

%%% =======================================================================
%%% initial oxygen
x0 = 14.6287;

% parameters
r = 0.0046;

beta = 16.6388;
b = 13.4256;
n = 2.6626;

dn = 0.6045; % natural death rate
dbs = 6.7686; % death due to bs antibiotics
gamma = 0.8976; % fractional reduction of bs antibiotics in killing attack

ep = 1.2124;
mu = 200*23*60*24; % 1/5 min

k = 10^10;
eta = 3.1611e-4; 
q = 3.2747e-5;

frac = 0.8659;

lambda = mu*x0;
% lambda = 11.6901e+07;

p = [x0,frac,beta,r,...
     eta,dbs,dn,gamma,...
     ep,q,b,n];


%%% Treament simulation stuff in here =====================================
t_start = Inf;
t_end = Inf;
treat_true = 0;


%%% =======================================================================

% solve ode's
x0 = p(1);
frac = p(2);
c0 = frac*N0;
% c0 = 1e-2*0;
f0 = (1 - frac)*N0;

y0 = [c0; f0; x0];
tspan = [0:0.01:100];
[t, y] = ode15s(@(t,y) cf_eqs(t,y,p), tspan, y0);

%%% relative abundances
Ct = y(:,1)./(y(:,1) + y(:,2));
Ft = y(:,2)./(y(:,1) + y(:,2));

%%% this stuff will find the time between exacerbations ===================
%%% =======================================================================

tol = 1e-2;

% swtchpts = find(abs(Ft - Ct) < tol); % find point where they switch
% swtimes = t(swtchpts) % time when they switch

try
    swtchpts = find(Ft > Ct);
    swtimes = t(swtchpts(1))
catch
    % do nothing if error
end



%%% =======================================================================

figure()
hold on; box on;
plot(t,Ct,'Linewidth',2)
plot(t,Ft,'Linewidth',2)
% plot(tdata,cdata,'bx', 'LineWidth',2)
% plot(tdata,fdata,'rx', 'LineWidth',2)
xlabel('Time (days)')
ylabel('Relative Abundance')
title('Climax and Attack Populations')
xline(swtimes)
% xline(t_c)
% legend('C model','F model','C data','F data','Location','e')
legend('C model','F model','Population switch','Location','e')

figure()
hold on; box on;
plot(t,log10(y(:,1)),'Linewidth',2)
plot(t,log10(y(:,2)),'Linewidth',2)
xlabel('Time (days)')
ylabel('Absolute Abundance')
title('Climax and Attack Populations')
legend('C model','F model','Location','e')

% figure()
% hold on; box on;
% plot(t,y(:,3),'Linewidth',2)
% xlabel('Time (days)')
% ylabel('Oxygen (\muM)')
% title('Oxygen')

%%% Functions =============================================================

%%% broad spectrum antibiotic function
function dbs = BrSpec(t,p)
    global t_b
    if t < t_b
        dbs = p(6);
    else 
        dbs = 0;
    end
end

%%% cf ode function
function yp = cf_eqs(t,y,p)
global k lambda t_c mu 
global t_start t_end treat_true

beta = p(3);
r = p(4);
eta = p(5);
dbs = BrSpec(t,p);
dn = p(7);
gamma = p(8);
ep = 0;
q = p(10);
b = p(11);
n = p(12);

%%% total death rates
dc = dn + dbs;
df = dn + gamma*dbs;
% [t dc df]

c = y(1);
f = y(2);
x = y(3);

yp = zeros(3,1);

yp(1) = (beta*x^n/(b^n + x^n))*c*(1 - (c + f)/k) - dc*c;
yp(2) = (r + beta*(1 - x^n/(b^n + x^n)))*f*(1 - (f + c)/k) - df*f - ep*f - q*f*x;
yp(3) = lambda - mu*x - eta*(c)*x;

% hold on
% scatter(t,(beta*x^n/(b^n + x^n)),'bo')
% scatter(t,(r + beta*(1 - x^n/(b^n + x^n))),'rx')
% scatter(t,(beta*x^n/(b^n + x^n))*(1 - (c + f)/k) - dc,'bo')
% scatter(t,(r + beta*(1 - x^n/(b^n + x^n)))*(1 - (f + c)/k) - df - ep - q*x,'rx')
end

%%% objective function for cf_fitter
function J = cf_err(p,tdata,cdata,fdata)
global N0 

x0 = p(1);
frac = p(2);

c0 = frac*N0;
f0 = (1 - frac)*N0;

y0 = [c0; f0; x0];
[t,y] = ode15s(@cf_eqs,tdata,y0,[],p);

Ct = y(:,1)./(y(:,1) + y(:,2));
Ft = y(:,2)./(y(:,1) + y(:,2));

errx = Ct - cdata;
erry = Ft - fdata;

J = errx'*errx + erry'*erry;
% J = errx'*errx;
% J = erry'*erry;
end

%%% Function to return two error vectors (and ode solution in rel. abund.)
function [sol,C_err,F_err] = err_vec(p,tdata,cdata,fdata)
global N0
% solve ode
x0 = p(1);
frac = p(2);

c0 = frac*N0;
f0 = (1 - frac)*N0;

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