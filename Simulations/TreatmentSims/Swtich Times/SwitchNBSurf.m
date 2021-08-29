%%% find the switch times for different ranges of parameters
%%% this one a surface for n and b
%%% to get to 0.5 relative abundance

% C' = (beta*x^n/(b^n + x^n))*c*(1 - (c + f)/k) - dc*c;
% F' = (r + beta*(1 - x^n/(b^n + x^n)))*f*(1 - (f + c)/k) - df*f - ep*f - q*f*x;
% X' = lambda - mu*x - eta*(c)*x;


%%% 6/15/2020

close all;

%%% =======================================================================

% fixed parameters
global k lambda t_b t_c N0 mu

% global parameters for treatment
global t_start t_end treat_true
 
N0 = 6.7e8;
t_b = 19*0;
t_c = 33*Inf;

%%% =======================================================================
res = 400;

N = linspace(0,6,res);
B = linspace(12,20,res);
results = zeros(length(N),length(B));

%%% parameters
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

for i = 1:length(N)
    for j = 1:length(B)
        [i j] 
        
        n = N(i);
        b = B(j);

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

        tol = 5e-1;

    %     swtchpts = find(abs(Ft - Ct) < tol); % find point where they switch

        try
            swtchpts = find(Ft>Ct);
            swtimes = t(swtchpts); % time when they switch

            results(i,j) = swtimes(1);
        catch
%             results(i,j) = NaN;
            %%% keep track of what happens to f
            %%% if f is extinct results = NaN
            %%% if f survives but doesnt cross restuls = 0
            if Ft(end) > 1e-3
                results(i,j) = NaN;
            elseif Ft(end) < 1e-3
                results(i,j) = NaN;
            end
        end
    end
end

%%% =======================================================================

[X,Y] = meshgrid(N,B);

hold on; box on;
% surf(X,Y,Fvals')
contourf(X,Y,results')
xlabel('Slope factor - n')
ylabel('Half-saturation value - b')
h2 = colorbar;
ylabel(h2, 'Days to population switch')
shading interp
caxis([0, 60])
title('Time to population switch, n and b')



figure()
hold on; box on;
plot(t,Ct,'Linewidth',2)
plot(t,Ft,'Linewidth',2)
% plot(tdata,cdata,'bx', 'LineWidth',2)
% plot(tdata,fdata,'rx', 'LineWidth',2)
xlabel('Time (days)')
ylabel('Relative Abundance')
title('Climax and Attack Populations')
% xline(t_b)
% xline(t_c)
% legend('C model','F model','C data','F data','Location','e')
legend('C model','F model','Location','e')

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