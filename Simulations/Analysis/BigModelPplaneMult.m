%%% plot nullclines of the model that the data is fit to, plus numerical
%%% steady states and eigenvalues
%%% in this r_c(X) = beta*x^n/(b^n + x^n)
%%% model is

%%% c' = rc*c*(1-c-f) - d*c
%%% f' = rf*f*(1-c-f) - d*f - q*lambda*f/(mu + eta*k*c)

%%% rc = (beta*lambda^n)/(lambda^n + b^n*(mu + eta*k*c)^n)
%%% rf = (r + beta*(1 - (lambda^n)/(lambda^n + b^n*(mu + eta*k*c)^n)))

%%% this one is multiple trajectories
%%% 5/23/21

close all;

%%% Parameters ============================================================
global k beta r d b mu eta lambda q n
k = 10^10;

beta = 16.6;
r = 0.4;
d = 0.6;
b = 13.4;
n = 2.6;

mu = 200*23*60*24;
eta = 3.1e-3;

lambda = 9.2e7;
q = 3e-1;

%%% =======================================================================
%%% growth rates
rc = @(c) (beta*lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n);
rf = @(c) (r + beta*(1 - (lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n)));

%%% nullcline functions
Cp = @(c,f) (beta*lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n).*(1-c-f) - d;
Fp = @(c,f) (r + beta*(1 - (lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n))).*(1-c-f) - d - q*lambda./(mu + eta*k*c);

%%% ODE solver ============================================================
%%% ode stuff
c0 = 0.1;
f0 = 0.001;

y0 = [c0; f0];
tspan = [0 800];

[t, y] = ode15s(@(t,y) cf_eqs(t,y,p), tspan, y0);





%%% plots =================================================================
%%% lotta runs here =======================================================

runs = 50;

for i = 1:runs
    i
    c0 = unifrnd(0,1);
    f0 = unifrnd(0,1);
    y0 = [c0; f0];

    tspan = [0 600];
    [t, y] = ode15s(@(t,y) cf_eqs(t,y,p), tspan, y0);

    hold on; box on;
    plot((y(:,1)),(y(:,2)),'-.','Linewidth',1.5)
    scatter((y(end,1)),(y(end,2)),'Linewidth',2)
    xlabel('C')
    ylabel('F')
end
%%% =======================================================================

% time series
c = y(:,1);
f = y(:,2);

interval = [0,1];
hold on
fimplicit(Cp, interval,'b','Linewidth',2)
fimplicit(Fp, interval,'r','Linewidth',2)
xlabel('C')
ylabel('F')

%%% ode lines
plot(c,f,'-.','Linewidth',2)
scatter(c(end),f(end),'kx','Linewidth',4)

xlim([0 1])
ylim([0 1])



%%% functions =============================================================

%%% cf ode function
function yp = cf_eqs(t,y,p)
global k beta r d b mu eta lambda q n

c = y(1);
f = y(2);

yp = zeros(2,1);

yp(1) = (beta*lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n).*c*(1-c-f) - d*c;
yp(2) = (r + beta*(1 - (lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n))).*f*(1-c-f) - d*f - q*lambda*f./(mu + eta*k*c);

% hold on
% scatter(t,(beta*lambda/(mu + eta*c))*(1 - c - f) - dc,'b')
% scatter(t,r*(k*f*(mu + eta*c)/(q*lambda) - 1)*(1 - c - f) - df,'r')
end

%%% steady state residual function
function F = SStates(x)
global k beta r d b mu eta lambda q n

c = x(1);
f = x(2);

F(1) = (beta*lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n).*c*(1-c-f) - d*c;
F(2) = (r + beta*(1 - (lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n))).*f*(1-c-f) - d*f - q*lambda*f./(mu + eta*k*c);

end
