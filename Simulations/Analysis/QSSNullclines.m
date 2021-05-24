%%% plot nullclines of simple scaled quasi-steady state oxygen model 
%%% also get numerical value of coexistence steady state and eigenvalues
%%% in this, r_c(X) = beta*X/(b + X)

%%% 5/21/21

%%% c' = (beta*lambda/((lambda+b*mu) + b*eta*k*c))*c*(1 - c - f) - dc*c
%%% f' = r*f*(1 - c - f) - df*f - q*lambda*f/(mu + eta*k*c)

%%% this probably has a type in either the jacobian function or c and f
%%% excclusion terms

close all;

%%% Parameters ============================================================
global k beta r d b mu eta lambda q 
k = 10^10;

beta = 16.0;
r = 20.0;
d = 0.2;
b = 12.4;

mu = 200*23*60*24;
eta = 2.2e-1;

lambda = 9.9e7;
q = 5e-1;

%%% nullcline curves
Cp = @(c,f) (beta*lambda./((lambda + b*mu) + b*eta*k*c)).*(1-c-f) - d;
Fp = @(c,f) r*(1-c-f) - d - q*lambda./(mu + k*eta*c);

%%% ODE solver ============================================================
%%% ode stuff
c0 = 0.08399;
f0 = 0.01333;

y0 = [c0; f0];
tspan = [0 800];

[t, y] = ode15s(@(t,y) cf_eqs(t,y,p), tspan, y0);

% time series
c = y(:,1);
f = y(:,2);



%%% find steady states numerically ========================================

%%% numerical coexistence value
x0 = y0;
fun = @SStates;
xc = fsolve(fun,x0);


%%% coexistence eigenvalues
Jx = my_jac(xc);
vx = eigs(Jx);

%%% exclusion equilibria
conly = -(d*lambda - beta*lambda + b*d*mu)/(b*d*k*eta + beta*lambda);
fonly = (-q*lambda - d*mu + r*mu)/(r*mu);

%%% exclusion eigenvalues
Jc = my_jac([conly,0]);
vc = eigs(Jc);

Jf = my_jac([0,fonly]);
vf = eigs(Jf);

%%% extinction
Jxx = my_jac([0,0]);
vxx = eigs(Jxx);

% [vxx vx vc vf]
[eigs(my_jac([c(end),f(end)])) vx]




%%% plots =================================================================

interval = [0.,1.];

%%% plot nullclines
hold on
% xline(0,'Linewidth',2)
% yline(0,'Linewidth',2)
fimplicit(Cp,interval,'b','Linewidth',2)
fimplicit(Fp,interval,'r','Linewidth',2)
% scatter(conly,0,'b','Linewidth',2)
% scatter(0,fonly,'r','Linewidth',2)
% scatter(xc(1),xc(2),'kx','Linewidth',4)
xlabel('C')
ylabel('F')

%%% ode lines
plot(c,f,'-.','Linewidth',2)
scatter(c(end),f(end),'kx','Linewidth',4)



%%% functions =============================================================

%%% cf ode function
function yp = cf_eqs(t,y,p)
global k beta r d b mu eta lambda q 

c = y(1);
f = y(2);

yp = zeros(2,1);

yp(1) = (beta*lambda/((lambda+b*mu) + b*eta*k*c))*c*(1 - c - f) - d*c;
yp(2) = r*f*(1 - c - f) - d*f - q*lambda*f/(mu + eta*k*c);

% hold on
% scatter(t,(beta*lambda/(mu + eta*c))*(1 - c - f) - dc,'b')
% scatter(t,r*(k*f*(mu + eta*c)/(q*lambda) - 1)*(1 - c - f) - df,'r')
end

%%% coexistence residual function
function F = SStates(x)
global k beta r d b mu eta lambda q 

F(1) = (beta*lambda./((lambda + b*mu) + b*eta*k*x(1))).*(1 - x(1) - x(2)) - d;
F(2) = r*(1 - x(1) - x(2)) - d - q*lambda./(mu + k*eta*x(1));

end

function J = my_jac(x)
global k beta r d b mu eta lambda q 

c = x(1);
f = x(2);

J = zeros(2,2);

bottom = b*c*k*eta + lambda + b*mu;

J(1,1) = -b*(1-c-f)*k*beta*eta*lambda/bottom^2 - b*lambda/bottom;
J(1,2) = -beta*lambda/bottom;
J(2,1) = -r + k*q*eta*lambda/((c*k*eta+mu)^2);
J(2,2) = -r;

end

