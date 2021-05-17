%%% messy model oxygen model
%%% 5/13/21

%%% C' = (beta*x/(b + x))*c*(1 - (c + f)/k) - dc*c;
%%% F' = (r + beta*(1 - x/(b + x)))*f*(1 - (c + f)/k) - df*f  - (q*lambda*f)/(mu + eta*c);
%%% X' = lambda - mu*x - eta*x*c;

close all;

%%% =======================================================================

% fixed parameters
global k N0
 
N0 = 6.7e8;
k = 10^10;

%%% Ode parameters ========================================================

beta = 12.3;
r = 0;
d = 0.2;
q = 12.2e-10;

lambda = 9.7e7;
mu = 200*23*60*24;
eta = 2.2e-3;
b = 2.4;

p = [beta,r,d,q,...
     lambda,mu,eta,b];

frac = 0.85;
c0 = frac*N0;
f0 = (1 - frac)*N0;
x0 = 14;
y0 = [c0; f0; x0];

tspan = [0 180];
[t, y] = ode15s(@(t,y) cf_eqs(t,y,p), tspan, y0);


figure()
hold on; box on;
plot(t,log10(y(:,1)),'Linewidth',2)
plot(t,log10(y(:,2)),'Linewidth',2)
xlabel('Time (days)')
ylabel('Absolute Abundance')
title('Climax and Attack Populations')
legend('C model','F model')

figure()
plot(t,y(:,3),'Linewidth',2)
xlabel('Time (days)')
ylabel('Oxygen (\muM)')
title('Oxygen')

%%% equilibria ============================================================

% time series
c = y(:,1);
f = y(:,2);

% rf term
A = @(c)(beta*b*(mu + eta*c))/(lambda + b*(mu + eta*c));

% this is the same as the climax only ss
Jf0 = @(c) A(c)*(1-c/k) - d - q*lambda/(mu + eta*c);

% Jacobian of just F
Jf = @(c,f) A(c)*(1-c/k) - 2*A(c)*f/k - d - q*lambda/(mu + eta*c);

% Jf0(c(end)) < 0 & Jf(c(end),f(end)) < 0
Jf(c(end),f(end));

% coexistence 1
rad1 = sqrt(1 + 4*b*q/d);
rad2 = sqrt(d^2 + 4*d*b*q);

cx1 = (lambda + lambda*rad1 - 2*b*mu)/(2*b*eta);
fx1 = -((1 + rad1)*beta*lambda + b*(3*d*k*eta + k*eta*rad2 - 2*beta*(k*eta + mu)))/(2*b*beta*eta);

% coexistence 2
cx2 = (lambda - lambda*rad1 - 2*b*mu)/(2*b*eta);
fx2 = -((1 - rad1)*beta*lambda - b*(-3*d*k*eta + k*eta*rad2 + 2*beta*(k*eta + mu)))/(2*b*beta*eta);




%%% funcitons =============================================================

%%% cf ode function
function yp = cf_eqs(t,y,p)
global k

beta = p(1);
r = p(2);
d = p(3);
q = p(4);
lambda = p(5);
mu = p(6);
eta = p(7);
b = p(8);

dc = d;
df = d;

c = y(1);
f = y(2);
x = y(3);

yp = zeros(3,1);

yp(1) = (beta*x/(b + x))*c*(1 - (c + f)/k) - dc*c;
yp(2) = (r + beta*(1 - x/(b + x)))*f*(1 - (c + f)/k) - df*f  - (q*lambda*f)/(mu + eta*c);
yp(3) = lambda - mu*x - eta*x*c;

% hold on
% scatter(t, ((beta*x/(b + x))*(1 - (c + f)/k) - dc),'b');
% scatter(t, (r + beta*(1 - x/(b + x)))*(1 - (c + f)/k) - df  - (q*lambda)/(mu + eta*c),'rx');
end

%%% jacobian function
function J = jac(y,p)
global k

c = y(1);
f = y(2);
x = y(3);

beta = p(1);
r = p(2);
d = p(3);
q = p(4);
lambda = p(5);
mu = p(6);
eta = p(7);
b = p(8);

J = zeros(3,3);

J(1,1) = (beta*x/(b + x))*(1 - (c + f)/k) - d - beta*c*x/(k*(b+x));
J(1,2) = -beta*c*x/(k*(b+x));
J(1,3) = -(beta*x/(b + x)^2)*c*(1 + (c + f)/k) + beta*c*(1 - (c + f)/k)/(b + x);

% J(2,1) = 


end