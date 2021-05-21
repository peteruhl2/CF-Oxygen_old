%%% Allee effect and quasi-steady state oxygen, scaled
%%% 5/18/21

%%% c' = (beta*lambda/(mu + eta*c))*c*(1 - c - f) - dc*c
%%% f' = r*f*(f*(mu + eta*c)/(q*lambda) - 1)*(1 - c - f) - df*f

%%% =======================================================================
%%% parameters
global k
k = 10^10;

beta = 18.1;
r = 1.5;
d = 0.8;
q = k*1e-3;
lambda = 9.7e7;
mu = 200*23*60*24;
eta = 2.2e-4;

p = [beta, r, d, q,...
     lambda, mu, eta];

%%% ode stuff
c0 = 0.6;
f0 = 0.5;

y0 = [c0; f0];
tspan = [0 80];

[t, y] = ode15s(@(t,y) cf_eqs(t,y,p), tspan, y0);

% time series
c = y(:,1);
f = y(:,2);

% %%% phase plane
% hold on; box on
% plot(c,f,'Linewidth',2)
% scatter(c(end),f(end))
% xlabel('C')
% ylabel('F')

% time series plot
figure()
hold on; box on;
plot(t,(y(:,1)),'Linewidth',2)
plot(t,(y(:,2)),'Linewidth',2)
xlabel('Time (days)')
ylabel('Absolute Abundance')
title('Climax and Attack Populations')
legend('C model','F model')






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

dc = d;
df = d;

c = y(1);
f = y(2);

yp = zeros(2,1);

yp(1) = (beta*lambda/(mu + eta*c))*c*(1 - c - f) - dc*c;
yp(2) = r*(k*f*(mu + eta*c)/(q*lambda) - 1)*f*(1 - c - f) - df*f;

hold on
scatter(t,(beta*lambda/(mu + eta*c))*(1 - c - f) - dc,'b')
scatter(t,r*(k*f*(mu + eta*c)/(q*lambda) - 1)*(1 - c - f) - df,'r')
end