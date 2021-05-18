%%% reduced constant oxygen model with Allee effect for F
%%% 5/18/21

%%% c' = beta*c*(1 - c - f) - dc*c
%%% f' = r*f*(k*f/q - 1)*(1 - c - f) - df*f

%%% looks like no bistability here

%%% =======================================================================
%%% parameters
global k
k = 1;

beta = 12.1;
r = 6.5;
d = 9.2;
q = k*0.1;

p = [beta, r, d, q];

%%% ode stuff
c0 = 0.7;
f0 = 0.8;

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


%%% Equilibria ============================================================

%%% exclusions
% c only
cs = (beta - d)/beta;

% f only
D1 = sqrt(r^2 - 4*d*q*r - 2*q*r^2 + (q^2)*(r^2));
fs1 = (r + q*r - D1)/(2*r);
fs2 = (r + q*r + D1)/(2*r);

% coexistence
cx = 1 - d/beta - q*(r+beta)/r;
fx = q*(r+beta)/r;














%%% funcitons =============================================================

%%% cf ode function
function yp = cf_eqs(t,y,p)
global k

beta = p(1);
r = p(2);
d = p(3);
q = p(4);

dc = d;
df = d;

c = y(1);
f = y(2);

yp = zeros(2,1);

yp(1) = beta*c*(1 - c - f) - dc*c;
yp(2) = r*f*(k*f/q - 1)*(1 - c - f) - df*f;
end