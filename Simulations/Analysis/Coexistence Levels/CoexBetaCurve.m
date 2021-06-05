%%% plot nullclines of the model that the data is fit to, draw out the path
%%% of the coexistence ss in the phase plane as a parameter varies
%%% this one is for beta

%%% in this r_c(X) = beta*x^n/(b^n + x^n)
%%% model is

%%% c' = rc*c*(1-c-f) - d*c
%%% f' = rf*f*(1-c-f) - d*f - q*lambda*f/(mu + eta*k*c)

%%% rc = (beta*lambda^n)/(lambda^n + b^n*(mu + eta*k*c)^n)
%%% rf = (r + beta*(1 - (lambda^n)/(lambda^n + b^n*(mu + eta*k*c)^n)))

%%% this one is multiple trajectories
%%% 5/23/21

close all;
global k beta r d b mu eta lambda q n

%%% array for results =====================================================
Beta = linspace(0,30,100);
Epts = zeros(length(Beta),2);
V = zeros(length(Beta), 2);

options = optimset('Display','off');
fun = @SStates;

%%% Parameters ============================================================

k = 10^10;

beta = 16.6;
r = 0.004;
d = 0.6;
b = 13.4;
n = 2.6;

mu = 200*23*60*24;
eta = 3.1e-3;

lambda = 9.6e7;
q = 3e-1;

%%% =======================================================================
%%% growth rates
rc = @(c) (beta*lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n);
rf = @(c) (r + beta*(1 - (lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n)));

%%% nullcline functions
Cp = @(c,f) (beta*lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n).*(1-c-f) - d;
Fp = @(c,f) (r + beta*(1 - (lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n))).*(1-c-f) - d - q*lambda./(mu + eta*k*c);


hold on
for i = 1:length(Beta)
    %%% new value of q
    beta = Beta(i);
    
    %%% have to bring the nullcline functions in here
    %%% growth rates
    rc = @(c) (beta*lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n);
    rf = @(c) (r + beta*(1 - (lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n)));

    %%% nullcline functions
    Cp = @(c,f) (beta*lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n).*(1-c-f) - d;
    Fp = @(c,f) (r + beta*(1 - (lambda^n)./(lambda^n + (b^n)*(mu + eta*k*c).^n))).*(1-c-f) - d - q*lambda./(mu + eta*k*c);
    
    % calculate nullcline intersection and eigenvalues
    x0 = [0.5;0.5];
    xco = fsolve(fun,x0,options);
    Epts(i,:) = xco;
    V(i,:) = cfeigs(xco);
    
    fimplicit(Cp, interval,'b','Linewidth',2)
    fimplicit(Fp, interval,'r','Linewidth',2)
    
    if max(V(i,:)) > 0 % unstable
        scatter(xco(1),xco(2),'kx','Linewidth',4)
    else 
        scatter(xco(1),xco(2),'cx','Linewidth',4)
    end
end

%%% Coesistence point =====================================================
x0 = [0.5;0.5];
xco = fsolve(fun,x0,options);
xf = fsolve(fun,[0,1],options);
xc = fsolve(fun,[1,0],options);
x0 = fsolve(fun,[0.01,0.01],options);

%%% eigenvalues
v = cfeigs(xco);
v0 = cfeigs([0,0]);
vc = cfeigs(xc);
vf = cfeigs(xf);

%%% trajectory
c0 = 0.9;
f0 = 0.9;

y0 = [c0; f0];
tspan = [0 800];

[t, y] = ode15s(@(t,y) cf_eqs(t,y), tspan, y0);

% time series
c = y(:,1);
f = y(:,2);

%%% plots =================================================================

%%% plot null clines
interval = [0,1];
hold on; box on;
fimplicit(Cp, interval,'b','Linewidth',2)
fimplicit(Fp, interval,'r','Linewidth',2)
xlabel('C')
ylabel('F')

%%% plot path of eq point
figure()
hold on; box on;
plot(Epts(:,1),Epts(:,2),'Linewidth',2)
xlabel('C')
ylabel('F')
xlim([0 1])
ylim([0 1])
% plot(c,f,'-.','Linewidth',2)
% scatter(c(end),f(end),'kx','Linewidth',2)
title('Path of equilibrium point for 0 < beta < 30')

%%% equilib population values as function of b
figure()
hold on; box on;
plot(Beta,Epts(:,1),'Linewidth',2)
plot(Beta,Epts(:,2),'Linewidth',2)
xlabel('b')
ylabel('Scaled populations at equilibrium')
legend('C','F')
title('Equilbrium values as function of maximum growth rate (beta))')



%%% test line changing color
xx = Epts(:,1);
yy = Epts(:,2);

figure()
clf
% Color changes along x axis
% sp(1) = subplot(3,1,1);
% col = [xx(1:end-1)', NaN]; % Must end in NaN to filling a solid
col = [Beta(1:end-1), NaN]; % Must end in NaN to filling a solid
patch(xx,yy,col, 'EdgeColor', 'interp', 'LineWidth', 3)
axis tight
box on
grid on
title('Path of equilibrium point for 0 < b < 30')
xlim([0 1])
ylim([0 1])
xlabel('C')
ylabel('F')
h = colorbar;
ylabel(h, 'beta')

%%% functions =============================================================

%%% cf ode function
function yp = cf_eqs(t,y)
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

%%% function that returns eigenvalues of a point (c,f)
function v = cfeigs(x)
global k beta r d b mu eta lambda q n

c = x(1);
f = x(2);

%%% big terms
A = (c*k*eta + mu);
B = (lambda^n)/(lambda^n + b^n*A^n);
D = (b^n)*(1-c-f)*k*n*beta*eta*(lambda^n)*(A^(n-1))/((lambda^n + (b^n)*A^n)^2);


J = zeros(2,2);

J(1,1) = -d - c*D - c*beta*B + (1-c-f)*beta*B;
J(1,2) = -c*beta*B;
J(2,1) = f*k*q*eta*lambda/(A^2) + f*D - f*(r + beta*(1-B));
J(2,2) = -d - q*lambda/A + (1-c-f)*(r + beta*(1-B)) - f*(r + beta*(1-B));

v = eigs(J);

end