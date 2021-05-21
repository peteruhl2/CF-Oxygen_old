%%% messy model oxygen model
%%% 5/13/21

%%% C' = (beta*x/(b + x))*c*(1 - (c + f)/k) - dc*c;
%%% F' = (r + beta*(1 - x/(b + x)))*f*(1 - (c + f)/k) - df*f  - q*x*f;
%%% X' = lambda - mu*x - eta*x*c;

close all;

%%% =======================================================================

% fixed parameters
global k N0

N0 = 6.7e8;
k = 10^10;

%%% Ode parameters ========================================================

beta = 16.3838;
r = 16;
d = 0.6;
q = 1.0e-5;

lambda = 9.7e7;
mu = 200*23*60*24;
eta = 4.2e-3;
b = 12.4;

n = 1;

p = [beta,r,d,q,...
     lambda,mu,eta,b,n];
 
F0 = 10.^linspace(1,9,5);
C0 = 10.^linspace(1,9,5);

for i = 1:length(C0)
    for j = 1:length(F0)
        frac = 0.85;
        c0 = C0(i);
        f0 = F0(j);
        x0 = 14;
        y0 = [c0; f0; x0];

        tspan = [0 600];
        [t, y] = ode15s(@(t,y) cf_eqs(t,y,p), tspan, y0);

        hold on; box on;
        plot(log10(y(:,1)),log10(y(:,2)),'Linewidth',2)
        scatter(log10(y(end,1)),log10(y(end,2)),'Linewidth',2)
%         plot((y(:,1)),(y(:,2)),'Linewidth',2)
%         scatter((y(end,1)),(y(end,2)),'Linewidth',2)
        xlabel('C')
        ylabel('F')

        if 1 == 1
            scatter(0,0,'Linewidth',2)
        end
    end
end

figure()
hold on; box on;
plot(t,log10(y(:,1)),'b','Linewidth',2)
plot(t,log10(y(:,2)),'r','Linewidth',2)
xlabel('Time (days)')
ylabel('Absolute Abundance')
title('Climax and Attack Populations')
legend('C model','F model')

% figure()
% plot(t,y(:,3),'Linewidth',2)
% xlabel('Time (days)')
% ylabel('Oxygen (\muM)')
% title('Oxygen')

%%% equilibria ============================================================

% time series
c = y(:,1);
f = y(:,2);
x = y(:,3);

v = eig(jac(y(end,:),p))
v2 = eig(jac(extinct,p))

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
n = p(9);

dc = d;
df = d;

c = y(1);
f = y(2);
x = y(3);

yp = zeros(3,1);

yp(1) = (beta*x^n/(b^n + x^n))*c*(1 - (c + f)/k) - dc*c;
% yp(2) = (beta*(1 - x/(b + x)))*f*(1 - (c + f)/k) - df*f  - q*x*f;
yp(2) = (r*(1 - x^n/(b^n + x^n)))*f*(1 - (c + f)/k) - df*f  - q*x*f;
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
n = p(9);

% coefficient terms
A = (1 - (c+f)/k);
B = x^n/(b^n + x^n);
E = (n*x^(2*n-1))/((b^n+x^n)^2);
D = (n*x^(n-1))/((b^n+x^n));

J = zeros(3,3);

J(1,1) = -d + beta*A*B + c*beta*B/k;
J(1,2) = -c*beta*B/k;
J(1,3) = -c*A*E + c*A*beta;

J(2,1) = -f*r*(1-B)/k;
J(2,2) = -d - q*x + A*r*(1-B) - f*r*(1-B)/k;
J(2,3) = -f*q + f*A*r*(E - D);

J(3,1) = -eta*x;
J(3,2) = 0;
J(3,3) = -mu - eta*c;

end
