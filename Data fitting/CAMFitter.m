%%% fitter with good parameter values
%%% 3/9/2020

close all;

data = xlsread('C:\Users\peter\OneDrive\Desktop\cyst fib\julia stuff\ODEs\Data fitting\cf data','Sheet1');
tdata = data(:,1);
cdata = data(:,2)/100;
fdata = data(:,3)/100;
% cdata = data(:,2);
% fdata = data(:,3);

%%% =======================================================================

% fixed parameters
global k lambda t_treat N0

N0 = 6.7e8;
lambda = 1.32;
t_treat = 28.;

%%% =======================================================================
% do optimization here
options = optimset('MaxFunEvals',5000,'Display','iter');
% options = optimset('MaxFunEvals',5000);

% parameters to fit
r = 24.1532;

beta = 23.1363; % try < 16
b = 3.81e-01;
n = 1.3858;

d = 6.0;

ep = 0.6357;
mu = 1.787e7; % 1/5 min

k = 10^10;
eta = 5.7659e-5;
q = 0.000246;

frac = 0.8901;
c0 = frac*N0;
f0 = (1 - frac)*N0;
x0 = 200;

lambda = mu*x0;

p = [b,n,q,ep,frac,beta,r,mu,eta,d];


A = []; b = []; Aeq = []; Beq = [];
lb = zeros(5,1);
ub = [1 5 2000 200 1];


tic
% [p,fval,flag,output] = fmincon(@cf_err,p,A,b,Aeq,Beq,lb,ub,[],options,tdata,cdata,fdata);
[p,fval,flag,output] = fminsearch(@cf_err,p,options,tdata,cdata,fdata);
toc


% solve ode's
y0 = [c0; f0; x0];
tspan = [0 40];
[t, y] = ode15s(@(t,y) cf_eqs(t,y,p), tspan, y0);
J = cf_err(p,tdata,cdata,fdata)
p

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

figure()
plot(t,y(:,3),'Linewidth',2)
xlabel('Time (days)')
ylabel('Oxygen')
title('Oxygen')

%%% =======================================================================
% %%% plot rc(w)
% w = y(:,3);
% rc = (p(1)*w.^p(3))./(p(2)^p(3) + w.^p(3));
% 
% figure()
% plot(t,rc)
% xlabel('Time (days)')
% ylabel('Climax growth rate')
% title('Climax growth rate')

%%% Functions =============================================================

%%% cf ode function
function yp = cf_eqs(t,y,p)
global k lambda t_treat

b = p(1); 
n = p(2);
q = p(3);
ep = 0; 
beta = p(6);
r = p(7);
mu = p(8);
eta = p(9);
d = p(10);

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
% erry = Ft - fdata;

% J = errx'*errx + erry'*erry;
J = errx'*errx;
end