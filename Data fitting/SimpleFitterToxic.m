%%% fitting simple CF odes with toxic oxygen
%%% 2/19/2020

close all;

data = xlsread('C:\Users\peter\OneDrive\Desktop\cyst fib\julia stuff\ODEs\Data fitting\cf data','Rescaled');
tdata = data(:,1);
% cdata = data(:,2)/100;
% fdata = data(:,3)/100;
cdata = data(:,2);
fdata = data(:,3);

%%% =======================================================================

% fixed parameters
global lambda rcmin rfmin t_treat alpha beta

lambda = 0.22;
rcmin = 0;
rfmin = 0;
t_treat = 28.;
alpha = 1.0;
beta = 1.0;

%%% =======================================================================
% do optimization here
% options = optimset('MaxFunEvals',5000,'Display','iter','UseParallel',true);
% options = optimoptions(options,'UseParallel',true);
options = optimset('MaxFunEvals',5000,'Display','iter');
% options = optimset('MaxFunEvals',5000);

% parameters to fit
Ec = 9; % try < 16
Ac = 0.0128;
nc = 1.0;

rf = 11;

d = 0.1000;
dc = d;
df = 0.15;

ep = 0.5743;
mu = 1.2079;
keta = 0.6281;
q = 0.0;

c0 = 0.9313;
f0 = 0.0561;
w0 = 0.0812;

p = [Ec,rf,dc,df,ep,q,c0,f0,w0];

A = []; b = []; Aeq = []; Beq = []; 
lb = [0 0 0 0 0 0 0.95 0 0.8*cdata(1) 0.08*fdata(1) 0.1]; 
ub = [Inf Inf Inf Inf Inf Inf Inf Inf 1.3*cdata(1) 1.3*fdata(1) 0.22];


tic
% [p,fval,flag,output] = fmincon(@cf_err,p,A,b,Aeq,Beq,lb,ub,[],options,tdata,cdata,fdata);
[p,fval,flag,output] = fminsearch(@cf_err,p,options,tdata,cdata,fdata);
toc


% solve ode's
y0 = [p(6), p(7), p(8)];
tspan = [0 40];
[t, y] = ode15s(@(t,y) cf_eqs(t,y,p), tspan, y0);
J = cf_err(p,tdata,cdata,fdata)
p

% figure(); hold on; box on;
% plot(t,y(:,1),'b','Linewidth',2)
% plot(t,y(:,2),'r','Linewidth',2)
% plot(tdata,cdata,'bx', 'LineWidth',2)
% plot(tdata,fdata,'rx', 'LineWidth',2)
% xlabel('Time (days)')
% ylabel('Population density')

hold on; box on;
plot(t,y(:,1),'Linewidth',2)
plot(t,y(:,2),'Linewidth',2)
plot(tdata,cdata,'bx', 'LineWidth',2)
plot(tdata,fdata,'rx', 'LineWidth',2)
xlabel('Time (days)')
ylabel('Population density')
title('Climax and Attack Populations')
legend('C model','F model','C data','F data','Location','w')

% figure()
% plot(t,y(:,3),'Linewidth',2)
% xlabel('Time (days)')
% ylabel('Oxygen')
% title('Oxygen')

%%% =======================================================================
%%% plot rc(w)
w = y(:,3);
rc = (p(1)*w.^p(3))./(p(2)^p(3) + w.^p(3));

% figure()
% plot(t,rc)
% xlabel('Time (days)')
% ylabel('Climax growth rate')
% title('Climax growth rate')

%%% Functions =============================================================
[Ec,rf,dc,df,ep,q,c0,f0,w0];
%%% cf ode function
function yp = cf_eqs(t,y,p)
global lambda rcmin t_treat alpha beta

Ec = p(1);
rf = p(2);
ep = 0; mu = 1; keta = 0.6281;
q = p(6);

dc = (3); df = (4);

if t >= t_treat
    ep = p(5);
end

c = y(1);
f = y(2);
w = y(3);

yp = zeros(3,1);

yp(1) = Ec*c*(1 - c - alpha*f) - dc*c;
yp(2) = rf*f*(1 - f - beta*c) - df*f - ep*f - q*f*w;
yp(3) = lambda - mu*w - keta*c*w;

% hold on
% scatter(t,(Ec*w^nc/(Ac^nc + w^nc)+ rcmin),'bx')
% scatter(t,(Ef*Af^nf/(Af^nf + w^nf) + rfmin),'ro')

end

%%% objective function for cf_fitter
function J = cf_err(p,tdata,cdata,fdata)
y0 = [0.8138, 0.1234, 0.15];
y0(1) = p(7); y0(2) = p(8); y0(3) = p(9);
[t,y] = ode15s(@cf_eqs,tdata,y0,[],p);

% odata = 0.22*ones(length(tdata),1);

errx = y(:,1) - cdata;
erry = y(:,2) - fdata;
% erro = y(:,3) - odata;

% J = errx'*errx + erry'*erry + erro'*erro;
J = errx'*errx + erry'*erry;
end