%%% numerically find steady states of the big model using fsolve
%%% model is 

%%% c' = rc*c*(1-c-f) - d*c
%%% f' = rf*f*(1-c-f) - d*f - q*lambda*f/(mu + eta*k*c)

%%% rc = (beta*lambda^n)/(lambda^n + b^n*(mu + eta*k*c)^n)
%%% rf = (r + beta*(1 - (lambda^n)/(lambda^n + b^n*(mu + eta*k*c)^n)))

%%% this one is multiple trajectories
%%% 5/24/21

close all;

%%% Parameters ============================================================
global k beta r d b mu eta lambda q n
k = 10^10;

beta = 16.6;
r = 0.004;
d = 0.6;
b = 13.4;
n = 2.6;

mu = 200*23*60*24;
eta = 3.1e-4;

lambda = 9.2e7;
q = 3e-5;



%%% find steady states numerically
fun = @SStates;
options = optimset('Display','off');

y0 = [0.8,0.4];
pts = 20;

C = linspace(0,1,pts);
F = linspace(0,1,pts);

hold on
xlim([0 1])
ylim([0 1])
xlabel('C')
ylabel('F')

% for i = 1:length(C)
%     for j = 1:length(F)
%         [i j]
%         
%         x0 = [C(i), F(i)];
%         xc = fsolve(fun,x0,options);
%         
%         scatter(xc(1),xc(2))
% 
%     end
% end



xc = fsolve(fun,[1,0]);

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