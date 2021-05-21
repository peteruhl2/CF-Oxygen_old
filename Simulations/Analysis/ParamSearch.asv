%%% parameter seach for the MessyModel
%%% 5/18/21

%%% C' = (beta*x/(b + x))*c*(1 - (c + f)/k) - dc*c;
%%% F' = (r + beta*(1 - x/(b + x)))*f*(1 - (c + f)/k) - df*f  - q*x*f;
%%% X' = lambda - mu*x - eta*x*c;

%%% =======================================================================

% fixed parameters
global k N0
 
N0 = 6.7e8;
k = 10^10;

%%% Ode parameters ========================================================

%%% parameters to serach
Q = 10.^linspace(-10,-1);
B = linspace(1,20);
Beta = linspace(1,20);

results = zeros(length(Q),length(B));

for i = 1:length(Q)
    for j = 1:length(Beta)
        
    beta = Beta(j);
    r = 18;
    d = 0.6;
    q = Q(i);

    lambda = 9.7e7;
    mu = 200*23*60*24;
    eta = 2.2e-4;
    b = 12.4;

    p = [beta,r,d,q,...
         lambda,mu,eta,b];

    % square root terms
    rad1 = sqrt(d + 4*b*q)/sqrt(d);
    rad2 = sqrt(d)*sqrt(d+4*b*q);

    %%% exclusion values
    cx = (-d*k*lambda +k*beta*lambda - b*d*k*mu)/(b*d*k*eta + beta*lambda);
    fx = k*(-q*lambda^2 - d*lambda*mu - b*q*lambda*mu - b*d*mu^2 + b*beta*mu^2)/(b*beta*mu^2);
    xx = (-b*d*k*eta - beta*lambda)/(d*k*eta - k*beta*eta - beta*mu);

    % coexistence 1
    cx1 = (lambda + lambda*rad1 - 2*b*mu)/(2*b*eta);
    fx1 = -((1 + rad1)*beta*lambda + b*(3*d*k*eta + k*eta*rad2 - 2*beta*(k*eta + mu)))/(2*b*beta*eta);
    xx1 = (-d + rad2)/(2*q);

    % coexistence 2
    cx2 = (lambda - lambda*rad1 - 2*b*mu)/(2*b*eta);
    fx2 = -((1 - rad1)*beta*lambda - b*(-3*d*k*eta + k*eta*rad2 + 2*beta*(k*eta + mu)))/(2*b*beta*eta);
    xx2 = (-d + rad2)/(2*q);


    %%% stability =============================================================
    vc = eig(jac([cx,0,xx],p));
    vf = eig(jac([0,fx,lambda/mu],p));
    vx = eig(jac([0,0,lambda/mu],p));
    vx1 = eig(jac([cx1,fx1,xx1],p));
    vx2 = eig(jac([cx2,fx2,xx2],p));

        [vc<0 vf<0 vx<0 vx1<0 vx2<0]
    S = sum([vc<0 vf<0 vx<0 vx1<0 vx2<0]);

    if sum(S>2) > 1
        results(i,j) = 1;
    end
          
       
    end
end


find(results>0)



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
J(1,3) = -(beta*x/((b + x)^2))*c*(1 - (c + f)/k) + beta*c*(1 - (c + f)/k)/(b + x);

J(2,1) = -(beta*f/k)*(1 - (c+f)/k);
J(2,2) = -d - q*x + (1 - (c+f)/k)*(1 - x/(b+x))*beta - (beta*f/k)*(1 - x/(b+x));
J(2,3) = -q*f + beta*f*(1 - (c+f)/k)*(x/((b+x)^2) - 1/(b+x));

J(3,1) = -eta*x;
J(3,2) = 0;
J(3,3) = -eta*c - mu;

end