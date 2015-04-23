function [Tout,Xout,DXout, info] = vdpsolver(T0,Tfinal,X0,DX0,tol,A,Mu,omega)
% VDPSOLVER gives an approximation to the Van Der Pol equation.
%
% This function uses the vectorized Runge-Kutta-Fehlberg method to provide
% a numerical solution to the Van Der Pol equation which is of the form:
% y'' - mu*(1-y^2)*y' + y - A*sin(omega*t) = 0
% where t0 <= t <= tf and y(0) = a0, y'(0) = a1.
%
% INPUT(S):
%
% T0: The initial starting t value.
%
% TFINAL: The final t value.
%
% X0: A number such that y(t0) = X0.
%
% DX0: A number such that y'(t0) = DX0
%
% TOL: The amount of tolerance we want to approximate our differential 
% equation to (we use this for our Runge-Kutta-Fehlberg method).
%
% A: Amplitude of our 'forcing' function.
%
% MU: A scaling parameter.
%
% OMEGA: The angular velocity of our forcing function.
%
% OUTPUT(S):
%
% TOUT: A list of mesh-points that were used to get within the right 
% tolerance level in our Runge-Kutta-Fehlberg method. 
% 
% XOUT: An array that lists all the approximations to y(t) at the mesh 
% points, i.e. y(ti) = Xout(i).
%
% DXOUT: An array that lists all the approximations to y'(t) at the mesh 
% points, i.e. y'(ti) = DXout(i).
%
% INFO: A string of characters that notifies the user how long the whole
% algorithm took to run.
%
%
% Written by: Mario Perales
% April 21, 2015

%--------------------%
% Initial Parameters %
%--------------------%
tic;
alpha = [X0;DX0];
tspan = [T0 Tfinal];
hmin = eps;
hmax = 3/(Mu+.1);
stepsize = [hmin hmax];

% We need to create a vector that just takes the parameters and outputs
% a vector with the parameters @(t,x). We do this by calling a function
% which does just that.

FunFunc = @(t,x) vdp(t,x,Mu,A,omega);

[w,Tout] = RKFv(FunFunc, tspan, alpha, tol, stepsize);

Xout = w(1,:); % Seperating our double column vector. %
DXout = w(2,:); % Seperating our double column vector. %

Xout = Xout'; % Just to make it look nicer...%
DXout = DXout';

figure
subplot(211)
plot(Tout,w(1,:)) % plot t vs x
subplot(212)
plot(Tout(1:end-1),diff(Tout)) % plot t vs h (step size)

time = toc;
str = num2str(time);
info = ['The time it took for this method to run was ' str ' seconds!'];

end

function xp = vdp(t,x,Mu,A,omega)
% on return: a vector of differential equations.
%
% Using the u substitution where x(1) = y(t) and x(2) = y'(t) we get the
% following differential equation vector that will be returned and is a
% valid vector that works with our differential equation methods.
 
    xp(1,1) = x(2);
    xp(2,1) = Mu*(1-x(1)^2)*x(2) - x(1) + A*sin(omega*t);
    
end

function [w,t,flg] = RKFv(FunFcnIn, Intv, alpha, tol, stepsize)
% On input: 
%   FunFcnIn is the name of function to be integrated
%   interv is the interval to be integrated over
% 
%   The problem: y' = f(t,y), y(a) = alpha, a<= t <= b
%   where Intv = [a b], and a call to FunFcnIn with 
%   argument (t, y) returns f(t,y).
%
%   RKF uses the Runge-Kutta-Fehlberg method to solve
%   the above problem, to a given tolerance. 
%   hmin and hmax are the minimum and maximum step sizes.
% 
% this is the vector version of RKF.
%
% On output
%   t contains the (unequi-spaced) mesh points, w the
%   function values at these points. 
%   flg is success flag. If flg == 0 method succeeded,
%   and if flg ~=0 method failed.
%
% Written by Ming Gu for Math 128A, Fall 2008
% 
[FunFcn,msg] = fcnchk(FunFcnIn,0);
if ~isempty(msg)
    error('InvalidFUN',msg);
end
flg  = 0;
a    = Intv(1);
b    = Intv(2);
hmin = stepsize(1);
hmax = stepsize(2);
h    = min(hmax,b-a);
if (h<=0)
    msg = ['Illegal hmax or interval'];
    error('InvalidFUN',msg);
end
w    = alpha;
t    = a;
c    = [0;      1/4; 3/8;       12/13;      1;  1/2];
d    = [25/216  0    1408/2565  2197/4104 -1/5  0   ];
r    = [1/360   0   -128/4275  -2197/75240 1/50 2/55];
KC   = [zeros(1,5); 
        1/4,       zeros(1,4); 
        3/32       9/32       zeros(1,3);
        1932/2197 -7200/2197, 7296/2197, 0,         0;
        439/216   -8          3680/513  -845/4104,  0;
          -8/27    2         -3544/2565  1859/4104 -11/40];
K    = zeros(length(alpha),6);
%
% main loop
%
while (flg == 0)
    ti   = t(end);
    wi   = w(:,end);
    for j = 1:6
        K(:,j) = h*FunFcn(ti+c(j)*h,wi+K(:,1:5)*KC(j,:)');
    end
%
% accept approximation
%
    R = max(eps,norm(K*r')/h);
    if (R <= tol)
        t = [t; ti+h];
        w = [w, wi+K*d'];
    end
%
% reset stepsize.
%
    delta = 0.84*(tol/R)^(1/4);
    delta = min(4,max(0.1,delta));
    h     = min(hmax, delta * h);
    h     = min(h, b-t(end));
    if (abs(h)<eps)
        return;
    end
    if (h < hmin)
        flg = 1;
        return;
    end
end
end
