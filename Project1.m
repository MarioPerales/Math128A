function [rts, info] = cubic(C)
% CUBIC finds the roots of a polynomial of degree three or less.
%
% This function uses Muller's Method for finding roots of a generalized
% function. In addition, it also uses Horner's method to check for
% multiplicity of the roots. 
%
% This version of Muller's Method is taken from the book:
% Title: Numerical Analysis (9th Edition)
% Author(s): Richard L. Burden & J. Douglas Faires
% Page(s): 97 - 98
%
% INPUT(S):
%
% C: An array of integers that satisfies the property that length(C) = 4. 
% Any higher would actually give you the roots of the polynomial, but it 
% would at most return three roots. This function will find roots of
% quadratics by simply making the first entry of the array 0.
% 
% *NOTE* We do not need the user to input three initial guesses. This
% function correctly checks the bounds for the roots to be in.
%
% OUTPUT(S):
%
% RTS: An array that contains all the roots of the polynomial given, in no
% particular order. Repeated roots are printed as many times as they
% repeat.
%
% INFO: A string of characters that notify the user what went wrong with
% the code or any other relevant info. Currently holds how long it took to
% find the roots of a function, as well as outputting a message if the 
% method is unsuccessful.
%
% *IMPORTANT NOTE* - Currently has issues when the coeffecients are large.
%
% Written by: Mario Perales
% March 3, 2015

%--------------------%
% Initial Parameters %
%--------------------%
tic;
info = [];
rts=[];
tolerance = .00000001; % We want six digits of accuracy (but let's shoot for more!) %
degree = length(C) - 1;
NMaxMul = 10000; % Max Iterations for Muller's Method %
doubleFlag = 0; % Multiplicity flags %
tripleFlag = 0;
numOfRoots = 1;

% What degree is our polynomial? %

i = 1; 
while i < 4;
    if(C(1,i) == 0);
        degree = degree - 1;
        i = i + 1;
    else
        break;
    end;
end;

% The bounds where the roots are guranteed to be in. For more information %
% on obtaining these bounds, see:
% http://www.cs.iastate.edu/~cs577/handouts/polyroots.pdf %

rho1 = min(degree * abs(C(1,length(C))/C(1,length(C) - 1)), (abs(C(1,length(C))/(C(1,length(C) - degree))))^(1/degree));
rho2 = 1 + max(max(abs(C(1,length(C))/C(1, length(C) - degree)), abs(C(1,length(C) - 1)/C(1, length(C) - degree))), max(abs(C(1,length(C)-2)/C(1, length(C) - degree)), abs(C(1,length(C) - 3)/C(1, length(C) - degree))));

% Setting up the first three guesses near the bounds so that we are 
% guruanteed to converge. %

rho = min(rho1, rho2);
oom = log(abs(rho))./log(10); % oom = Order of Magnitude, helps to correctly place the three guesses %
a = 10^(oom - oom*(1/4));
b = 10^(oom - oom*(1/2));
c = 10^(oom);

if(oom == 0);
    a = rho - rho*(1/4);
    b = rho - rho*(1/2);
    c = rho;
end;  

%-------------------%
% Muller's Method   %
%-------------------%

h1 = b - a;
h2 = c - b;
delta1 = (polyval(C,b) - polyval(C,a))/(h1);
delta2 = (polyval(C,c) - polyval(C,b))/(h2);
d = ((delta2)-(delta1))/((h2) + (h1));
i = 3;

while i <= NMaxMul
        b = delta2 + h2*d;
        D = sqrt(b^2 - 4 * polyval(C,c) * d);
        if abs(b-D) < abs(b+D);
            E = b + D;
        else
            E = b - D;
        end
        
        h = -2 * polyval(C,c)/E;
        p = c + h;
        if abs(h) < tolerance;
            rts(numOfRoots, 1) = p;
            numOfRoots = numOfRoots + 1;
            break;
        end;
        
        a = b;
        b = c;
        c = p;
        h1 = b - a;
        h2 = c - b;
        delta1 = (polyval(C,b) - polyval(C,a))/(h1);
        delta2 = (polyval(C,c) - polyval(C,b))/(h2);
        d = (delta2 - delta1)/(h2+h1);
        i = i + 1;
end;

% We found one root, now we must check its multiplicity. %

dividedPoly = deconv(C,[1, -p]); % Synthetic Division %

%------------------%
% Horner's Method  %
%------------------%

y = dividedPoly(1);
z = dividedPoly(1);

j = 2;

while j < length(dividedPoly);
    y = y * p + C(j);
    z = p * z + y;
    j = j + 1;
end;

y = p * y + C(length(dividedPoly));

if(y == 0); % We check if dividing the polynomial by the same root gives us a remainder, thus implying no multiplicity. %
    doubleFlag = 1; % We hit a double root. %
    rts(numOfRoots, 1) = p;
    numOfRoots = numOfRoots + 1;
    
    % Repeat process as before. %
    dividedPoly2 = deconv(dividedPoly,[1, -p]);
    
    %------------------%
    % Horner's Method  %
    %------------------%
    
    y = dividedPoly2(1);
    z = dividedPoly2(1);

    j = 2;
    while j < length(dividedPoly2);
        y = y * p + dividedPoly1(j);
        z = p * z + y;
        j = j + 1;
    end;

    y = p * y + dividedPoly(length(dividedPoly));
    
    if(y == 0);
        rts(numOfRoots, 1) = p;
        numOfRoots = numOfRoots + 1;
        tripleFlag = 1; % We hit three roots. %
    end;
end;

if(i > NMaxMul);
    info = 'Well, this is awkward...We did not converge :(';
end;    
    

% If we didn't find any multiple roots, we now find roots of the factored
% out polynomial the same way we did before.                               

if(tripleFlag ~= 1 && doubleFlag ~= 1);
    % We find the bounds of this new polynomial so we don't converge to our
    % previous root
    
    rho1 = min(degree * abs(dividedPoly(1,length(dividedPoly))/dividedPoly(1,length(dividedPoly) - 1)), (abs(dividedPoly(1,length(dividedPoly))/(dividedPoly(1,length(dividedPoly) - degree + 1))))^(1/degree));
    rho2 = 1 + max(max(abs(dividedPoly(1,length(dividedPoly))/dividedPoly(1, length(dividedPoly) - degree + 1)), abs(dividedPoly(1,length(dividedPoly) - 1)/dividedPoly(1, length(dividedPoly) - degree + 1))), abs(dividedPoly(1,length(dividedPoly) - 2)/dividedPoly(1, length(dividedPoly) - degree + 1)));

    rho = min(rho1, rho2);
    oom = log(abs(rho))./log(10);
    a = 10^(oom - oom*(1/4));
    b = 10^(oom - oom*(1/2));
    c = 10^(oom);
    
    if(oom == 0);
        a = rho - rho*(1/2);
        b = rho - rho*(1/3);
        c = rho;
    end;  
    

    %-------------------%
    % Muller's Method   %
    %-------------------%
    
    h1 = b - a;
    h2 = c - b;
    delta1 = (polyval(dividedPoly,b) - polyval(dividedPoly,c))/h1;
    delta2 = (polyval(dividedPoly,c) - polyval(dividedPoly,b))/h2;
    d = (delta2-delta1)/(h2 + h1);
    i = 3;

    while i <= NMaxMul
        b = delta2 + h2*d;
        D = (b^2 - 4 * polyval(dividedPoly,c) * d)^(1/2);
        
        if abs(b-D) < abs(b+D);
            E = b + D;
        else
            E = b - D;
        end
        
        h = -2 * polyval(dividedPoly,c)/E;
        p1 = c + h;

        if abs(h) < tolerance;
            rts(numOfRoots, 1) = p1;
            numOfRoots = numOfRoots + 1;
            break;
        end;
        
        a = b;
        b = c;
        c = p1;
        h1 = b - a;
        h2 = c - b;

        delta1 = (polyval(dividedPoly,b) - polyval(dividedPoly,c))/h1;
        delta2 = (polyval(dividedPoly,c) - polyval(dividedPoly,b))/h2;
        d = (delta2 - delta1)/(h2+h1);
        i = i + 1;
    end;
    
    dividedPoly2 = deconv(dividedPoly,[1, -p1]);
    
    %------------------%
    % Horner's Method  %
    %------------------%

    y = dividedPoly2(1);
    z = dividedPoly2(1);

    j = 2;
    while j < length(dividedPoly2);
        y = y * p1 + dividedPoly2(j);
        z = p1 * z + y;
        j = j + 1;
    end;

    y = p1 * y + dividedPoly2(length(dividedPoly2));

    if(y == 0);
        rts(numOfRoots, 1) = p1;
        numOfRoots = numOfRoots + 1;
        tripleFlag = 1; % We hit three roots. %
    end;
    
end;

% This is the case where we have factored out our polynomial as much as we
% can, thus making the final root trivial to find. %

if(doubleFlag ~= 1 && tripleFlag ~= 1 && C(1,1) ~= 0);
    rts(numOfRoots,1) = -(dividedPoly2(1,2));
end;

time = toc;
str = num2str(time);

if(isempty(info) == 1);
    info = ['The time it took to find these roots was ' str ' seconds'];
end;

end
