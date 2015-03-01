function [rts, info] = cubic23604571(C)
info = [];
tolerance = .000001; % We want six digits of accuracy. %
rts = [];

NMaxBis1 = 60;
NMaxNewton = 25;
numOfRoots = 1;

i = length(C) - 1;

% How many roots do we have? %
while i > 0;
    if(C(1,i) == 0);
        i = i - 1;
    else
        rts = zeros(i, 1);
        degree = i;
        break;
    end;
end;

% Find max range% 


% approxRoot = -1;
% N = 10

% Laguerre's Method %

% while(i <= N);
%     if(abs((C(1,1) * approxRoot^3 + C(1,2) * approxRoot^2 + C(1,3) * approxRoot + C(1,4))) > tolerance)
%         fApproxRoot = C(1,1) * approxRoot^3 + C(1,2) * approxRoot^2 + C(1,3) * approxRoot + C(1,4);
%         dfApproxRoot = 3 * C(1, 1) * approxRoot^2 + 2 * C(1,2) * approxRoot + C(1,3); 
%         ddApproxRoot = 6 * C(1,1) * approxRoot + 2 * C(1,2);
%     
%         G = dfApproxRoot/fApproxRoot;
%         H = G^2 - ddApproxRoot/fApproxRoot;
%     
%         positivea = degree/(G + sqrt((degree - 1)*(degree * H - G^2)));
%         negativea = degree/(G - sqrt((degree - 1)*(degree * H - G^2)));
% 
%         a = max(positivea, negativea);
%     
%         approxRoot = approxRoot - a;
%         i = i + 1;
%         
%     else
%         rts(numOfRoots, 1) = approxRoot;
%         numOfRoots = numOfRoots + 1;
%         break;
%     end;
% end;
% end
        
    

% % Newton's Method (taken from page 66 of Numerical Analysis Textbook [9th Ed.] by Burden et Faires) %

% i = 1;
% 
% approxRoot = 1; % Just a place holder %
% 
% while i <= NMaxNewton;
%     fApproxRoot = C(1,1) * approxRoot^3 + C(1,2) * approxRoot^2 + C(1,3) * approxRoot + C(1,4);
%     dfApproxRoot = 3 * C(1, 1) * approxRoot^2 + 2 * C(1,2) * approxRoot + C(1,3);
%     
%     p = approxRoot - fApproxRoot/dfApproxRoot;
%     
%     if abs(p - approxRoot) < tolerance;
%         rts(numOfRoots, 1) = p;
%         numOfRoots = numOfRoots + 1;
%         % Might be useful later :) %
%         dividedPoly = deconv(C,[1, -p]);
%         
%         % We must now check the multiplicity of the roots (see p. 82) %
%         dFP = 3 * C(1,1) * p^2 + 2 * C(1,2) * p + C(1,3);
%         if(abs(dFP) < tolerance*10); %*Make note of the error here*%
%             rts(numOfRoots, 1) = p;
%             numOfRoots = numOfRoots + 1;
%             
%             ddFP = 6 * C(1,1) * p + 2 * C(1,2);
%             if(ddFP < tolerance*100 && length(rts) == 3);
%                 rts(numOfRoots, 1) = p;
%                 numOfRoots = numOfRoots + 1;
%             break;
%             end;
%         break;    
%         end;
%     break;
%     end;
%     i = i + 1;
%     approxRoot = p;  
% end;

% Bisection Method %
% a = something
% b = something
% i = 1;
%     FA = C(1,1) * a^3 + C(1,2) * a^2 + C(1,3) * a + C(1,4);
% 
% while i <= NMaxBis1
%    
%    p = a + (b - a)/2;
%    FP = C(1,1) * p^3 + C(1,2) * p^2 + C(1,3) * p + C(1,4);
%     
%    if FP == 0 || abs((b-a)/2) < tolerance;
%        rts(numOfRoots, 1) = p;
%        numOfRoots = 1 + numOfRoots;
%        break;
%    end
%    
%    i = i + 1;
%    
%    if (FA * FP) > 0;
%        a = p;
%        FA = FP;       
%    else
%        b = p;
%    end
% 
% end

% %Muller's Method%

% fa = C(1,1) * a^3 + C(2,1) * a^2 + C(3,1) * a + C(4,1);
% fb = C(1,1) * b^3 + C(2,1) * b^2 + C(3,1) * b + C(4,1);
% fc = C(1,1) * c^3 + C(2,1) * c^2 + C(3,1) * c + C(4,1);
% 
% h1 = b - a;
% h2 = c - b;
% delta1 = (fb - fa)/h1;
% delta2 = (fc - fb)/h2;
% d = (delta2-delta1)/(h2 + h1);
% i = 3;
% 
% while i <= NMaxMul
%         b = delta2 + h2*d;
%         D = (b^2 - 4 * fc * d)^(1/2);
%         
%         if abs(b-D) < abs(b+D);
%             E = b + D;
%         else
%             E = b - D;
%         end
%         
%         h = -2 * fc/E;
%         p = c + h;
%         
%         if abs(h) < tolerance;
%             numOfRoots;
%             rts(numOfRoots, 1) = p;
%             break;
%         end;
%         
%         a = b;
%         b = c;
%         c = p;
%         h1 = b - a;
%         h2 = c - b;
%         delta1 = (fb-fa)/h1;
%         delta2 = (fc - fa)/h2;
%         d = (delta2 - delta1)/(h2+h1);
%         i = i + 1;
% end;




