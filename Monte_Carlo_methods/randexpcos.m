function X = randexpcos(m,n)

% generates X ~ K*exp(cos(x))


R = m*n;
X = zeros(1,R);
reject = (1:R);
while R > 0
   V = rand(1,R); 
   % >> Now generate Xg:
   Xg = quantileF(V);
   X(reject) = Xg;
   U = rand(size(Xg));
   % >> Now define rejection criterium:
   fXoverMgX = f(Xg) / g(Xg);
   reject = reject(U > fXoverMgX);
   R = length(reject);
end
X = reshape(X,m,n);
end

function y = f(x)
    y = NaN(size(x)); 
    mask = x >= 0 & x <= 1; 
    y(mask) = exp(cos(2 * pi * x(mask))); 
end

function y = g(x)
    y = NaN(size(x)); 
    mask = x >= 0 & x <= 1; 
    y(mask) = (1 + 10 * abs(x(mask) - 0.5)) / 2;
end

function Q = quantileF(p)
    Q = arrayfun(@(x) calcQuantile(x), p);
end

function q = calcQuantile(x)
    if x >= 0 && x <= 1/2
        q = (6 - sqrt(36 - 70 * x)) / 10;
    elseif x > 1/2 && x <= 1
        q = (4 + sqrt(16 - 10 * (5 - 7 * x))) / 10;
    else
        q = NaN;
    end
end