studnumber = 601156;
if exist('realrandom') ~= 1, realrandom = false; end
if realrandom==false, rand('state',studnumber); randn('state',studnumber); end

x = linspace(0,1,1000);
y_f = arrayfun(@f, x);  
y_g = arrayfun(@g, x);
f_min = (1/exp(1)) * ones(1,1000);
g_min = (1/2) * ones(1,1000);
G = arrayfun(@(z) CDF(z), x);
Q = arrayfun(@(z) quantileF(z), x);
check1 = arrayfun(@(z) CDF(quantileF(z)), x);
check2 = arrayfun(@(z) quantileF(CDF(z)), x);

figure(1)
plot(x, y_f, 'r-', 'LineWidth',2 ); 
hold on 
plot(x,y_g,'b-','LineWidth',2);
plot(x,f_min, 'r--', 'LineWidth',0.1)
plot(x,g_min, 'b--', 'LineWidth',0.1)
xlabel('x');
% add 1/e to y axis
currentTicks = yticks; 
newValue = 1/exp(1); 
yticks(unique([currentTicks, newValue])); 
legend('f(x)','g(x)', 'f_{min}', 'g_{min}', 'Location', 'north')


figure(2)
plot(x, G, 'k-', 'LineWidth',2 ); 
hold on 
line([0.5 0.5], ylim, 'Color', 'black', 'LineStyle', '--');
line(xlim, [0.5 0.5], 'Color', 'black', 'LineStyle', '--');
xlabel('x');
ylabel('G_{X}')

figure(3)
plot(x, Q, 'm', 'LineWidth',2 ); 
hold on 
line([0.5 0.5], ylim, 'Color', 'black', 'LineStyle', '--');
line(xlim, [0.5 0.5], 'Color', 'black', 'LineStyle', '--');
xlabel('p');
ylabel('Q')

figure(4)
plot(x, check1, 'y', 'LineWidth',1);
hold on
plot(x,check2, 'k--')
legend('G_{X}(Q(p))', 'Q(G_{X}(x))','Location', 'north')
xlabel('x or p')
ylabel('composition')

%Parameters estimation
m=1;
n=100000;
X = randexpcos(m,n);
mu = mean(mean(X));
%Variance estimation knowing the mean
deviations = X - (1/2);
variance_app1 = zeros(m,1);
for i = 1:m
    deviations_squared = deviations(i,:).^2;
    variance_app1(i,1) = sum(deviations_squared)/(n-1);
end
variance1 = mean(variance_app1);
%Variance estimation without knowning the mean
variance_app2 = var(X,1,2);
variance2 = mean(variance_app2);



function y = f(x)
    if x >= 0 && x <= 1
        y = exp(cos(2 * pi * x));  
    else
        y = NaN; 
    end
end

function y = g(x)
    if x >= 0 && x <= 1
        y = (1 + 10 * abs(x - 0.5)) / 2;  
    else
        y = NaN;  
    end
end

function G = CDF(x)
    if x >= 0 && x <= 1/2
        G = (12/7)*x - (10/7)*x^2; 
    elseif x > 1/2 && x <= 1
        G = (10/7)*x^2 - (8/7)*x + 5/7;  
    else
        G = NaN;  
    end
end

%version that can work also with vectors
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