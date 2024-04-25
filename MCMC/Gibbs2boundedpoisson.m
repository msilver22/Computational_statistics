% This simulation study allows to set three parameters, for instance, run
%   >> nsamples = 1000; upperbound = Inf; mumin = 0; Gibbs2boundedpoisson
% Use the command clear to remove all current variables from workspace
% nsamples is the number of Gibbs samples
% upperbound is the maximum value in each observation (the maximum count)
% mumin is the minimum expected value, to be added to the observed

if exist('nsamples') ~= 1, nsamples = 10000; end
if exist('upperbound') ~= 1, upperbound = 100; end
if exist('mumin') ~= 1, mumin = 0.5; end
figsize = [1200 480];
papersize = figsize/96; figpaperpos = [0 0 papersize]; figpos = [0 0 figsize];


b = floor(upperbound); a = mumin;
%or
%a = 0;
%b = Inf;

X1 = 1;
xx1 = zeros(1,nsamples); xx2 = zeros(1,nsamples);
rand('state',601156);
for s=1:nsamples
   X2 = min(b,poissonnoise(X1+a,true));
   X1 = min(b,poissonnoise(X2+a,true));
   xx1(s) = X1; xx2(s) = X2;
end

figure(1)
%for the case b = Inf
vmax = max(max(xx1),max(xx2));
xmax = min(b,vmax);
%xmax = min(b,100);
x = (0:xmax);
pX1hat = zeros(size(x)); pX2hat = zeros(size(x));
for k=0:xmax, pX1hat(k+1) = sum(xx1==k); pX2hat(k+1) = sum(xx2==k); end
pX1hat = pX1hat/nsamples; pX2hat = pX2hat/nsamples;
FX1hat = cumsum(pX1hat); FX2hat = cumsum(pX2hat);
%bar(x,FX1hat,'b')
plot(x,FX1hat,'g--','LineWidth',3)
hold on
plot(x,FX2hat,'linewidth',1,'color','r')
% bar(x,FX2hat,0.5,'r')
title('Empirical CDF for X1 and X2')
xlabel('values')
legend('F_{X_{1}} hat', 'F_{X_{2}} hat','Location','north')

figure(2)
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)
bar(x,pX1hat,0.7,'k')
title('Empirical PMF X1')
xlabel('values')

figure(3)
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)
bar(x,pX2hat,0.7,'k')
title('Empirical PMF X2')
xlabel('values')

figure(4)
%plot(sort(xx1),(1:nsamples)/nsamples,'b')
plot(xx1(1:30),'b-', 'LineWidth',1)
title('X_{1}')
xlabel('samples')
ylabel('X_{1}')
%hold on
%plot(sort(xx2),(1:nsamples)/nsamples,'r--')
%hold off

figure(5)
%plot(sort(xx1),(1:nsamples)/nsamples,'b')
plot(xx2(1:30),'b-', 'LineWidth',1)
title('X_{2}')
xlabel('samples')
ylabel('X_{2}')

figure(6)
%plot(sort(xx1),(1:nsamples)/nsamples,'b')
plot(xx2(1:length(xx2)),'b-', 'LineWidth',1)
title('X_{2}')
xlabel('samples')
ylabel('X_{2}')


