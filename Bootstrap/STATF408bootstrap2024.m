screensize = get(0,'screensize');
figwidth0 = screensize(3); figsize = [figwidth0 420];
papersize = figsize/96; figpaperpos = [0 0 papersize];
if exist('samplesize') ~= 1, samplesize = 20; end
studnumber = 601156;
if exist('realrandom') ~= 1, realrandom = false; end
if realrandom==false, rand('state',studnumber); randn('state',studnumber); end

samplesize = 2000;
B = 10000; %number of resamples
n = samplesize; % samplesize
mu = 3; theta = exp(mu); %true parameter
ns = 1000; % number of samples
basic = zeros(ns,2);
percentile = zeros(ns,2);
T = zeros(ns,1);
tstaralfa = zeros(ns,2);
alfa = 0.05;
qq = round([alfa/2 1-alfa/2]*B);
for s=1:ns
   X = exp(randn(n,1)+mu); %DGP
   thetahat = median(X); 
   T(s) = thetahat;
   r = ceil(rand(n,B)*n);
   Xstar = X(r);
   thetahatstar = sort(median(Xstar));
   tstaralfa(s,1:2) = thetahatstar(qq);
end
   
basic = 2*T-tstaralfa(1:ns,[2 1]);
percentile = tstaralfa;

figure(1)
[i1, b1] = plotasblocks(basic(1:100,1));
[i2, b2] = plotasblocks(basic(1:100,2));
nn = length(i2);
ii = [i1; i2(nn:-1:1)];
bb = [b1; b2(nn:-1:1)];
fill(ii,bb,'y')
hold on
plotasblocks(T(1:100),'k','linewidth',2)
plot([1 100],[theta theta],'r','linewidth',2)
hold off
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)
legend('Basic IC', 'theta hat', 'True median')
title('Basic bootstrap CIs for the first 100 samples')
xlabel('samples')

figure(2)
[i1, b1] = plotasblocks(basic(1:ns,1));
[i2, b2] = plotasblocks(basic(1:ns,2));
nn = length(i2);
ii = [i1; i2(nn:-1:1)];
bb = [b1; b2(nn:-1:1)];
fill(ii,bb,'y')
hold on
plotasblocks(T,'k')
plot([1 ns],[theta theta],'r','linewidth',2)
hold off
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)
legend('Basic IC', 'theta hat','True median')
title('Basic bootstrap confidence intervals')
xlabel('samples')

figure(3)
[i1, b1] = plotasblocks(percentile(1:ns,1));
[i2, b2] = plotasblocks(percentile(1:ns,2));
nn = length(i2);
ii = [i1; i2(nn:-1:1)];
bb = [b1; b2(nn:-1:1)];
fill(ii,bb,[1 0.7 1])
hold on
plotasblocks(T,'k')
plot([1 ns],[theta theta],'r','linewidth',2)
hold off
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)
legend('Percentile IC','theta hat', 'True median')
title('Percentile bootstrap confidence intervals')
xlabel('samples')

figure(4)
[i1, b1] = plotasblocks(percentile(1:100,1));
[i2, b2] = plotasblocks(percentile(1:100,2));
nn = length(i2);
ii = [i1; i2(nn:-1:1)];
bb = [b1; b2(nn:-1:1)];
fill(ii,bb,[1 0.7 1])
hold on
plotasblocks(T(1:100),'k','linewidth',2)
plot([1 100],[theta theta],'r','linewidth',2)
hold off
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)
legend('Percentile IC', 'theta hat', 'True median')
title('Percentile bootstrap CIs for the first 100 samples')
xlabel('samples')

coverbasic = sum((basic(1:ns,1)<theta)&(theta<basic(1:ns,2)))/ns
coverpercentile = sum((percentile(1:ns,1)<theta)&(theta<percentile(1:ns,2)))/ns
