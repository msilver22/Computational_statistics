
if exist('samplesize') ~= 1, samplesize = 100; end
if exist('modelsize') ~= 1, modelsize = 50; end
if exist('realrandom') ~= 1, realrandom = false; end
%if exist('studnumber') ~= 1, studnumber = 12345; end
if exist('studnumber') ~= 1, studnumber = 000601156; end
if exist('betamax') ~= 1, betamax = 10; end
if exist('maxit') ~= 1, maxit = 100; end
if exist('figsize') ~= 1, figsize = [1280 420]; end
if exist('papersize') ~= 1, papersize = figsize/96; end
if exist('figpaperpos') ~= 1, figpaperpos = [0 0 papersize]; end

forestgreen = [34 139 34]/255; fg = forestgreen;
linedotopt = {'-o','linewidth',3,'markersize',8,'color'};
linedotoptg = [linedotopt(:)',{fg},{'markerfacecolor'},{'w'}];
linedotoptr = [linedotopt(:)',{'r'},{'markerfacecolor'},{'w'}];
linedotoptb = [linedotopt(:)',{'b'},{'markerfacecolor'},{'w'}];
bulletcolor = 'b';
bullet = {'o','markersize',6,'linewidth',2,'color',bulletcolor,...
          'markerfacecolor',bulletcolor};

seed = studnumber;
if realrandom==false, randn('state',seed); rand('state',seed); end

samplesize = 1000;
n = samplesize; m = modelsize;
x = sort(rand(n,1));
p = (cos(5*x.^2)+1)/2; U = rand(n,1); Y = double(U<p); % DGP
X = ones(n,m); for k=1:m-1, X(1:n,k+1) = x.^k; end
% alternative: DGP within specified model
% beta = round(rand(m,1)*2*betamax-betamax);
% theta = X*beta; p = 1./(1+exp(-theta)); U = rand(n,1); Y = double(U<p); % DGP
theta = log(p./(1-p));
U = rand(n,1); Y = double(U<p);

% Model selection
submodels = cell(1,m); submodels{1} = 1;
for k=2:m, submodels{k} = (1:k); end
pp = (1:m);

% CpGLM = icGLM(X,Y,'bin','Cp',submodels)*stdev^2/n;
AICGLM = icGLM(X,Y,'bin','AIC',submodels);
BICGLM = icGLM(X,Y,'bin','BIC',submodels);

figure(1)
plot(pp,AICGLM,linedotoptb{:})
title('AIC'); 
xlabel('modelsize');              
ylabel('ic'); 

figure(2)
plot(pp,AICGLM,linedotoptb{:})
title('AIC'); 
xlabel('modelsize');              
ylabel('ic'); 
xlim([1 15]);
ylim([-0.7 -0.58]);

figure(3)
plot(pp,BICGLM,linedotoptb{:})
title('BIC'); 
xlabel('modelsize');              
ylabel('ic');  

figure(4)
plot(pp,BICGLM,linedotoptb{:})
title('BIC'); 
xlabel('modelsize');              
ylabel('ic');  
xlim([1 15]);
ylim([-0.7 -0.6]);

% [minCp pCp] = min(CpGLM);
[maxAIC, pAIC] = max(AICGLM);
[maxBIC, pBIC] = max(BICGLM);

% pCp
pAIC
pBIC

%IRLS with BIC
S = submodels{pBIC};
Xp = X(1:n,S);
[betahat, phat, thetahat, niter] = IRLS(Xp,Y,'bin',maxit);
% thetahat = X*betahat; phat = 1./(1+exp(-thetahat));
figure(5)
plot(x,p,'r','linewidth',3)
hold on
plot(x,phat,'m','linewidth',3)
plot(x,Y,bullet{:})
hold off 
xlabel('x');              
ylabel('p');  
legend('true model', 'IRLS estimation', 'data')
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)

%IRLS with AIC
S = submodels{pAIC};
Xp = X(1:n,S);
[betahat, phat, thetahat, niter] = IRLS(Xp,Y,'bin',maxit);
% thetahat = X*betahat; phat = 1./(1+exp(-thetahat));
figure(6)
plot(x,p,'r','linewidth',3)
hold on
plot(x,phat,'m','linewidth',3)
plot(x,Y,bullet{:})
hold off
xlabel('x');              
ylabel('p');  
legend('true model', 'IRLS estimation', 'data')
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)

