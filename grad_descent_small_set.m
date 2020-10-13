clc
clear
%
load ../data_neurons
nNeuron = 15;
[a b] = sort(mx,'descend');
goodIDs = b(1:nNeuron);
cr0 = crn - diag(diag(crn));cr0 = cr0(goodIDs,goodIDs);mx = mx(goodIDs);
nK    = 8;eO = 1e10;
%
lthet = log(rand(nNeuron,nK));
thet = exp(lthet);
filen = strcat('data_theta_k',num2str(nK));
load(filen)
% % %
mu = 1./(thet+1);mu = prod(mu')';
cX = ones(nNeuron,nNeuron);
for k=1:nK
    t = thet(:,k) + thet(:,k)' + 1;cX = cX.*(1./t);
end
cX = cX - diag(diag(cX));
%
for iter=1:500000
    thet = exp(lthet);
    [g1 g2 eTot e1 e2] = grads(thet,mx,cr0,nNeuron,nK);
    grd = (g1+g2).*thet;
    lthet = lthet - 50*grd;
    erx(iter) = eTot;
    if mod(iter,10000) == 0
        [iter/10000 log10([norm(grd) eTot e1 e2])]
    end
end
%
thet = exp(lthet);
save(filen,'lthet','eTot')
subplot(1,3,1)
hold on
plot(erx/erx(1))
%
mX = 1./(thet+1);mX = prod(mX')';
cX = ones(nNeuron,nNeuron);
for k=1:nK
    t = thet(:,k) + thet(:,k)' + 1;cX = cX.*(1./t);
end
cX = cX - diag(diag(cX));
%
plot([1:100:2000],ones(20,1)*0.995,'k--')
subplot(1,3,2)
loglog(mx,mX,'ko')
hold on
plot([1e-3 1],[1e-3 1],'r--')
%
subplot(1,3,3)
loglog(cr0,cX,'ko')
hold on
plot([1e-5 1e-1],[1e-5 1e-1],'r--')

