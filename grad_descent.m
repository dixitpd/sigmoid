clc
clear
%
load microbiome
zrx = find(mean(microbiome')==0);
microbiome(zrx,:) = [];
microbiome = 1.0*(microbiome > 0);
mx = mean(microbiome')';
cr0 =  (microbiome*microbiome')/size(microbiome,2);cr0 = cr0 - diag(diag(cr0));
nBact = size(microbiome,1);
%
nK    = 8;
lthet = 1.1*log(repmat((1./mx.^(1/nK) - 1),1,nK).*(1 + 0.5*(rand(nBact,nK))));
thet = exp(lthet);
filen = strcat('data_thet_k',num2str(nK));
load(filen)
% % %
mu = 1./(thet+1);mu = prod(mu')';
cX = ones(nBact,nBact);
for k=1:nK
    t = thet(:,k) + thet(:,k)' + 1;cX = cX.*(1./t);
end
cX = cX - diag(diag(cX));
%
for iter=1:10000
    thet = exp(lthet);
    [g1 g2 eTot e1 e2] = grads(thet,mx,cr0,nBact,nK);
    grd = (g1+g2).*thet;
    lthet = lthet -3*(grd);
    erx(iter) = eTot;
    if mod(iter,100) == 0
        [iter/100 log10([norm(grd) eTot e1 e2])]
    end
end
%
thet = exp(lthet);
save(filen,'lthet','eTot')
% subplot(1,3,1)
% hold on
% plot(erx/erx(1))
%
mX = 1./(thet+1);mX = prod(mX')';
cX = ones(nBact,nBact);
for k=1:nK
    t = thet(:,k) + thet(:,k)' + 1;cX = cX.*(1./t);
end
cX = cX - diag(diag(cX));
%
%plot([1:100:2000],ones(20,1)*0.995,'k--')
subplot(1,2,1)
loglog(mx,mX,'ko')
hold on
plot([1e-4 1],[1e-4 1],'r--')
%
subplot(1,2,2)
c1 = reshape(cr0,nBact*nBact,1);
c2 = reshape(cX,nBact*nBact,1);
loglog(c1,c2,'ko')
hold on
plot([1e-6 1],[1e-6 1],'r--')

