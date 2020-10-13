clc
clear
%
nNeuron = 15;
load ../data_neurons
[a b] = sort(mx,'descend');
goodIDs = b(1:nNeuron);
crn = crn(goodIDs,goodIDs);
mx  = mx(goodIDs);
%
%
nSt = 2^nNeuron;
% for i=0:nSt-1
%     ix = de2bi(i);ix = fliplr(ix);
%     l = length(ix);
%     if l < nNeuron
%         ix = [zeros(1,nNeuron-l) ix];
%     end
%     stx(i+1,:) = ix;
% end
% 'Done with the states'
% 
% J = sparse(nNeuron,nNeuron);
% 'Initialized J'
% %
load Coupling_15
for iter=1:100
    s1 = zeros(nNeuron,nNeuron);s2 = 0;
    for i=1:nSt
        px = exp(-stx(i,:)*J*stx(i,:)');
        s1 = s1 + px*(stx(i,:)'*stx(i,:));
        s2 = s2 + px;
    end
    pr = s1/s2;
    
    gr = crn-pr;
    J = J - (gr./crn)*0.05;
    if mod(iter,10) == 0
        norm((crn-pr)./crn)
    end
end
save Coupling_15 J stx
loglog(crn,pr,'ko')
hold on
plot([1e-3 1],[1e-3 1],'r--')

