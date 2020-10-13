function [g1 g2 eTot e1 e2] = grads(thet,mx,cr0,nNeuron,nK)




%
% Evaluate functions
%
if nK>1
mu = 1./(thet+1);mu = prod(mu')';
else
    mu = 1./(thet+1);
end
cX = ones(nNeuron,nNeuron);
for k=1:nK
    t = thet(:,k) + thet(:,k)' + 1;cX = cX.*(1./t);
end
cX = cX - diag(diag(cX));

%
% Gradient of e1 wrt theta and lambda
%
t1 = -diag(sparse(mu.*(mu-mx)));
t2 = 1./(1+thet);
g1 = t1*t2;
e1 = mu-mx;e1 = e1'*e1;
% 


%
% Gradient of e2 wrt theta and lambda
%
t1 = cX - cr0;
for k=1:nK
    t2 = -cX./(1+thet(:,k) + thet(:,k)');
    t3 = (t1.*t2)*ones(nNeuron,1);
    g2(:,k) = t3;
end



%
% Error
%
e2 = sum(sum(t1.*t1));
eTot = e1 + e2;


end

