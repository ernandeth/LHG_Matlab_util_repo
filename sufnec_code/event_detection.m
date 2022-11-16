function [events th] = event_detection(X,lambda)

Npts = length(X);

baseline_data = baseline_est_EM(X);
mu = mean(baseline_data); % mean of the baseline 
    
Z = zeros(Npts,1);
Z(1) = mu;

for j=2:Npts % EWMA statistic
    Z(j) = lambda*X(j) + (1-lambda)*Z(j-1); 
end

th = compute_threshold_EM(Z);
events = Z>=th;

function th = compute_threshold_EM(data)

c=2;
n=length(data);
% FCM
[notused,U]=fcm(data,c,[2 100 1e-5 0]);
[notused,label]=max(U);
for i=1:c
    mu(i)=mean(data(find(label==i)));
    sigma(i)=var(data(find(label==i)));
    w(i)=sum(label==i)/n; 
end

iter=0;
max_it=200;
while iter<max_it
    for i=1:c
        h(:,i)=w(i)*gaussian_pdf(data,mu(i),sigma(i));
    end
    sumh=sum(h,2);
    h=h./repmat(sumh,1,c);
    w=sum(h)/n;
    for i = 1:c
        arg=(data-mu(i)).^2;
        tmp=sum(h(:,i));
        sigma(i)=arg'*h(:,i)/tmp;
    end
    mu=(data'*h)./sum(h);
    iter = iter+1;
    logL(iter)=sum(log(sumh));
    if iter>1 
        if logL(iter)-logL(iter-1)<1e-2
            logLout=logL(iter);
            break;
        end
    end
end

[smu,ind] = sort(mu);
th=(smu(1)+smu(2))/2;

function baseline_data = baseline_est_EM(data)

c=2;
n=length(data);
% FCM
[notused,U]=fcm(data,c,[2 100 1e-5 0]);
[notused,label]=max(U);
for i=1:c
    mu(i)=mean(data(find(label==i)));
    sigma(i)=var(data(find(label==i)));
    w(i)=sum(label==i)/n; 
end

iter=0;
max_it=200;
while iter<max_it
    for i=1:c
        h(:,i)=w(i)*gaussian_pdf(data,mu(i),sigma(i));
    end
    sumh=sum(h,2);
    h=h./repmat(sumh,1,c);
    w=sum(h)/n;
    for i = 1:c
        arg=(data-mu(i)).^2;
        tmp=sum(h(:,i));
        sigma(i)=arg'*h(:,i)/tmp;
    end
    mu=(data'*h)./sum(h);
    iter = iter+1;
    logL(iter)=sum(log(sumh));
    if iter>1 
        if logL(iter)-logL(iter-1)<1e-2
            logLout=logL(iter);
            break;
        end
    end
end

[smu,ind] = sort(mu);
baseline_data = data(data<=smu(1));

function f=gaussian_pdf(x,mu,sigma)
f=exp(-(x-mu).^2/(2*sigma))/sqrt(sigma*2*pi);