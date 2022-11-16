function events=detect_event(X)
% function events=detect_event(X)
% 
% X: data
% events: binary vector (1 active, 0 non active)

X = X(:);
Npts=length(X);

%{
% lambda: smooting parameter
 baseline_data = baseline_est_EM(X);
mu = mean(baseline_data); % mean of the baseline 
Z1 = zeros(Npts,1);
Z1(1) = mu;
for j=2:Npts
    Z1(j) = lambda*X(j) + (1-lambda)*Z1(j-1); 
end
%}

th1 = MixModel2EM(X);
events = X>=th1;

function th = MixModel2EM(data)

c=2;
n=length(data);
% inicializacion usando FCM.
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
th=0.5*(smu(1)+smu(2));
% th=smu(1);

function baseline_data = baseline_est_EM(data)

c=2;
n=length(data);
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
