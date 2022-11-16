function sigmap=sigma3(n)
load('/home/luisgo/sigma2.mat');
sigmap(ceil(length(sigma(:,1,1))/n),ceil(length(sigma(1,:,1))/n),length(sigma(1,1,:)))=0;
for i=1:length(sigma(:,1,1))
    for j=1:length(sigma(1,:,1))
        for k=1:length(sigma(1,1,:))
            sigmap(ceil(i/n),ceil(j/n),ceil(k/n))=sigmap(ceil(i/n),ceil(j/n),ceil(k/n))+sigma(i,j,k)/n^3;
        end
    end
end
