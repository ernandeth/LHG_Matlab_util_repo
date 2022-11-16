function out = doaliasing(x1 , N, R)

FF = FourierMat(N);
% downsample the encoding matrix
Fdown = zeros((N^2)/R, N^2);
for n=1:N/R
    start1= (n-1)*N +1; fin1 = (n-1)*N +N;
    start2= (n-1)*N*R +1;  fin2 = (n-1)*N*R +N;
    Fdown(start1 :fin1, :) = FF( start2 : fin2, :);
end

x2 = Fdown*x1;
x3 = [zeros(N,N/R)  x2  zeros(N,N/R)];
out = FF * x3;

return

