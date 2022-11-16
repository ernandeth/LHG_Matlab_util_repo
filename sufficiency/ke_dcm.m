clear;
load sim5;
%dz/dt = sigma*A*z + C*u + N(0, 1/20), A is the connection matrix, C is identity matrix I, z is time series of each node, u is external inputs, sigma controls within node temporal
%smoothing and neural lag between nodes;
%use a poisson process to control the state(0 or 1) in u(t), the external input
A = reshape(net(22, :, :), Nnodes, Nnodes)'; C = eye(Nnodes); M = 1;
x = zeros(Nnodes, 1); 
% control = poissrnd(0.2, [1000 1]);

uold = poissrnd(0.2, [Nnodes 1]); count = zeros(Nnodes, 1);
while (M < 200)
%     control = poissrnd(0.2, 1);
%       u = double(poissrnd(0.2, [Nnodes 1]) & 1);
%       xdot = A * x(:, end) + C * u + normrnd(0, 1/20, [Nnodes 1]);
%       x = [x (x(:, end) + xdot)];
%       M = M+1;
    



%**************************************************************************************************
    unew = zeros(Nnodes, 1);
    
    %change excitation u every time, make sure u's independent, and silent time should be 4 times than the active time, as in the paper's realization
    for i = 1 : Nnodes
        if (uold(i) == 0) %silent state
            
            if (count(i) < 4) %silent state lasts 4 samples
                unew(i) = 0; count(i) = count(i) + 1;
            else
                unew(i) = poissrnd(0.2, 1); count(i) = 0;
            end
            
        elseif (uold(i) ~= 0) %active state lasts 1 sample
            
            unew(i) = double(poissrnd(0.2, 1) & 1);
        end
           
    end
    
    
    xdot = A * x(:, end) + C * unew + normrnd(0, 1/20, [Nnodes 1]);
    x = [x (x(:, end) + xdot)];
    M = M+1;
    uold = unew;
end

xdata = x';
save('DCM', 'x');
%See what x is
haemody = [ 0 0 1 5 8 9.2 9 7 4 2 0 -1 -1 -0.8 -0.7 -0.5 -0.3 -0.1 0 ];
xconv = zeros(size(temp, 1) + 18, Nnodes);
for j = 1 : Nnodes
    xconv(:, j) = conv(temp(:, j), haemody);
end

xdata = xconv;
figure; plot(xdata(:, 1), '--r'); hold on; plot(xdata(:, 2)); legend('Series 2', 'Series 5'); hold off;
% figure; plot(xconv(:, 1), '--r'); hold on; plot(xconv(:, 2)); legend('Series 1', 'Series 2'); hold off;

% xdata = xconv;
%Apply threshold method
xdata = xdata - min(min(xdata));
Nmat = zeros(Nnodes, Nnodes); Smat = zeros(Nnodes, Nnodes); 
sigma = (std(xdata,[],1))';

thresh = zeros(size(sigma, 1), 9);
for i = 1 : 9
    thresh(:, i) = (1+0.5*(i-1))*sigma;
end

Nmat1 = zeros(Nnodes, Nnodes); Smat1 = zeros(Nnodes, Nnodes);
i = 2;
% for i = 1 : size(thresh, 2)  % try different thresholds
    
%     if thresh(i) == 2*sigma
        temp = zeros(size(xdata, 1), Nnodes);
        ts_new = zeros(size(xdata, 1), Nnodes); bt = ones(size(xdata, 1), 1);
        for j = 1 : Nnodes
            ts_new(:, j) = bt .* (xdata(:, j) > thresh(j, i));  %binarize the data
            temp(:, j) = ts_new(:, j);
        end

        for k = 1 : Nnodes-1   % calculate N/S value
            %for w = 1 : Nnodes-k
            w = k+1
                Nmat(k, w) = 1 - sum(ts_new(:, k) .* (1-ts_new(:, w)))./sum(1-ts_new(:, w));
                Smat(k, w) = sum(ts_new(:, k) .* ts_new(:, w)) ./ sum(ts_new(:, k));
            %end
        end

        %normalize the values
        for k = 1 : Nnodes-1
            %for w = 1 : Nnodes-k
            w = k+1;
                Nmat1(k, w) = (Nmat(k, w) - Nmat(w, k)) ./ (Nmat(k, w) + Nmat(w, k));
                Smat1(k, w) = (Smat(k, w) - Smat(w, k)) ./ (Smat(k, w) + Smat(w, k));
            %end
        end
    
    
        figure; imagesc(Nmat1); colormap gray; colorbar; title(['Nmat normalized, threshold = sigma * ', num2str(1+0.5*(i-1))]);
        figure; imagesc(Smat1); colormap gray; colorbar; title(['Smat normalized, threshold = sigma * ', num2str(1+0.5*(i-1))]);
%     end

% end
