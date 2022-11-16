% running regression
p=3;
N=vnum-2;
doNuisance = 1;
twindow = 100;

if vnum>twindow+1
    nn = N  - twindow+1;
else
    nn=1;
end

% override to use the whole data for analysis
% nn=1;

ref1 = ref(nn:N);
%ref1 = ref1-mean(ref1);


X = [ref1  ones(N-nn+1,1)];
C = [ 1 0];

if doNuisance
    X = [ref1 junkreg(nn:N) ones(N-nn+1,1)];
    C = [ 1 0 0];
end


xtx_inv = pinv(X);
beta_est = xtx_inv*allASL(nn:N, :) ;

cx = C * xtx_inv;
cxt = xtx_inv' * C';

residual(nn:N, :) = allASL(nn:N,:) - X*beta_est(:,:);
clean(nn:N, :) = allASL(nn:N,:) - X(:,2)*beta_est(2,:);



% this was actually slower!

% % compute RSS
% inds = find(slmask);
% for pix=inds
% %         residual(nn:N, pix) = allASL(nn:N,pix) - X*beta_est(:,pix);
% %         clean(nn:N, pix) = allASL(nn:N,pix) - X(nn:N,2)*beta_est(2,pix);
%         var_est(pix) = (residual(nn:N,pix))' * (residual(nn:N,pix)) / (N-p);
%         var_con(pix) = cx * var_est(pix) * cxt;
% end

var_est = sum(residual.^2, 1) / ((N-nn)-p);
var_con = (cx*cxt) * var_est;

% update T score map
corrMap = beta_est(1,:) ./ sqrt(var_con);
