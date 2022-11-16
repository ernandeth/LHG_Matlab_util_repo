function N = ncsty(a,b)
% wrong (luis's) deffinition
% Necessity = 1-P(B|~A);
%  = P(~B | ~A)
%
% p_ab = sum(a.*b)/length(a)
% p_bnota = sum( (~a) .*b ) / length(a) ;
% p_nota = sum(~a) / length(a);
% N =  1 - p_bnota / p_nota;
%
% Pearl's defin.
% PNS = p_y_given_x - p_y_given_not_x;
%      --------------------------------
%         p_y_given_x;

t = [sum(a&b),sum(a&~b);
    sum(~a&b),sum(~a&~b)];
p_y_given_x = t(1,1)/sum(t(:,1));
p_y_given_not_x = t(1,2)/sum(t(:,2));
PNS = p_y_given_x - p_y_given_not_x;
N = PNS/p_y_given_x;

% 
% if N<0
%     N=nan;
% end
% 
% disp('---NECESSITY------')
% disp('         a          not a')
% disp(['    b    ',num2str(t(1,1)),'  ',num2str(t(1,2)),'    ',num2str(sum(t(1,:)))])
% disp(['not b    ',num2str(t(2,1)),'  ',num2str(t(2,2)),'    ',num2str(sum(t(2,:)))])
% disp(['         ',num2str(sum(t(:,1))),'  ',num2str(sum(t(:,2))),'    ',num2str(sum(t(:)))])
% disp('---------------------------')
