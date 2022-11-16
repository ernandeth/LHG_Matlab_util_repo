function S = sfcy(a,b)
% Sufficiency = P(B|A)

t=[sum(a&b),sum(a&~b);sum(~a&b),sum(~a&~b)];
p_y_given_x=t(1,1)/sum(t(:,1));
p_y_given_not_x=t(1,2)/sum(t(:,2));
PNS=p_y_given_x-p_y_given_not_x;
S=PNS/(1-p_y_given_not_x);

if S<0
    S=nan;
end

disp('---SUFFICIENCY------')
disp('         a          not a')
disp(['    b    ',num2str(t(1,1)),'  ',num2str(t(1,2)),'    ',num2str(sum(t(1,:)))])
disp(['not b    ',num2str(t(2,1)),'  ',num2str(t(2,2)),'    ',num2str(sum(t(2,:)))])
disp(['         ',num2str(sum(t(:,1))),'  ',num2str(sum(t(:,2))),'    ',num2str(sum(t(:)))])
disp('---------------------------')


