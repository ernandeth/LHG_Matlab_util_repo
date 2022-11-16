function [N,S]=nec_suf(x,y)

[sx,x] = proc_signal(x);
[sy,y] = proc_signal(y);

% x and y are signals of events

t = [sum(x&y),sum(~x&y);
    sum(x&~y),sum(~x&~y)];

sumrow=sum(t);

p_y_given_x = t(1,1)/sumrow(1);
p_y_given_notx = t(1,2)/sumrow(2);
p_noty_given_notx = t(2,2)/sumrow(2);


PNS = p_y_given_x-p_y_given_notx;
bound1 = max(0,p_y_given_x-p_y_given_notx);
bound2 = min(p_y_given_x,p_noty_given_notx);

% disp('-----------------')
% disp(['p(y|x)=',num2str(p_y_given_x)])
% disp(['p(y|~x)=',num2str(p_y_given_notx)])


if PNS<bound1
    PNS=0;
    N=0;
    S=0;
    
%     disp(['PN=0'])
%     disp(['PS=0'])
%     disp('-----------------')
    return
elseif PNS>bound2
    PNS=1;
    N=1;
    S=1;
    
%     disp(['PN=0'])
%     disp(['PS=0'])
    return
end

N = PNS/p_y_given_x;
S = PNS/(1-p_y_given_notx);

if isnan(N) | isinf(N)
    N=0;
end

if isnan(S) | isinf(S)
    S=0;
end

% disp('-----------------')
% disp(['PN=',num2str(N)])
% disp(['PS=',num2str(S)])
% disp('-----------------')

function [y,ey] = proc_signal(y)

y= y-mean(y);
y = detrendfcn(y,0);
y = mrtemporalfilter(y);
y= y-mean(y);    
[ey,th] = event_detection(y,1);

function y = mrtemporalfilter(y)

sigma = 2;
size = 10;
x = linspace(-size/2,size/2,size);
mask = exp(-x.^2/(2*sigma^2));
mask = mask/sum(mask); 
y = filtfilt(mask,1,y);

function y=detrendfcn(x,doplot)

axx = (1:length(x))';
p = polyfit(axx,x,3);
yh = polyval(p,axx);
y = x-yh;

if doplot
plot(x,'g')
hold on
plot(yh,'g')
hold on
plot(y,'r')
end
