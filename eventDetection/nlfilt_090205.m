function yfilt = nlfilt_090205(y, Nnbrs, h)
% function yfilt = nlfilt_090205(y, Num_nghbrs, h)

%close all
DEBUG=0;
if nargin==0
    y = synthdata;  % make some fake data(below) for testing
    load ysynth
    yold = y;
    DEBUG=1;
end

y = y(:)'; % make sure it's a row vector


% h = 1;
%sigma = 3; % width of Gaussian kernel(std)
% Nnbrs = 5 ; % width of the non-local filter window


% zero pad the data
y = [zeros(1,Nnbrs) y zeros(1,Nnbrs)];

t = 1:length(y);
% compute the self-similatity between each two points (N(t1, t2))

% Gaussian smoothing:
% Gkernel = normpdf(t, 1, sigma);
% tmp = conv(Gkernel,y);
% y = tmp(1:length(y));

weights = zeros(1, length(y));
edistance = zeros(length(y));
ytmp = zeros(length(y));


for t1 = Nnbrs+1:length(y)-Nnbrs

    % compute edistance with other voxels
    for t2 = Nnbrs+1:length(y)-Nnbrs

        edistance(t1, t2) = sqrt(sum((y(t1-Nnbrs:t1+Nnbrs) - y(t2-Nnbrs : t2+Nnbrs)).^2));

    end

    % smooth the edistance
    %    tmp = conv(Gkernel, edistance(t1,:));
    %    weights = tmp(1:length(y)) ;

    weights = exp( -edistance(t1,:) / h^2 );
    % normalize weights:
    weights = weights / sum(weights);

    % apply edistance weighting to data at each time point
    ytmp(t1,:) = weights .* y;

    if DEBUG
        % some plots for debugging:
%         % plot(Gkernel,'r');
%         plot(squeeze(edistance(t1,:)),'g');
%         hold on;
%         plot(weights*100, 'b');
%         axis([1 length(y) -5 10]);
%         hold off; drawnow
    end
end

% average across all columns (time points)
yfilt = sum(ytmp,1); %./ [1:length(y)];

%take out the padding:
yfilt = yfilt(Nnbrs+1:end-Nnbrs);
y = y(Nnbrs+1 : end-Nnbrs);

yfilt = yfilt * sum(abs(y))/ sum(abs(yfilt));  % normalization
yfilt = yfilt';

if DEBUG
    subplot(211); hold off;  plot(y); hold on ;  plot(yfilt,'r'); plot(ysynth,'k');
    subplot(212); hold off; plot(abs(fft(y))); hold on; plot(abs(fft(yfilt)),'r');
    save work
end



return

%%
function y=synthdata
TP = 10 ; % this is the number of POSTIVES
s2 = 5;
Nnbrs1=6;
Nnbrs2=16;

h = inline('t.^(Nnbrs1-1).*exp(-t)/gamma(Nnbrs1) - 0.16*t.^(Nnbrs2-1).*exp(-t)/gamma(Nnbrs2)','Nnbrs1','Nnbrs2','t');
T = 200;
t = [0:2:2*T-2]';
hh = h(Nnbrs1,Nnbrs2,t);   hh = hh - mean(hh); hh = hh/norm(hh);
H = toeplitz( [hh(1); zeros(T-1,1)], hh')';  %  *** there was a transpose missing here !  I don't agree lets talk about this next meeting ***
TP = 20;
%%%%%% Define the event sequence %%%%%%%%%
x=zeros(T,1);
% note that I'm sticking in a long activation period too.
onsets = (T-15)* rand(TP,1) +1;
x(round(onsets)) = 1;
x(40:70) = 1;

n = sqrt(s2)*randn(T,1);
ysynth = H*x;
%plot(ysynth,'k'); hold on
y = ysynth+ n;
save ysynth
return