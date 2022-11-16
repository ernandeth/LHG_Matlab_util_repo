for n=1:2:size(raw,2)-1
    subs = raw(:,n+1) -raw(:,n);
end
msubs = mean(subs,2);

plot(abs(fft(msubs)));
