tdel=ones(20,6);
for n=1:6
    tdel(:,n) = 0.2 + n*0.3;
end

tdel = tdel(:);

ttag = tdel;
ttag(:) = 2;

tadj = 6-ttag-tdel - 0.035;

save t_tags.txt ttag -ascii
save t_delays.txt tdel -ascii
save t_adjusts.txt tadj -ascii


