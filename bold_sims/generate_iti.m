mean_iti = 15
DURATION = 600



iti = (mean_iti/2)*randn(size(times));
iti = iti + abs(min(iti)) + 1;
iti = iti*1000
save iti.dat iti -ASCII

iti(2:end) = iti(2:end) +500 
times = cumsum(iti);

save times.dat times -ASCII


!unix2dos times.dat times_dos.dat
!unix2dos iti.dat iti_dos.dat

