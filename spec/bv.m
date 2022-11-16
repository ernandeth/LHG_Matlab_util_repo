tmp = read_cplx('P24576.7',2048,16);
tmp2 = abs(fft(mean(tmp')));
a = tmp2;

tmp = read_cplx('P25088.7',2048,16);
tmp2 = abs(fft(mean(tmp')));
a = [a ; tmp2];

tmp = read_cplx('P25600.7',2048,16);
tmp2 = abs(fft(mean(tmp')));
a = [a ; tmp2];

tmp = read_cplx('P26112.7',2048,16);
tmp2 = abs(fft(mean(tmp')));
a = [a; tmp2];

