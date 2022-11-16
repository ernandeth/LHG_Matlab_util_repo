b=[0 21.44 85.7 193 343 536 772 1051]; %(s/mm2)
!ls -1 P* > files.txt
fid=fopen('files.txt','r');

shiftdrift = 0.16;

%    water   mI    Cho   Cr    NAA
T2 = [70     270   350   380   305]     /1000;
wo = [4.7    3.59  3.21  3.04  2.02]; 

T2 = T2 * 5000/2048;

wo = [4.7 - wo + shiftdrift] * 128 * 2048/5000;
wo(1) = 0;

parms = [T2'  wo'];

raw_data=[];
cleandata = [];

for i=1:size(b,2)

	name=fgetl(fid);
	raw=mean(read_cplx( name,2048,16),2);
	%raw_data=[raw_data; raw'];

	guess =lsqcurvefit('t2_decay', [16 0.001],[1:2048]', abs(raw));
	m0 = guess(1);
	R2 = guess(2);
	water = t2_decay([m0 R2],  [1:2048]);
    data = abs(raw')- abs(water);
    fdata=(abs(fft(data')))';
    fdata = fdata(1:1024);

    
    
    Mofit = lsqcurvefit('NMRfunc',[20 1 2 3 4 ], [1:1024], fdata, [],[],[],parms); 
    cdat = NMRspect(Mofit, wo, T2, [1:1024]);
    cleandata = [cleandata; cdat];
    if i==1
        clf
        plot(fdata)
        hold on
        plot(cdat, 'r');
        axis([50 170 0 1000])
        pause
    end
end


fclose(fid);


slopes = [];
for i=1:size(cleandata,2)
		
	sl= - regression(b', log(cleandata(:,i)));
	slopes=[slopes sl];

end

subplot 211, plot(slopes),axis([0 300 0 0.003]), ...
    title('Diffusion coefficent spectrum (mm2/s)');
subplot 212, plot(cleandata(1,:)), axis([0 300 0 700]),...
    title('NMR Spectrum');


