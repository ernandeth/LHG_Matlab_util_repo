b=[0 21.44 85.7 193 343 536 772 1051];
!ls -1 P* > files.txt
fid=fopen('files.txt','r');

raw_data=[];
water_data = [];

for i=1:size(b,2)

	name=fgetl(fid);
	raw=mean(read_cplx( name,2048,16),2);
	raw_data=[raw_data; raw'];

	guess =lsqcurvefit('t2_decay', [16 -0.01],[1:2048]', abs(raw));
	m0 = guess(1);
	R2 = guess(2);
	water = t2_decay([m0 R2],  [1:2048]);
	water_data = [water_data ; water];

end

data = abs(raw_data)- abs(water_data);
fdata=(abs(fft(data')))';

fclose(fid);

slopes=[];
for i=1:size(fdata,2)
		
	sl=regression(b', log(fdata(:,i)));
	slopes=[slopes sl];

end



