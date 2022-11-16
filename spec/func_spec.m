function fdata = func_spec(raw_data)
% Fit and subtract the water signal:
water_data = [];
for i=1:size(raw_data,2)
	fprintf('\rFID N. %d', i);
	guess =lsqcurvefit('t2_decay', [16 -0.01],[1:2048]', abs(raw_data(:,i) ));
	m0 = guess(1);
	R2 = guess(2);
	water = t2_decay([m0 R2],  [1:2048]);
	water_data = [water_data ; water];

end
whos
data = abs(raw_data)- abs(water_data');
fdata=fftshift(fft(data')))';

return


