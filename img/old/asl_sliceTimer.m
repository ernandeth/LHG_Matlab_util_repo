function asl_sliceTimer(root, TR, Ttag)
% interpolates the data to the END of the TR
% Assumes The first slice is acquired at TR+Tag
%

[data h] = read_nii_img(root);
Nslices = h.dim(4);
Nframes = size(data,1);
Npix = h.dim(2) * h.dim(3);
Tslice = (TR-Ttag)/Nslices;
TRtimes = TR * [1:Nframes];
outdata = zeros(size(data));

for sl = 1:Nslices
	fprintf('\rslice timing : ... %d', sl);
	AQtimes= TR*(0:Nframes-1) + Ttag + (sl-1)*Tslice;
	for p=1:Npix
		ts = data(:, Npix*(sl-1) + p);	
		outTS = interp1(AQtimes, ts, TRtimes,'sinc');
		outdata(:, Npix*(sl-1) + p) = outTS;	
	end

end
warning off

write_nii(['a' root], outdata, h, 0)
%keyboard
warning on
return
