x = D(:,2);
for i=1:5
	for j=1:5
		for k=1:4
			[ i j k  ]
			ThePix = sub2ind([xdim, ydim, zdim],i, j, k);
			tc = allY(:,ThePix);

			tc = detrend(tc);
			tc = tc - mean(tc);
			tc = tc/max(tc);

			plot(tc); hold on; stem(x,'r'); hold off
			pause
		end; end; end;