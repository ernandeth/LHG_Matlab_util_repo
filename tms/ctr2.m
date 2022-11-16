function ctr = stim_ctr(data)
% Computes the center of the stimulation area by calculating the 
% of the points involved weighted average


	sz = size(data);
	ctr = zeros(1,3);
	total_weight = sum(data(:,4));

	for i=1:sz(1),
		if total_weight > 0 ,
			weight = data(i,4) / total_weight;
			else weight=0;
		end;
		ctr(1,1) = ctr(1,1) + data(i,1) * weight;
		ctr(1,2) = ctr(1,2) + data(i,2) * weight;
		ctr(1,3) = ctr(1,3) + data(i,3) * weight;
	end;

	center_of_Stimulation = ctr
	
	plot3(ctr(3), ctr(1), ctr(2),'*r');

return;
