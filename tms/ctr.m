function ctr = stim_ctr()
% Computes the center of the stimulation area by calculating the 
% of the points involved weighted average

	name =  input('Select file to analyze: ','s');
	data = read_mat(name, 4);
	sz = size(data);
	ctr = zeros(1,3);
	total_weight = sum(data(:,4));
data(:,4)
	for i=1:sz(1),
		weight = data(i,4) / total_weight;
		ctr(1,1) = ctr(1,1) + data(i,1) * weight;
		ctr(1,2) = ctr(1,2) + data(i,2) * weight;
		ctr(1,3) = ctr(1,3) + data(i,3) * weight;
	end;

	%ctr = ctr/sz(1);
	plot3(ctr(3), ctr(1), ctr(2),'*r');
return;
