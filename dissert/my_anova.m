function result = my_anova
% reads data from file to matrix and performs anova calculations os it

	name = 'sadfas';

	while	name ~isempty(name)
		name = input('Select ANOVA data file: ','s');
		if ~isempty(name)
			cols = input('How many columns: ');
			reps = input('How many reps. per measurement in data set: ');
			data = read_mat(name,cols)

			sz = size(data);
			data2 = [data(:,1) data(:,2)];
			for i=2: sz(2)/2
				data2 = [data2;
					data(:,2*i-1) data(:,2*i)];
			end


			sprintf('ANOVA returns:\n')
			result=anova2(data2,reps)
		end;
	end;
return;

