
function result = remove(data, mini, maxi)
% function result = remove(data, mini, maxi)
% Remove data points whose distance to the center are out of range
	sz = size(data);
	rows = sz(1);
	cols = sz(2);
	j=1;

	for i=1:rows,
      if (data(i,1) < maxi)  & ( data(i,1) > mini) ...
            & (data(i,2) < maxi)  & ( data(i,2) > mini) ...
            & (data(i,3) < maxi)  & ( data(i,3) > mini)

			result(j,:) = data(i,:);
			j=j+1;
		end;
	end;
return
