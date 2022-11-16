function result = translation(data)

% result = translation(data)
%
% Performs the translation by subtraction of the position cordinates from reference
% point coordinates
% Both Input and output are six-column matrices
% Last edit: 1-5-98
% Luis Hernandez

	sz = size(data);
	rows = sz(1);
	cols = sz(2);
	result = zeros(rows/2,cols);
	for i= 1:rows/2,
		result(i,1:3) = data(i*2-1, 1:3) - data(i*2 , 1:3);
		result(i,4:6) = data(i*2,4:6);
	end

return

