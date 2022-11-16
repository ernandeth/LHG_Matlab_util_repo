
% GRAPH.M 
% Procedure to plot however many files containing three columns

cols = 3

name = 'junk';
while ~isempty(name),
	name = input('select file: ','s');
	if ~isempty(name),
		hold on
		grid on
		mat = read_mat(name,cols);
		plot(mat(:,1), mat(:,2),'b' );
		plot(mat(:,1), mat(:,3),'r' );
		axis ([-1 1 -1 1]);
	end
end
