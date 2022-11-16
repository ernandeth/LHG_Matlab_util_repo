function result = sse_fun(A)
% Definition of function to be minimized in least squares algorithm
	global to;
	global from;

	result = (to - A*from);
	result = result.^2;

return;
