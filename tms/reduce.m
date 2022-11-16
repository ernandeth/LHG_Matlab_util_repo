
function result = reduce(data, ntrvl)

% function result = reduce(data)
% Take every nth data point and throw out the rest.
% (Keep the fourth point out of the set of four)
	sz = size(data)

   
   for i=1 : sz(1)/ntrvl
      result(i, : ) = data( ntrvl*i , : ) ;
	end;
   
   return
