function outdata = slicesmooth( data , weight)

paddata = [zeros(size(data,1),1) data zeros(size(data,2),1)];
paddata = [zeros(1,size(paddata,2)); paddata ; zeros(1,size(paddata,2))];

outdata = data ...
    + weight *paddata(3:end, 2:end-1) + ...
    + weight *paddata(1:end-2, 2:end-1) + ...
    + weight *paddata(2:end-1, 1:end-2) + ...
    + weight *paddata(2:end-1, 3:end);


outdata = outdata /(4*weight+1) ;


return
