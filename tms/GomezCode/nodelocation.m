function nodeijk=nodelocation(i,j,k)
global dx dy dz;
nodeijk(1)=dx*(i-1);
nodeijk(2)=dy*(j-1);
nodeijk(3)=dz*(k-1);
