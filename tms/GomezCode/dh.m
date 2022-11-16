
function h=dh(r1,r2,rs,js,dir)
h=0;




for source=1:length(rs)


	  resp=squeeze(rs{source}(2:end,:)+rs{source}(1:end-1,:))/2;
         arg1=r1(dir)-resp(:,dir)+sqrt((r1(1)-resp(:,1)).^2+(r1(2)-resp(:,2)).^2+(r1(3)-resp(:,3)).^2);
         arg2=r2(dir)-resp(:,dir)+sqrt((r2(1)-resp(:,1)).^2+(r2(2)-resp(:,2)).^2+(r2(3)-resp(:,3)).^2);
         arg1=arg2./arg1;
         arg1=log(arg1);
         arg1=arg1.*js{source}(1:end-1,dir);
         h=h+sum(arg1(arg1<Inf));


end