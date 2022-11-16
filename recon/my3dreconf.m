function [obj im1]=my3dreconf(nechos,timepoint,xres,centric)
%nechos=32;
offset=nechos*(timepoint-1);
discards=0;

tmp=hanning(nechos);
tmp=ones(nechos,1);

zvec=[nechos/2:-1:-nechos/2+1]
zshift=0;

if(nargin<3)
    xres=64;
end
if(nargin<4)
    centric=0;
end

if(centric)
tmp1=nechos/2:-1:1;
tmp2=nechos/2+1:nechos;
encvals=zeros(1,nechos);
encvals(1:2:end)=tmp2;
encvals(2:2:end)=tmp1;

%encvals=[9,8,10,7,11,6,12,5,13,4,14,3,15,2,16,1];
  for n=1+discards:nechos+discards;  
    eval(sprintf('im1_mag(:,:,nechos+1-encvals(n))=readim(''sl1.%03d'');',n+offset)); 
    eval(sprintf('im1_phs(:,:,nechos+1-encvals(n))=readim(''sl1.phs.%03d'');',n+offset)); 
    im1_phs(:,:,nechos+1-encvals(n))=im1_phs(:,:,nechos+1-encvals(n))/1000;
    im1(:,:,nechos+1-encvals(n))=im1_mag(:,:,nechos+1-encvals(n)).*exp(-i.*im1_phs(:,:,nechos+1-encvals(n))).*tmp(nechos+1-encvals(n)-discards);
  end;
else
  for n=(1+discards):(nechos+discards);  
  
   if(exist('sl1.1.001','file')==0)
      eval(sprintf('im1_mag(:,:,n)=readim(''sl1.%03d'',[xres xres]);',n+offset)); 
      eval(sprintf('im1_phs(:,:,n)=readim(''sl1.phs.%03d'',[xres xres]);',n+offset));
      %2*pi*.0625.*zvec(n).*zshift;
      im1_phs(:,:,n)=im1_phs(:,:,n)/1000;    %+2*pi*.0625.*zvec(n).*zshift.*ones(xres);
      im1(:,:,n)=im1_mag(:,:,n).*exp(-i.*im1_phs(:,:,n)).*tmp(n-discards);
    else
      for nn=1:8
        eval(sprintf('im1_mag(:,:,n,nn)=readim(''sl1.%d.%03d'',[xres xres]);',nn,n+offset));
        eval(sprintf('im1_phs(:,:,n,nn)=readim(''sl1.%d.phs.%03d'',[xres xres]);',nn,n+offset));
      end
      im1_phs(:,:,n,nn)=im1_phs(:,:,n,nn)/1000;
      im1coils(:,:,n,nn)=im1_mag(:,:,n,nn).*exp(-i.*im1_phs(:,:,n,nn)).*tmp(n-discards);

    end
  end
   if(exist('sl1.1.001','file')~=0)
      im1=mean(im1coils,4);
  end

end

% 
% for x=1:xres
%     for y=1:xres
%         im1_phs(x,y,:)=interp1(1:nechos,squeeze(im1_phs(x,y,:)),1+zshift:nechos+zshift,'linear');
%         im1_mag(x,y,:)=interp1(1:nechos,squeeze(im1_mag(x,y,:)),1+zshift:nechos+zshift,'cubic');
%     end
% end
% 
% im1_phs(find(isnan(im1_phs)))=0;
% im1_mag(find(isnan(im1_mag)))=0;
% im1=im1_mag.*exp(-i.*im1_phs);




obj=abs(fftshift(ifft(fftshift(im1,3),[],3),3));

%figure,tile3d(obj)


docorrect=0;
if(docorrect)
    TE=25;
    corrfact=exp([TE:TE:nechos*TE]/350);
    for n=1:nechos;  im1corr(:,:,n)=im1(:,:,n)./corrfact(n); end;
    corrobj=abs(fftshift(ifft(fftshift(im1corr,3),[],3),3));
    figure,tile3d(corrobj)
end


% my_max_pez = (int)(2*PI*1000000*my_zres/(2*my_fovz/10*26754));  /* area in G*us/cm     my_fovz in mm*/
