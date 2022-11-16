function map=fmap(nsl,xres,mapdel);

if(nargin<3)
    mapdel=2.5e-3;
end
if(nargin<2)
    xres=64;
end

   for n=1:nsl
       eval(sprintf('tp0mag=readim(''sl%d.000'',[xres xres]);',n));
       eval(sprintf('tp0phs=readim(''sl%d.phs.000'',[xres xres])./1000;',n));
       eval(sprintf('tp1mag=readim(''sl%d.001'',[xres xres]);',n));
       eval(sprintf('tp1phs=readim(''sl%d.phs.001'',[xres xres])./1000;',n));
       tp0=tp0mag.*exp(i*tp0phs);
       tp1=tp1mag.*exp(i*tp1phs);
  
       mymask=tp0mag>0.1*max(tp0mag(:));
       %map(:,:,n)=unwrap2d(angle(tp1.*conj(tp0)),mymask,xres/2,xres/2)./mapdel./(2*pi);
       map(:,:,n)=(angle(tp1.*conj(tp0)).*mymask)./mapdel./(2*pi);
 end
   
       
   
