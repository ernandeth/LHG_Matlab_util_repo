function [xdatafit fit] =fun_bat2(xdata,t,y)

    [xdatafit] = fmincon(@batfun, xdata, [],[],[],[],[0 0 0 0],[10 10 10 10],[],optimset('MaxFunEvals',100000,'MaxIter',100000,'TolX', 1e-8 ));

    [SSE fit] = batfun(xdatafit);

   

    

    function  [sse, x]=batfun(xdata)
         x=xdata(2).*(t-xdata(1)).^xdata(3).*exp(-(t-xdata(1))./xdata(4));
         x(x<0)=0;
    
         sse = sum((y - x).^2);
        
    end
end