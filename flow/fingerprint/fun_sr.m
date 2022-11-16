function x=fun_sr(xdata,t)

    x=xdata(2)-xdata(3).*exp(-t./xdata(1));

end