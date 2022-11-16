function padborder(hAxs,pct)

% EXAMPLE
% padborder(hAxs,pct);

for i = 1:length(hAxs)
    
    hAx = hAxs(i);
    oa = axis(hAx);

    xrng = (100+pct)*diff(oa(1:2))/100;
    yrng = (100+pct)*diff(oa(3:4))/100;

    xmid = mean(oa(1:2));
    ymid = mean(oa(3:4));

    set(hAx,'ylim',ymid+[-yrng yrng]./2);
    set(hAx,'xlim',xmid+[-xrng xrng]./2);
    
end