function [cphs,cpeak]=res_phs(small)
% Usage ... [phs,peak]=function res_phs(s)
% 
% Caution s contains 2 columns: one is index, two is phase data

% for lp0 = 1:8,
sp = -sign(diff(small(:,2))); % differentiate and take sign
sp2 = sp;
flg = 1;
% loop through data and mark pos to neg zero crossings
for lp = 1:length(sp)
  if sp2(lp) == -1
    if flg == 0
      sp2(lp) = 1;
    end
    flg = 0;
  else
    flg = 1;
  end
end
% get indices of all zero crossings
id = find(sp2-1);
sp3 = sp2;
% set phase between crossings for each point
for lp = 1:length(id)-1
  sp3(id(lp):id(lp+1)) =  ([id(lp):id(lp+1)] - id(lp))./(id(lp+1)-id(lp)); 
end

cphs(1,1:length(small)-1) = sp3';
cphs(1,length(small)) = 1;
cpeak(1,1:length(small)-1) = (1-sp2')./2;
cpeak(1,length(small)) = 0;
%end

