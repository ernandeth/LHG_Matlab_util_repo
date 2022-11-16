function out = nan2zero(in)

if any(isnan(in))
    out = zeros(1,length(in));
else
    out = in;
end

end
