function value = pulse(t, start,N_cycles, amplitude, w)

duration = N_cycles;
if w~=0
    duration = N_cycles/w;
end

value=zeros(size(t));
stimtime=(t>=start & t<=start+duration);
if any(stimtime)
    value(stimtime)=amplitude;
end



if w~=0
    wave = cos(2*pi*w*t);
    value = value.*wave;

end