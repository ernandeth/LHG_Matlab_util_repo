function result = shuffle(input)
% function result = shuffle(input)
%
% re-order the data in a vector in a random order.
%
len = length(input);
idx = randperm(len);
result = input(idx);
return