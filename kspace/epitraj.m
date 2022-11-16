x1 = linspace(-8, 8, 16).';
x2 = linspace(8, -8, 16).';
y = linspace(-8, 8, 16);

y = repmat(y, [16 1]);
y = y(:);

x = x1;
for i=2:16
    if ~mod(i, 2)
        x = cat(1, x, x2);
    else
        x = cat(1, x, x1);
    end
end

plot(x,y, 'LineWidth', 2);
axis([-10 10 -10 10]);
