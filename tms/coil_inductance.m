% some calculations for design of TMS coils and drivers:

% inductance of a bunch of concentric loops

N = 15;
r = 0.03;
d = 0.02;

L = 1e6 * (N*r)^2 /((2*r+2.8*d)*1e5)


% a strainght wire : inductance of the cable

l = [0.01:0.01:1];
L = l.*(log(4*l/0.001)-1)*200e-9;

plot(l,L), grid on


% dependence of B- field on radius

B = zeros(10,20);
for p=1:20
    for r=1:10
        
        B(r,p) = 2*pi*r / sqrt(r^2 + p^2);
    end
end

plot(B)
xlabel('Radius')
ylabel('B-field')

plot(B')
xlabel('penetration depth')
ylabel('B-field')

    