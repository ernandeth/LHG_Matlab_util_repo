pf = fopen('myVSI_5500.rho','rb');
fseek(pf, 64, 'bof');
old_rho = fread(pf,5500 , 'int16','ieee-be');
fclose(pf);

t=[1:5500]*0.008;
plot(t, old_rho,'r')


pf = fopen('myVSI_5500.theta','rb');
fseek(pf, 64, 'bof');
old_theta = fread(pf,5500 , 'int16','ieee-be');
fclose(pf);

hold on
plot(t, old_theta,'g')


pf = fopen('myVSI_5500.grad','rb');
fseek(pf, 64, 'bof');
old_grad = fread(pf,5500 , 'int16','ieee-be');
fclose(pf);

plot(t, old_grad, 'b')


% 
figure
t=[1:11000]'*0.004;
new_rho=load('myVSI_11000.rho.txt');
new_theta=load('myVSI_11000.theta.txt');
new_grad=load('myVSI_11000.grad.txt');

plot(t, new_rho,'r')
hold on
plot(t, new_theta,'g')
plot(t, new_grad, 'b')