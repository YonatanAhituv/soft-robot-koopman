%%% Van der Pol oscillator
%%%    x'=(landa-y^2)*x-p^2*y
%%%    y'=x
%%% ref: Balanov, Alexander, et al. Synchronization: From Simple to Complex. Springer, 2009.
clear all;close all;
num_train = 20;
num_valid = 10;

data = struct;
data.train = {};
for j = 1:num_train+num_valid
data_out = struct;
landa=rand();  % bifurcation parameter
p =rand();   % determines the frequency of oscilations

y=rand();  % initial condtion
x=rand();  % initial condtion
out(1)=y;  % initial condtion


dt = 0.01;   
t = -0.01;
T_start=2;
T_end=5000;
for i= T_start:T_end
   t = t + dt;
   time(i)=t;
   u(i) = rand();
   x=dt*((landa-y.^2)*(x-p.^2*y + 5*u(i))) + x;
   y = dt*x + y;

   out(i)=y;   % keeping the output in differernt time steps
end
data_out.u = u';
data_out.x = out';
data_out.y = out';
data_out.t = time';
data.train = [data.train data_out];
end
data.val = data.train(num_train:end);
data.train(num_train:end) = [];

save('datafiles/vanderpol.mat','-struct','data');
figure;
hold on

plot (time, out)  % time domain
xlabel('time');ylabel ('y')
figure
plot(out(1:end-1),diff(out))  % phase space
xlabel('y');ylabel ('x')
