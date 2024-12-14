clear;
close all;
data = struct;
data.train = {};
n_dims = 2;
for j = 1:20
    % Simulation parameters
    dt = 0.01;  % time step
    t_end = 10; % total simulation time
    t = 0:dt:t_end; % time vector
    data_out = struct;
    data_out.t = transpose(t);
    for k = 1:n_dims
        % Initial condition
        x0 = rand; % random initial condition between 0 and 1
        x = zeros(size(t)); % preallocate state vector
        x(1) = x0; % set initial condition
        
        % Input u: sinusoidal function
        u = rand*10*sin(rand*2*pi*t);
        u_vals = zeros(size(t));
        
        % Simulation using Euler method
        for i = 1:length(t)-1
            dx = -(x(i)*sin(x(i))) + u(i); % compute derivative
            u_vals(i) = u(i);
            x(i+1) = x(i) + dx * dt; % update state
        end
        % Update data file
        data_out.x(:, k) = x; 
        data_out.y(:, k) = x; 
        data_out.u(:, k) = u_vals; 
    end
   

    data.train = [data.train, data_out]; % Append
    
    
    % Plot results
    figure;
    subplot(2, 1, 1);
    plot(t, x);
    title('State x over time');
    xlabel('Time (s)');
    ylabel('x');

    subplot(2, 1, 2);
    plot(t, u);
    title('Input u over time');
    xlabel('Time (s)');
    ylabel('u = sin(t)');
end
data.val = data.train(11:20);
data.train(11:20) = [];

save('datafiles/dummyData.mat','-struct','data');