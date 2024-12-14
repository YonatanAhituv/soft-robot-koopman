clear;
close all;
%%
%% Load model and val data
[ model_name , model_path ] = uigetfile('systems/fromData/*.mat' , 'Choose a model file...' );
koopman_model = load([model_path,model_name]).sysid_class;
[ datafile_name , datafile_path ] = uigetfile( 'datafiles/*.mat' , 'Choose a data file for validation and initial conditions...' );
valdata = load([datafile_path,datafile_name]).val{1,2};

% % Simple Dummy Data
% valdata = struct;
% valdata.x = [[0 1];
%              [2 3]
%              [4 5]];
% valdata.y = valdata.x;
% valdata.u = valdata.x+1;
% valdata.t = [0,1,2];
% 
% valdata = struct;
% valdata.x = [[-3];
%     [-2]
%     [-1]];
% valdata.y = [[0];
%     [1]
%     [2]];;
% valdata.u = [[3];
%     [4]
%     [5]];
% valdata.t = [0,1,2];

zeta = koopman_model.get_zeta(valdata);

lifted_zeta = koopman_model.lift.full(zeta.zeta(1,:)'); % Linear
% lifted_zeta = koopman_model.lift.full( [ zeta.zeta(1,:) , zeta.uzeta(1,:) ]' ); % Nonlinear

%% Add noise
for i = 1:size(valdata.u,1)
    valdata.u(i,1) = valdata.u(i,1) + rand*0.1;
    % valdata.u(i,2) = valdata.u(i,2) + rand*0.1;
end
for i = 1:size(valdata.x,1)
    valdata.x(i,1) = valdata.x(i, 1) + rand*0.1;
    % valdata.x(i,2) = valdata.x(i, 2) + rand*0.1;
end
for i = 1:size(valdata.y,1)
    valdata.y(i,1) = valdata.y(i,1) + rand*0.1;
    % valdata.y(i,2) = valdata.y(i,2) + rand*0.1;
end



% system = koopman_model.get_model(koopman_model.koopData); % Linear Model
% results = koopman_model.val_model(system,valdata);
system = koopman_model.get_NLmodel(koopman_model.koopData); % Nonlinear Model
results = koopman_model.val_NLmodel(system,valdata);


%%
%-------------------------------------------------------------------------%
%---------------------------- Plot Results -------------------------------%
%-------------------------------------------------------------------------%

% Define color(s)
cb_blue = [55,126,184] ./ 255;  % blue
cb_red = [255,0,0] ./ 255;  % red
figure;
hold on;
plot( results.sim.t , results.sim.y , 'color' , cb_red , 'LineWidth' , 2,'DisplayName','Koopman Approximation');
plot( results.real.t , results.real.y , 'color' , cb_blue , 'LineWidth' , 2,'DisplayName','Real');
legend;
% plot(  results.sim.t , results.sim.u, 'color' , cb_blue , 'LineWidth' , 2);
% plot(  results.sim.t , results.sim.z, 'color' , cb_blue , 'LineWidth' , 2);
xlabel('Time (s)' , 'Interpreter' , 'Latex' , 'FontSize', 20);
ylabel('$y$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 20);
hold off;
grid on; box on;





