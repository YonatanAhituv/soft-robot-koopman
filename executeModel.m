clear; 
close all;
%%
ksysid = load('systems/fromData/n-2_m-2_del-1_2024-10-23_13-05.mat').sysid_class;
kmpc = kmpc( ksysid ,...
        'horizon' , 25 ,...
        'input_bounds' , [ 0 , 10 ],... 
        'input_slopeConst' , [0.5e-2],... 
        'input_smoothConst' , [1e-1],... 
        'state_bounds' , [] ,...
        'cost_running' , 0.1 ,...
        'cost_terminal' , 100 ,...
        'cost_input' , 0 ); 
    
%-------------------------------------------------------------------------%
%-------------------------- Run Simulation -------------------------------%
%-------------------------------------------------------------------------%

% load in reference trajectory from file
[ reffile_name , reffile_path ] = uigetfile( 'ref-trajectories/*.mat' , 'Choose reference trajectory file...' );
temp = load( [reffile_path , reffile_name] );
ref = temp.ref;

% run simulation
y0 = [ 1 ,1 ];      % initial laser dot position
u0 = [ 0 , 0];  % initial input


% run_simulation: Runs a simulation of system under mpc controller
function results = run_simulation( obj , ref , y0 , u0)
    %run_trial: Runs a simulation of system under mpc controller.
    %   Tries to follow the trajectory in ref and impose the
    %   shape constraints in shape_bounds.
    %   Assume ref and shape_bounds have same sampling frequency as
    %   sytem, and that they are already scaled to be consistent 
    %   with the lifted model.
    %   ref - struct containing reference trajectory with fields:
    %       t - vector of timesteps
    %       y - each row is a desired point at the corresponding timestep
    %   x0 - [1,n] initial condtion
    %   u0 - [1,m] initial input
    
    % shorthand
    nd = obj.params.nd;
    Np = obj.horizon;
    
    % set value of initial conditions to zero if none provided
    if nargin < 3
        y0 = zeros( nd+1 , obj.params.n );
        u0 = zeros( nd+1 , obj.params.m );
    elseif nargin < 4
        y0 = kron( ones( nd+1 , 1 ) , y0 );
        u0 = zeros( nd+1 , obj.params.m );
    else
        y0 = kron( ones( nd+1 , 1 ) , y0 );
        u0 = kron( ones( nd+1 , 1 ) , u0 );
    end
    
    % resample and scale the reference trajectory
    ref_Ts = obj.resample_ref( ref );
    ref_sc = obj.scaledown.y( ref_Ts );
    
    % set initial condition
    initial.y = y0; initial.u = u0;
    [ initial , zeta0 ] = obj.get_zeta( initial );    % LINE NOT NEEDED
    
    % initialize results struct
    results = struct;
    results.T = [ 0 ];
    results.U = [ u0( end , : ) ];
    results.Y = [ y0( end , : ) ];
    results.K = [ 0 ];
    results.R = [ ref.y(1,:) ];
    results.X = [ y0( end , : ) ];
    results.Z = [ obj.lift.full( zeta0' )' ]; % lifted states
    
    k = 1;
    while k < size( ref_sc , 1 )
        
        % current time
        t = k * obj.params.Ts;
        
        % get current state and input with delays
        if k == 1
            current.y = obj.scaledown.y( y0 );   
            % current.u = u0;  
        elseif k < nd + 1
            y = [ y0( k : end-1 , : ) ; results.Y ];
            % u = [ u0( k : end-1 , : ) ; results.U ];
            current.y = obj.scaledown.y( y );
            % current.u = obj.scaledown.u( u ); 
        else
            y = results.Y( end - nd : end , : );
            % u = results.U( end - nd : end , : );
            current.y = obj.scaledown.y( y ); 
            % current.u = obj.scaledown.u( u ); 
        end
        current.u = u0;  
        % isolate the reference trajectory over the horizon
        if k + Np <= size( ref_sc , 1 )
            refhor = ref_sc( k : k + Np , :);
        else
            refhor = ref_sc( k : end , : );     % repeat last entry
        end 
        
        % get optimal input over horizon
        [ U , z ] = obj.get_mpcInput( current , refhor );
        
        % if a solution was not found, break out of while loop
        % if any( isnan(U) )
        %     break;
        % end
        
        % isolate input for this step (may need to make it U(1,:)
        % u_kp1_sc = U( 2 , : );
        
        % scaleup the input for the results
        u_kp1 = obj.scaleup.u(u0)';
        
        % simulate the system over one time-step
        z_k = z;
        u_k_sc = obj.scaledown.u( results.U(end,:) );  % need to use previously calculated input NEED TO MAKE THIS CLEANER!!!
        z_kp1 = obj.model.A * z_k + obj.model.B * u_k_sc';
        x_kp1 = obj.model.C * z_kp1;
        y_kp1_sc = x_kp1;  % output juse scaled version of state since model was learned from observations
        y_kp1 = obj.scaleup.y( y_kp1_sc' )';  % scale output back up 
        
        % record updated results
        results.T = [ results.T ; t ];
        results.U = [ results.U ; u_kp1' ];
        results.Y = [ results.Y ; y_kp1' ];
        results.K = [ results.K ; k ];
        results.R = [ results.R ; obj.scaleup.y( ref_sc( k , : ) ) ];   % note that this is not scaled down
        results.X = [ results.X ; x_kp1' ];
        results.Z = [ results.Z ; z'  ]; % current lifted state
        
        k = k + 1;  % increment step counter
    end
end


sim = run_simulation(kmpc, ref , y0 , u0 );


%%
%-------------------------------------------------------------------------%
%---------------------------- Plot Results -------------------------------%
%-------------------------------------------------------------------------%

% Define color(s)
cb_blue = [55,126,184] ./ 255;  % blue

figure;
hold on;
plot( sim.R(:,1) , sim.R(:,2) , 'color' , [0.25 0.25 0.25 1] , 'LineWidth' , 2);
plot( sim.Y(:,1) , sim.Y(:,2) , 'color' , cb_blue , 'LineWidth' , 2);
xlabel('$x_1$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 20);
ylabel('$x_2$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 20);
hold off;
grid on; box on;





