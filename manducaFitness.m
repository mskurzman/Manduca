% EN01-04 Computational Modeling and Design
% Department of Computer Science
% Tufts University

% This function returns the distance covered by Manduca's head.
% Inputs:
%	Two 10-row matrices; each row is one of the 10 time segments.
%	- 'muscles' is a 10x4 matrix of muscle forces; all entries are 100 or 0
%	- 'legs' is a 10x5 matrix of legs locked. 1 means locked, 0 means not.
%	- 'record': if 0 or not supplied, then no recording.
%	  If 1, then save the end-of-time-segment values to a file.
%	  If 2, then save all values to a file.
% Returns:
%	- the distance covered by the center of mass.
%	- result_details: 10x5 matrix, with one row per timepoint (i.e., at
%	  t=10,...,90,100). The first 5 columns are the positions of the 5 legs.
%	  Then the COM position and velocity.
%
% The system is just 5 legs & 4 springs+muscles; there is no intrinsic
% difference between the head and the tail. However, since we start in the
% position leg #1 x=0 and leg #5 x=2000, then it makes sense to call leg #5
% the head. This means that increasing x values correspond to moving the head
% forwards.
function [COM_distance,resultDetails] = manducaFitness(legs, muscles, record)
    % Error checking our input parameters.
    nTimeSeg = size (legs,1);	% Number of time segments (usually 10).
    %assert (nTimeSeg == 10);	% Comment this out if desired.
    assert (isequal(size(legs), [nTimeSeg,5]));	   % 10 time slots x 5 legs.
    assert (isequal(size(muscles), [nTimeSeg,4])); % 10 time slots x 4 muscles.
    assert (isequal (abs(legs-.5), repmat (.5, nTimeSeg, 5))); % Legs are 1 or 0
    assert (isequal (abs(muscles-50), repmat (50, nTimeSeg, 4))); %Musc 100 or 0
    narginchk(2,3);	% Either 2 or 3 arguments; 'record' is optional.
    if (nargin==2)	% If 'record' is not supplied, default it to false.
	record=0;
    end;

    % If we are told to record, then open the output file to record into.
    if (record>0)
	fileID = open_outfile();
    end
    resultDetails=zeros(nTimeSeg,7);	% Allocate debug-return matrix.

    % Set up the time steps.
    timeInterval =[0 100];		% Start time & end time.
    stepSize = 100/nTimeSeg;		% Time length per mini-simulation.
    tpts = timeInterval(1):stepSize:timeInterval(2); % 1x11 vector of timepoints

    % Initialize the running state vector. We will update it after each sim.
    x0 = [0; 500; 1000; 1500; 2000];	% Positions  of the 5 legs at t=0.
    xprime0 = zeros(5, 1);		% Velocities of the 5 legs at t=0.

    % The distance covered by Manduca will be the difference between the head
    % of Manduca (at x=2000) at time =0 vs. t=100.
    manducaCOM_T0 = mean (x0);		% Location of the COM at t=0

    % The main simulation loop.
    endL = length(tpts)-1;		% Idx of next-to-last timepoint (ie., 10)
    for i = 1:endL			% For each time interval.
	timeInterval_ode = [tpts(i), tpts(i+1)]; % Start & end of this interval.
	lockedLegs_row = legs(i, :);	% Leg-locked row for this interval.
	muscles_row = muscles(i, :);	% Muscle row for this interval

	% Call fivePointManduca() to do the work.
	[t x] = fivePointManduca(x0, xprime0, muscles_row, ...
				 lockedLegs_row, timeInterval_ode);

	% The 'end' row gives the values at the end of the simulation interval.
	% We use them to make 5x1 column vectors for position & velocity of the
	% 5 legs, to use as initial conditions for the next time interval.
	x0 = x(end, 1:5)';		% position of the 5 legs at sim end
	xprime0 = x(end, 6:10)';	% velocity of ""

	% Save our end-of-sim results for this segment.
	if (record ~= 0)
	    if (record==2) start=1; else start=size(t); end
	    for r=start:size(t)
		fprintf (fileID, '%f,  %f,%f,%f,%f,%f, %d,%d,%d,%d,%d, %d,%d,%d,%d\n', ...
			 t(r), x(r,1),x(r,2),x(r,3),x(r,4),x(r,5), ...
			 lockedLegs_row(1),lockedLegs_row(2),lockedLegs_row(3),...
			 lockedLegs_row(4),lockedLegs_row(5),...
			 muscles_row(1),muscles_row(2),muscles_row(3),muscles_row(4));
	    end
	end
	resultDetails(i,1:5) = x(end,1:5);
	resultDetails(i,6) = mean(x(end,1:5));	% COM position
	resultDetails(i,7) = mean(x(end,6:10));	% COM velocity
    end	% of for each of the 10 timepoints.

    % 'x0' is left over from the final simulation interval. We grab the final
    % position of the COM to compute the distance traveled.
    COM_distance = mean(x0) - manducaCOM_T0;

    if (record)
	status = fclose(fileID);
	assert (status==0, 'Cannot close recording file manduca_output.txt');
    end
end

% fivePointManduca() simulates one time period. It just calls the library
% function ode45(), passing it our actual ODE function 'dfile()'.
% Inputs:
%	- initialPositions and initialVelocities are each a 5x1 vector of
%	  initial conditions (one per leg).
%	- muscles is the 1x4 vector of muscle forces.
%	- legLocked2 is the 1x5 vector of legs locked.
%	- timeInterval is the 1x2 vector of [start_time,end_time].
% Return [t,x] straight from ode45:
%	- 't' is a column vector of intermediate timepoints that ode45 chose.
%	- 'x' is a matrix with the same # of rows as 't'. Each such row has
%	  one column for each variable, and contains the variable values at
%	  the corresponding timepoint.
function [t, x] = fivePointManduca(initialPositions, initialVelocities, ...
				  muscles, legLocked2, timeInterval)
    global M12 M23 M34 M45 legLocked

    % If a leg is locked, its position is constant. Its velocity is thus 0 the
    % entire interval. dfile() makes sure that its acceleration is 0; we now
    % ensure that its velocity starts at 0.
    if (legLocked2(1))
	initialVelocities(1)=0;
    end
    if (legLocked2(1))
	initialVelocities(2)=0;
    end
    if (legLocked2(1))
	initialVelocities(3)=0;
    end
    if (legLocked2(1))
	initialVelocities(4)=0;
    end
    if (legLocked2(1))
	initialVelocities(5)=0;
    end

    % Set globals to communicate the leg-frozen & muscle-force choices to
    % dfile(). Since the integration variables are just leg positions &
    % velocity, dfile() wouldn't otherwise know the parameters.
    % Nesting dfile() inside fivePointManduca() would be a cleaner solution.
    legLocked=legLocked2;	% 1x5 row vector of legs locked.
    M12 = muscles(1);
    M23 = muscles(2);
    M34 = muscles(3);
    M45 = muscles(4);
    [t,x] = ode45(@dfile, timeInterval,[initialPositions; initialVelocities]);
end

% The ODE-implementation function dfile(), which gets passed to ode45().
% Inputs:
%	- 'x' is a 10x1 column vector of the variables at time 't'.
% Outputs:
%	- 'xprime', a 10x1 column vector of the first derivatives.
% The first 5 variables are leg position; then next 5 are leg velocity.

% In this intuitive description, we'll call the ten variables x1-x5 and v1-v5.
% The first 5 eqns: just say that the derivative of a position variable
% is its respective velocity variable.
%	x1' = locked1? 0:v1.
%	x2' = locked2? 0:v2.
%	x3' = locked3? 0:v3.
%	x4' = locked4? 0:v4.
%	x5' = locked5? 0:v5.
% The next 5 eqns are f=ma. Note the viscous-damping terms work in the same
% direction as the Hookian terms. So if you stretch a spring very quickly, it
% fights back more. Real-life metal springs don't do this much, but real-life
% polymers (and most biomaterials) do. Since our spring is actually soft tissue,
% it makes sense for it to have visco-elasticity.
%	v1' = (k(x2-x1-L0) + c(v2-v1) + M12)/m
%	v2' = (k(x3-x2-L0) - k(x2-x1-L0) + c(v3-v2) + c(v1-v2)+ M23-M12)/m
%	v3' = (k(x4-x3-L0) - k(x3-x2-L0) + c(v4-v3) + c(v2-v3)+ M34-M23)/m
%	v4' = (k(x5-x4-L0) - k(x4-x3-L0) + c(v5-v4) + c(v3-v4)+ M45-M34)/m
%	v5' =              (-k(x5-x4-L0)            + c(v4-v5)      - M45)/m
% To model the ground reaction forces on a locked leg:
% So if legs #4 and #5 are locked, then x4' and
% x5' are correctly forced to 0, which correctly forces x4 and x5 to be
% constant. 
% we would want to set v4'=0 and v5'=0 rather than doing f=ma on them.

function xprime = dfile(t, x)
    c = 2;	% Viscosity constant.
    m = 1;	% Mass.
    k = 1;	% Elastic spring constant.
    L0 = 500;	% Resting spring length.
    global M12 M23 M34 M45 legLocked

    xprime = zeros(10,1);
    xprime(1) = x(6); 
    xprime(2) = x(7); 
    xprime(3) = x(8); 
    xprime(4) = x(9); 
    xprime(5) = x(10); 

    xprime(6)= (1/m) * (- k * x(1) + k * x(2) - c * x(6) + c * x(7) - k* L0  + M12);

    xprime(7)= (1/m)* ( -k *x(2) + k * x(3) - c * x(7) + c * x(8) - k * L0 + M23) +...
	 (1/m) * (k * x(1) - k * x(2) + c * x(6) - c * x(7) +  k * L0  - M12);

    xprime(8)= (1/m)* ( -k *x(3) + k * x(4) - c * x(8) + c * x(9) - k * L0 + M34) +...
	 (1/m) * (k * x(2) - k * x(3) + c * x(7) - c * x(8) +  k * L0  - M23);

    xprime(9)= (1/m)* ( -k *x(4) + k * x(5) - c * x(9) + c * x(10) - k * L0 + M45) +...
	 (1/m) * (k * x(3) - k * x(4) + c * x(8) - c * x(9) +  k * L0  - M34);

    xprime(10)= (1/m)* (-k * x(5) + k * x(4) - c * x(10) + c * x(9) + k * L0 - M45);

  
    if (legLocked(1))
	xprime(1)=0;
	xprime(6)=0;
    end;
    if (legLocked(2))
	xprime(2)=0;
	xprime(7)=0;
    end;
    if (legLocked(3))
	xprime(3)=0;
	xprime(8)=0;
    end;
    if (legLocked(4))
	xprime(4)=0;
	xprime(9)=0;
    end;
    if (legLocked(5))
	xprime(5)=0;
	xprime(10)=0;
    end;
end

function fileID = open_outfile()
    full_name = 'manduca_output.txt';
    fprintf ('Recording to file %s.\n', full_name);
    [fileID,msg] = fopen(full_name,'w');
%    assert (fileID>3, 'Problem *%s* opening %s', msg, full_name);
    fprintf (fileID, 'Time x_leg1 x_leg2 x_leg3 x_leg4 x_leg5 lock1..5 musc1..4\n');
end
