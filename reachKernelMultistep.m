% scalars for lower and upper bounds are not symmetric.
% scales the generators only; also, introduced extra 
% generators that come from the uniform distribution over the n-sphere.
% when having multiple time steps, the scaling from the previous time step
% is preserved and new scalings are applied iteratively. 
% the first scaling is contained in the initial set, subsequent scalings
% are forced to lie within the previous reach set

% parameters of interest:
% - timeStep and tFinal
% - IC and CS (initial and constraint sets)
% - d and extraScale (number of random generators and their length)

% TIME
%----------------------------------------------------------------------
% potentially will need to get rid of the similar fields in the options
tStart=0; %start time
tFinal=5; %final time
timeStep=1; %time step size for reachable set computation
% number of linearization points
number_steps = (tFinal-tStart)/timeStep;
%----------------------------------------------------------------------

% INITIAL CONDITIONS
%----------------------------------------------------------------------
%IC = interval([-1;-1], [1;1]);
%IC = interval([-1;-1], [2;3]);
IC = interval([1;-2], [3.5;1]);


IC_z = zonotope(IC);
IC_z_generators = get(IC_z, 'Z');
IC_z_generators = IC_z_generators(:, 2:length(IC_z_generators));
IC_c = center(IC_z);

% add extra generators
d = 20;
dim = d+2;
extraScale = 10/d;
extraG = extraScale*generatePts(d, 2);
IC_z_generators = horzcat(IC_z_generators, extraG);
IC_mat = horzcat(IC_c, IC_z_generators);
IC_z = zonotope(IC_mat);
%----------------------------------------------------------------------

% CONSTRAINTS
%----------------------------------------------------------------------
%CS = interval([-1.1;-1.2], [1.2;1.3]);
%CS = interval([-2.1;-3.2], [4;5]);
CS = interval([-2;-3], [4;3]);

CS_z = zonotope(CS);
%----------------------------------------------------------------------

% MATRICES
%----------------------------------------------------------------------
% 90 degree rotation
A = [0   -1;
     1   0];
% transformation from continuous to discrete time
A_d = expm(A*timeStep);
%----------------------------------------------------------------------

% OPTIMIZATION
%----------------------------------------------------------------------
% keep track of initial set scalings
scaledICs = [];
%scaledICs = horzcat(scaledICs, IC_z);

% keep track of the reach sets
R_zs = [];

% initialize vector of coefficients
alphas = zeros(dim*number_steps, 1);

% initialize intermediate values for center and generator matrix
IC_c_copy = IC_c;
IC_z_generators_copy = IC_z_generators;


            
for s=1:number_steps   
        cvx_begin
            variable alpha_gx_low(dim)
            variable alpha_gx_high(dim)
            
            maximize sum(alpha_gx_high) - sum(alpha_gx_low)     

            subject to
            for i=1:dim
                alpha_gx_low(i) >= -1;
                alpha_gx_low(i) <= 1;
                alpha_gx_high(i) >= -1;
                alpha_gx_high(i) <= 1;
                % make sure that additional generators don't degenerate
                alpha_gx_low(i)+0.1 <= alpha_gx_high(i);
            end
            
            % scaled initial set is inside unscaled one, need this because
            % of extra generators added
            if s==1
                IC_c_copy + IC_z_generators_copy*(alpha_gx_high + alpha_gx_low)/2 ...
                 + abs(IC_z_generators_copy)*(alpha_gx_high - alpha_gx_low)/2  <= supremum(IC);
                 IC_c_copy + IC_z_generators_copy*(alpha_gx_high + alpha_gx_low)/2 ...
                 - abs(IC_z_generators_copy)*(alpha_gx_high - alpha_gx_low)/2 >= infimum(IC);
            end
            
            % scaled initial condition is within the constraint set
           IC_c_copy + IC_z_generators_copy*(alpha_gx_high + alpha_gx_low)/2 ...
           + abs(IC_z_generators_copy)*(alpha_gx_high - alpha_gx_low)/2  <= supremum(CS);
           IC_c_copy + IC_z_generators_copy*(alpha_gx_high + alpha_gx_low)/2 ...
           - abs(IC_z_generators_copy)*(alpha_gx_high - alpha_gx_low)/2 >= infimum(CS);
            
            % reach sets are within constraint set
           A_d*(IC_c_copy + IC_z_generators_copy*(alpha_gx_high + alpha_gx_low)/2) ...
               + abs(A_d*IC_z_generators_copy)*(alpha_gx_high-alpha_gx_low)/2 <= supremum(CS);
           A_d*(IC_c_copy + IC_z_generators_copy*(alpha_gx_high + alpha_gx_low)/2) ...
                - abs(A_d*IC_z_generators_copy)*(alpha_gx_high-alpha_gx_low )/2 >= infimum(CS); 

        cvx_end
        
       [IC_z_generators_copy, alphas, scaledICs, R_zs, IC_z_copy, IC_c_copy] = ...
           updateZonotopes(alpha_gx_low, alpha_gx_high, IC_z_generators_copy, IC_c_copy, alphas, dim, s, scaledICs, A_d, R_zs);   
end

% ACCUMULATE SCALARS - APPLY TO INITIAL SET
%----------------------------------------------------------------------
IC_z_mat =[];
for j=1:dim
    tempGenerator = IC_z_generators(:,j);
    for i=1:number_steps
         tempGenerator = alphas(generatorIndex(dim, i)+j)*tempGenerator;
    end
    IC_z_mat = horzcat(IC_z_mat, tempGenerator);
end
IC_z_mat = horzcat(IC_c, IC_z_mat);
IC_z_g = zonotope(IC_z_mat);
%----------------------------------------------------------------------


% PLOTS
%----------------------------------------------------------------------
figure;
hold on;

plot(IC, [1,2], 'y','lineWidth',2);
plot(CS, [1,2], 'g','lineWidth',2);
plot(IC_z_g, [1,2],'y','lineWidth',2);
pause(1);

% plot scalings of initial set
for i=1:number_steps
    plot(scaledICs(i), [1,2], 'r','lineWidth',2);
    pause(1);
    plot(R_zs(i), [1,2],'b','lineWidth',2);
    pause(1)
end
