% scales the generators only; also, introduced extra 
% generators that come from the uniform distribution over the n-sphere.
% when having multiple time steps, the scaling from the previous time step
% is preserved and new scalings are applied iteratively. 
% the first scaling is contained in the initial set, subsequent scalings
% are forced to lie within the previous reach set

% parameters of interest:
% - timeStep and tFinal
% - IC and CS (centered or skewed)
% - d and extraScale (number of random generators and their length)
% - alpha_gx(i) >= ...; alpha_gx(i) <= ... bounds on generator scalars

% TIME
%----------------------------------------------------------------------
% potentially will need to get rid of the similar fields in the options
tStart=0; %start time
tFinal= 3; %final time
timeStep= 1; %time step size for reachable set computation
% number of linearization points
number_steps = (tFinal-tStart)/timeStep;
%----------------------------------------------------------------------

% INITIAL CONDITIONS
%----------------------------------------------------------------------
IC = interval([-1;-1], [1;1]);
%IC = interval([-1;-1], [2;3]);
%IC = interval([1;-2], [3.5;1]);


IC_z = zonotope(IC);
IC_z_generators = get(IC_z, 'Z');
IC_z_generators = IC_z_generators(:, 2:length(IC_z_generators));
IC_c = center(IC_z);

% add extra generators
d = 10;
dim = d+2;
extraScale = 10/d;
extraG = extraScale*generatePts(d, 2);
IC_z_generators = horzcat(IC_z_generators, extraG);
IC_mat = horzcat(IC_c, IC_z_generators);
IC_z = zonotope(IC_mat);
%----------------------------------------------------------------------

% CONSTRAINTS
%----------------------------------------------------------------------
CS = interval([-1.1;-1.2], [1.2;1.3]);
%CS = interval([-2.1;-3.2], [4;5]);
%CS = interval([-2;-3], [4;3]);

CS_z = zonotope(CS);
%----------------------------------------------------------------------

% MATRICES
%----------------------------------------------------------------------
% 90 degree rotation
A = [0   -1;
     1   0];
% transformation from continuous to discrete time
%A_d = @(t) expm(A*t*timeStep);
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

    if s==1
        cvx_begin
            variable alpha_gx(dim)
            maximize sum(alpha_gx)     

            subject to
            for i=1:dim
                alpha_gx(i) >= 0.1 %100/d;
                %alpha_gx(i) <= 0.2 %10/d;
            end

            % scaled initial set is inside unscaled one
            IC_c_copy + abs(IC_z_generators_copy)*alpha_gx <= supremum(IC);
            IC_c_copy - abs(IC_z_generators_copy)*alpha_gx >= infimum(IC);

            % reach sets are within constraint set
            A_d*IC_c_copy + abs(A_d*IC_z_generators_copy)*generatorWeights(alpha_gx, dim, 1) <= supremum(CS);
            A_d*IC_c_copy - abs(A_d*IC_z_generators_copy)*generatorWeights(alpha_gx, dim, 1) >= infimum(CS);

        cvx_end
        
       [IC_z_generators_copy, alphas, scaledICs, R_zs, IC_z_copy] = ...
           updateZonotopes(alpha_gx, IC_z_generators_copy, IC_c_copy, alphas, dim, s, scaledICs, A_d, R_zs);
    
    else
            cvx_begin
                variable alpha_gx(dim)
                maximize sum(alpha_gx)     

                subject to
                for i=1:dim
                    alpha_gx(i) >= 0.1 %100/d;
                    % need this constraint to make sure that scalings are
                    % contained in the original sets
                    alpha_gx(i) <= 1 
                end

                % scaled initial condition is within the previous reach set
                IC_c_copy + abs(IC_z_generators_copy)*alpha_gx <= supremum(interval(IC_z_copy));
                IC_c_copy - abs(IC_z_generators_copy)*alpha_gx >= infimum(interval(IC_z_copy));

                % reach sets are within constraint set
                A_d*IC_c_copy + abs(A_d*IC_z_generators_copy)*generatorWeights(alpha_gx, dim, 1) <= supremum(CS);
                A_d*IC_c_copy - abs(A_d*IC_z_generators_copy)*generatorWeights(alpha_gx, dim, 1) >= infimum(CS);

            cvx_end
        
       [IC_z_generators_copy, alphas, scaledICs, R_zs, IC_z_copy] = ...
           updateZonotopes(alpha_gx, IC_z_generators_copy, IC_c_copy, alphas, dim, s, scaledICs, A_d, R_zs);
   end

   
end

alphas
size(alphas)

% PLOTS
%----------------------------------------------------------------------
figure;
hold on;

plot(IC, [1,2], 'y','lineWidth',2);
plot(CS, [1,2], 'g','lineWidth',2);
pause(1);

% plot scalings of initial set
for i=1:number_steps
    plot(scaledICs(i), [1,2], 'r','lineWidth',i);
    pause(2);
    plot(R_zs(i), [1,2],'b','lineWidth',i);
    pause(2)
end




