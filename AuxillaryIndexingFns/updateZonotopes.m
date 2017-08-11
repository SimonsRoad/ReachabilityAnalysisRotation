% scaling and reach set updating for implementation with lower and upper
% bound scalars on generators

function [IC_z_generators_copy, alphas, center_shift, scaledICs, R_zs, IC_z_copy, IC_c_copy] = ...
    updateZonotopes(alpha_gx_low, alpha_gx_high, IC_z_generators_copy, IC_c_copy, alphas, center_shift, dim, s, scaledICs, A_d, R_zs)
         % add scalars from the current iteration to the vector of scalars
        for j=1:dim
            alphas(generatorIndex(dim, s)+j) = (alpha_gx_high(j) - alpha_gx_low(j))/2;
        end
        
        % update cumulative center shift
        center_shift = center_shift + (A_d^(-1))^(s-1)*(IC_z_generators_copy*(alpha_gx_high + alpha_gx_low)/2);
        
        % update the center of the scaled initial set
        IC_c_copy = IC_c_copy + IC_z_generators_copy*(alpha_gx_high + alpha_gx_low)/2;
        
        % update the generator matrix to construct scaling of the initial set      
        IC_z_g =[];
        for i=1:dim
            IC_z_g = horzcat(IC_z_g, (alpha_gx_high(i) - alpha_gx_low(i))/2*IC_z_generators_copy(:,i));
        end
        IC_z_generators_copy = IC_z_g;
        
        % construct scalings of initial sets
        IC_z_copy_mat = horzcat(IC_c_copy, IC_z_generators_copy);
        IC_z_copy = zonotope(IC_z_copy_mat);
        scaledICs = horzcat(scaledICs, IC_z_copy);


        % construct a reach set
        reachSet = A_d*IC_z_copy;
        R_zs = horzcat(R_zs, reachSet);
        IC_c_copy = center(reachSet);
        
        % the next initial condition is the reach set from this
        % iteration
        IC_z_copy = reachSet;
        % make a generator matrix for the new reach set
        IC_z_generators_copy = get(IC_z_copy, 'Z');
        IC_z_generators_copy = IC_z_generators_copy(:, 2:end);   
        
end
