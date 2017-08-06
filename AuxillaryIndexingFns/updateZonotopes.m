function [IC_z_generators_copy, alphas, scaledICs, R_zs, IC_z_copy, IC_c_copy] = ...
    updateZonotopes(alpha_gx, IC_z_generators_copy, IC_c_copy, alphas, dim, s, scaledICs, A_d, R_zs)
         % add scalars from the current iteration to the vector of scalars
        for j=1:dim
            alphas(generatorIndex(dim, s)+j) = alpha_gx(j);
        end
        
        % update the generator matrix to construct a scaling of the initial set
        % later, will need to update center in the interval implementation as well
        IC_z_g =[];
        for i=1:dim
            IC_z_g = horzcat(IC_z_g, alpha_gx(i)*IC_z_generators_copy(:,i));
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
        % make a generator matrix for the reach set in the previous
        % step, scale the generators and construct the new initial set
        IC_z_generators_copy = get(IC_z_copy, 'Z');
        IC_z_generators_copy = IC_z_generators_copy(:, 2:length(IC_z_generators_copy));   
        
end