function alpha = generatorWeights(alpha_gx, dim, i)

        index = generatorIndex(dim, i);        
        alpha = alpha_gx(index+1:index+dim);
             
end