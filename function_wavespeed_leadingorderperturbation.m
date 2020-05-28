function c= function_wavespeed_leadingorderperturbation(kappa,phi)

if kappa > 1
    fun = @(c) c - sqrt(2*( (kappa - c*phi) - (kappa - c*phi)^2*log( (1/(kappa - c*phi))) -(kappa - c*phi)^2));
    
    
    c0 = [0.001 kappa-0.001]; % initial interval - may require updating depending on parameters
    c = fzero(fun,c0);
    
elseif kappa==1
    
    c=0;
    
elseif kappa < 1
    fun = @(c) c + sqrt(2*( (kappa - c*phi) - (kappa - c*phi)^2*log( (1/(kappa - c*phi))) -(kappa - c*phi)^2));
    
    c0 = [(-1+kappa-0.001) -0.001]; % initial interval - may require updating depending on parameters
    c = fzero(fun,c0);
        
end


end
