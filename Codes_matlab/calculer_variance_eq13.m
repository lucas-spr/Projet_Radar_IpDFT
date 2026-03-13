function var_th = calculer_variance_eq13(delta, A, M, sigma_n, H) %Implémentation de la fonction de l'équation 13
    % La limite de x/sin(x) en 0 est 1, donc 2*pi*delta / sin(pi*delta) -> 2
    if abs(delta) < 1e-10
        bloc1 = 2 / A;
    else
        bloc1 = (2 * pi * delta) / (A * sin(pi * delta));
    end
    
    bloc2 = (H - abs(delta)) / factorial(2*H - 1);
    
    bloc3 = 1;
    for h = 1:(H-1)
        bloc3 = bloc3 * (h^2 - delta^2);
    end
    
    n_comb = 4*H - 4;
    k_comb = 2*H - 2;
    C_term = nchoosek(n_comb, k_comb); 
    bloc4 = sqrt(C_term / (2 * M));
    
    bloc5 = sqrt( (2*(4*H - 3)*(delta^2 - abs(delta)) + 2*H^2 - 1) / (2*H - 1) );
    
    % Calcul de l'écart-type
    sigma_delta = bloc1 * bloc2 * bloc3 * bloc4 * bloc5 * sigma_n;
    
    % Calcul de la variance
    var_th = sigma_delta^2;
end