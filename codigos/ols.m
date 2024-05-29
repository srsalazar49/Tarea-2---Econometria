% Definimos una funcion MCO generica que permite calcular los coeficientes 
% de manera matricial junto a los residuos del modelo

function [beta, e] = ols(Y,X)

% Por formula, beta = (X'X)^(-1) (X'Y)
beta = (X' * X)\(X' * Y);

% Por formula, e = Y - X * beta_gorro
e = Y - (X * beta);

end

