function [beta_random, e_random] = fgls(Y,X,omega)


% Por formula, beta = (X'X)^(-1) (X'Y)
beta_random = inv(X' *inv(omega)* X)*(X'*inv(omega) * Y);


% Por formula, e = Y - X * beta_gorro
e_random = Y - (X * beta_random);

end
