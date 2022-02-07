function r = RMSE(x1,x2)

MSE = mean((x1-x2).^2);
r = sqrt(MSE);

end

