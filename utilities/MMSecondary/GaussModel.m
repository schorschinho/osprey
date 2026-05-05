function F = GaussModel(x,freq)

F = x(1)*exp(x(2)*(freq-x(3)).*(freq-x(3)))+x(4)*(freq-x(3))+x(5);
%NB Integral is x(1)*sqrt(pi/x(2))
end
