function fh = ExponentialDecay
    fh.fun = @ExpDecay;
    fh.jac = @ExpDecayJac;
end

% Dynamic forward models, Jacobian and gradient functions
function prediction = ExpDecay(x, t)
 
    metAmpl     = x(1,:);
    metDecay    = x(2,:);

    prediction = metAmpl .* exp(-t./metDecay);
    
end


function [jac] = ExpDecayJac(x, t)
    metAmpl     = x(1,:);
    metDecay    = x(2,:);   
   
    dYdmetAmpl = exp(-t./metDecay);
    dYdmetDecay = -t .* metAmpl .* exp(-t./metDecay);

    jac = cat(3,dYdmetAmpl,dYdmetDecay);
    jac = squeeze(jac);
end

