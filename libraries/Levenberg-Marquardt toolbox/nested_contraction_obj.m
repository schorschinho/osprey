function [fval_1,x2_new,fval_2]=nested_contraction_obj(x1,x2)
%%  the function f(x)=(x+a)^(1/1.01) is a contraction
%   example a) to use arguments that update themself within the routine,
%   and b) of a nested contraction approach.
global counter_global
count=0;
while 1
    count=count+1;
    x2_new=(x2+1/x1)^(1/1.01);
    fval_2=abs(x2_new-x2);
    x2=x2_new;
    if fval_2<1e-8
        break
    end
end
%   count
counter_global=counter_global+[1 count];
%   via line fitting
fval_1=x1-(x1+1/x2)^(1/1.01);
end