%%  PREP
clear;clc;close

%%  SLOW CONVERGENCE
%   This is an example of slow convergence for a fixed point problem, and
%   how to use the algorithms that come with this folder to solve them, in
%   particular when there is a nested contraction.
%
%   Lets use a simple fixed point problem that converges relatively slowly:
%   x=f1(x)=(x+a)^(1/1.01)
%   An iterative scheme would follow this logic:
%   x(i1)=(x(i1-1)+a)^(1/1.01)
%   Set a=1; then
abstol=1e-8;
a1=1;
x1(1)=1e2;
i1=1;
while 1
    i1=i1+1;
    x1(i1,:)=(x1(i1-1)+a1)^(1/1.1);
    if abs(x1(i1)-x1(i1-1))<abstol
        break
    end
end
figure
plot(x1)

%%  NESTED CONTRACTION
%   But lets say there are two contractions, 1 and 2, and lets assume that  
%   the "a1" parameter is the inverse of the fixed point of other contraction, 
%   or 1/a1=x2=f2(x2). Further, assume that 1/a2=x1=f1(x1). We could try to 
%   solve them jointly, but we use a nested approach to learn a) how nested 
%   contractions look like, and b) how to pass on arguments.
%
%   Our initial guess is always
x2Guess=1e9;
x1Guess=x2Guess;

%%  Lets start with an explicit nested iterative scheme:
x2=x2Guess;
x1=x1Guess;
iter_1=0;
iter_2=0;
tic
while 1
    %   interior contraction
    while 1
        iter_2=iter_2+1;
        %   interior iterative scheme
        x2_new=(x2+1/x1)^(1/1.01);
        fval_2=abs(x2_new-x2);
        x2=x2_new;
        if fval_2<1e-8
            break
        end
    end
    iter_1=iter_1+1;
    %   exterior iterative scheme
    x1_new=(x1+1/x2)^(1/1.01);
    fval_1=abs(x1_new-x1);
    x1=x1_new;
    if fval_1<1e-8
        break
    end
end
time_used=toc;
results=table(x1,x2,fval_1,fval_2,iter_1,iter_2,time_used);

%%  Lets see how the line fitting works for this problem. That is we use 
%   line fitting to determine x1 and an iterative scheme for x2. In a first
%   instance, lets not pass the argument through the LM algorithm but fix 
%   it through the function handle.
obj=@(x1Guess)nested_contraction_obj(x1Guess,x2Guess);
%   some minor tweaks for LM
opt.Display='notify';
opt.Jacobian='limit';
%   lets get rid of all stopping criteria
opt.FooTol=NaN; 
opt.RelFooTol=NaN;
opt.FooTol=NaN;
opt.RelTolX=NaN;
%   note that we need to square here
opt.AccTol=abstol^2;
%   Let us count how many times the objective function is called, and how
%   many times the interior contraction is called up
global counter_global
counter_global=0;
tic
[x1,~,fval_1,~,extra_arguments_isempty] =LevenbergMarquardt(obj,x1Guess,[],[],opt);
[~,x2,fval_2]=nested_contraction_obj(x1,x2Guess);
time_used=toc;
%   Note how extra_arguments_isempty is empty.
%   In order to get the second parameter we needed to run the file solo 
%   with an optimal x1.
iter_1=counter_global(1);
iter_2=counter_global(2);
results=[results; table(x1,x2,fval_1,fval_2,iter_1,iter_2,time_used)];

%%  Now lets pass the argument through the LM algorithm. All we need to do 
%   is create an extra argument, and give an initial guess into the option 
%   structure.
obj=@(x1Guess,x2Guess)nested_contraction_obj(x1Guess,x2Guess);
opt.extra_arguments={x2Guess};
counter_global=0;
tic
[x1,~,fval_1,~,extra_arguments] =LevenbergMarquardt(obj,x1Guess,[],[],opt);
%   Lets retrieve the second parameter:
x2=extra_arguments{1};
%   We have the second argument, but we dont have the convergence criterium
%   of the interior contraction. (We could pass it similar to the starting
%   value, but thats left to the user). So:
[~,~,fval_2]=nested_contraction_obj(x1,x2);
time_used=toc;
iter_1=counter_global(1);
iter_2=counter_global(2);
results=[results; table(x1,x2,fval_1,fval_2,iter_1,iter_2,time_used)];

%%  Lets compare the three approaches:
results.Properties.RowNames={'Nested iterative schemes','Line fitting','L-f with arg. passing'};
disp(results)

%   Now, x1 and x2 are in a fixed point with all three algorithms. That is 
%   implied by cross validation and by the value of the innovation fval_1 
%   and fval_2. The nested iterative scheme is the fastest in terms of time
%   but it also requires many calls of the function. This can be very 
%   expensive in more complicated problems. The line fitting is much more
%   efficient, and it gets even more efficient when we pass a good starting
%   point for x2.
