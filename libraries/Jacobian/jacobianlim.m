function  [J,h,func_evals,f_0,E,S,F,X]=jacobianlim(func,x,varargin)
%   JACOBIAN CALCULATION AS A LIMITING PROCESS
%   
%   evaluates the jacobian of a function func at point x as a converging 
%   limiting process. It respects bounds, can accept a baseline evaluation 
%   (e.g. f_0=func(x)), and will provide it for further evaluations to 
%   speed up computations. Further, the step size, no of evaluations, an 
%   estimate of the error size and no of steps taken are also returned. The
%   "inspiration" can be found in "NUMERICAL METHODS Using MATLAB" by John 
%   H. Mathews and Kurtis D. Fink. They just show it for a simple 
%   derivative, and I expanded it to a Jacobian while looking at
%   efficiency. This calculation can also be used to find simple
%   derivatives, btw.
%
%   syntax: 
%
%   [J,h,func_evals,f_0,E,S]=jacobianlim(func,x);
%   [J,h,func_evals,f_0,E,S]=jacobianlim(func,x,opt);
%   [J,h,func_evals,f_0,E,S]=jacobianlim(func,x,name,value);
%   [J,h,func_evals,f_0,E,S]=jacobianlim(func,x,opt,name,value);
%
%   func must accept one argument x which has n elements and it provides
%   one output f that has m elements. x is the point at which the jacobian
%   matrix is evaluated, and opt can change the default variables.
%   The available name-value pairs and fields of the options are explained 
%   below together with the default values.
%   Note that in the last syntax the name-value pair overrides any argument
%   in the option structure.
%
%   J           jacobian matrix of size m x n
%   h           step size, same size as x
%   func_evals  no of function evaluations
%   f_0         evaluation value at func(x)
%   E           matrix of estimates for errors
%   S           matrix of no of steps taken
%   
%
%   examples:
%
%   func=@(X)[sin(X(1));cos(X(2))];
%   x=[0 0];
%   J=jacobianlim(func,x);
%   J2=jacobianlim(func,x'); %yields identical results
%   opt.FinDiffRelStep=eps^(1/4);
%   J3=jacobianlim(func,x,opt);
%   J4=jacobianlim(func,x,'FinDiffRelStep',eps^(1/10));
%   J5=jacobianlim(func,x,opt,'FinDiffRelStep',eps^(1/10));
%
%   The last calculation is done with FinDiffRelStep=eps^(1/10) as the 
%   name-value pair overrides the structure.
%
%   author: Alexander.Dentler at gmail.com, October 9th 2015, all errors
%   are yours to deal with. 
%
%   version 3
     
%%  SET OPTIONS
f_0=[];                     %   value at evaluation point (can be user-supplied
lb=-inf(size(x));           %   lower boundary of domain of X
ub=inf(size(x));            %   upper boundary of domain of X
FinDiffRelStep=eps^(1/3);   %   relative differentiation step size
TypicalX=1;                 %   sets lower limit of differentiation step 
%                               size, similar to the optimization routines 
%                               in matlab 
DecreaseStepSize=10;        %   factor that reduces the step size in limiting process
MaxStepNo=3;                %   max no of steps for the limiting process
toler=1e-6;                 %   absolute tolerances
%%  DYNAMIC READ IN
if nargin>2    
    list=who;
    opt=[];arg_list=[];
    pos=1;flag_thereismore=1;
    %%  check if argument gives us a structure
    if flag_thereismore
        if isstruct(varargin{pos})
            opt=varargin{pos};
            pos=pos+1;
            flag_thereismore=nargin>(pos+1);
        end
    end
    %%  check for name-variable pairs
    if flag_thereismore
        if ((nargin-pos-1)/2)==fix((nargin-pos-1)/2)
            arg_list=varargin(pos:end);
        else
            error('No of arguments is off.')
        end
    end  
    %%  add option structure variables if they are part of the list
    if ~isempty(opt)
        for i1=1:numel(list)
            if isfield(opt, char(list{i1}))
                eval(horzcat(matlab.lang.makeValidName(char(list(i1))),'=opt.',char(list{i1}),';'));
            end
        end
    end
    %%  add name-value pair arguments if they are part of the list
    if ~isempty(arg_list)
        for i1=1:numel(arg_list)/2
            if ismember(arg_list{(i1-1)*2+1},list) && isnumeric(arg_list{(i1-1)*2+2})
                eval(horzcat(arg_list{(i1-1)*2+1},'=',num2str(arg_list{(i1-1)*2+2}),';'));
            end
        end
    end 
end
%%  PREP
n=numel(x);
func_evals=0;
if isempty(f_0)
    f_0=func(x);
    func_evals=func_evals+1;
end
m=numel(f_0);
J=NaN(m,n);
E=NaN(m,n);
S=ones(m,n);
F=NaN(m,2*MaxStepNo+1,n);
X=NaN(n,2*MaxStepNo+1,n);
%   size of increment before boundaries are considered
h=FinDiffRelStep.*max(abs(x),TypicalX);  
%%  ROUTINE FOR DIFFERENT PARAMETERS
for k=1:n
    %%  find final step size for given boundaries    
    if isfinite(lb(k)) || isfinite(ub(k))
        h(k)=min([h(k);abs(x(k)-lb(k))/2;abs(x(k)-ub(k))/2]);
        if h(k)==0
            error('evaluation on boundary not possible.')
        end
    end
    %%  prep for this loop
    D=NaN(m,MaxStepNo);
    D_r=NaN(m,MaxStepNo);
    D_l=NaN(m,MaxStepNo);
    Q=inf(m,MaxStepNo);
    R=zeros(m,MaxStepNo);
    Q_sided=inf(m,MaxStepNo);
    R_sided=zeros(m,MaxStepNo);
    %%  first approximation of derivative
    [D(:,1),D_l(:,1),D_r(:,1),x_up,f_up,x_down,f_down]=deriv(func,f_0,x,h(k),k);
    func_evals=func_evals+2;
    F(:,[1 MaxStepNo+1 end],k)=[f_down(:) f_0(:) f_up(:)];
    X(:,[1 MaxStepNo+1 end],k)=[x_down(:) x(:) x_up(:)];
    Q_sided(:,1)=abs(D_l(:,1)-D_r(:,1));
    R_sided(:,1)=2*Q_sided(:,1).*(abs(D_l(:,1))+abs(D_r(:,1))+eps);
    %   see if first step is already sufficient
    i1=1;
    ind= find(R_sided(:,i1)<toler & ~isnan(D(:,i1)));
    J(ind,k)=D(ind,i1);
    E(ind,k)=Q_sided(ind,i1);
    S(ind,k)=i1;
    %%  HIGHER LIMITS
    if any(isnan(J(:,k)))
        h(k)=h(k)/DecreaseStepSize;
        for i1=2:MaxStepNo
            %%  evaluation
            [D(:,i1),D_l(:,i1),D_r(:,i1),x_up,f_up,x_down,f_down]=deriv(func,f_0,x,h(k),k);
            func_evals=func_evals+2;
            F(:,[i1 end-i1+1],k)=[f_down(:) f_up(:)];
            X(:,[i1 end-i1+1],k)=[x_down(:) x_up(:)];
            Q(:,i1)=abs(D(:,i1)-D(:,i1-1));
            R(:,i1)=2*Q(:,i1).*(abs(D(:,i1))+abs(D(:,i1-1))+eps);
            Q_sided(:,i1)=abs(D_l(:,i1)-D_r(:,i1));
            R_sided(:,i1)=2*Q_sided(:,i1).*(abs(D_l(:,i1))+abs(D_r(:,i1))+eps);
            %%  current error tolerance is small enough
            ind= find( ~(Q(:,i1-1)<Q(:,i1)) & ( R(:,i1)<toler | R_sided(:,i1)<toler)  & ~isnan(D(:,i1)));%& isnan(J(:,order)) 
            J(ind,k)=D(ind,i1);
            E(ind,k)=Q(ind,i1);
            S(ind,k)=i1;
            %%  error is growing and we take the last iterations value
            ind= find( ~(Q(:,i1-1)>=Q(:,i1)) & isnan(J(:,k)) & ~isnan(D(:,i1-1)));
            J(ind,k)=D(ind,i1-1);
            E(ind,k)=Q(ind,i1-1);
            S(ind,k)=i1-1;
            %%  if error was growing, make sure future value growth reflects this
            Q(:,i1-1)=min(Q(:,i1-1),Q(:,i1));     
            if all(~isnan(J(:,k)))
                break
            else
                h(k)=h(k)/DecreaseStepSize;
            end
        end
    end
    %   collecting results
    ind=find(isnan(J(:,k)));
    J(ind,k)=D(ind,end);
    E(ind,k)=Q(ind,end);
    S(ind,k)=MaxStepNo;    
end
end

function [out,out_l,out_r,x_up,f_up,x_down,f_down]=deriv(func,f_0,x,h,order)
%   right point
x_up=x;
x_up(order)=x(order)+h;
%   right evaluation
f_up=func(x_up);f_up=f_up(:);
%   left point
x_down=x;
x_down(order)=x(order)-h;
%   left evaluation
f_down=func(x_down);f_down=f_down(:);
%   central derivative
out=diff([f_down(:) f_up]')'/(2*h);
%   left derivative
out_l=diff([f_down f_0(:)]')'/h;
%   right derivative
out_r=diff([f_0(:) f_up]')'/h;
end