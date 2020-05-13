function [J,h,func_evals,E,S,D]=jacobianext(func,x,varargin)

%   JACOBIAN CALCULATION AS A ROMBERG EXTRAPOLATION
%   
%   evaluates the jacobian of a function func at point x using romberg  
%   extrapolation. It respects bounds. Further, the step size, no of 
%   evaluations, an estimate of the error size and no of steps taken are 
%   also returned. The"inspiration" can be found in "NUMERICAL METHODS
%   Using MATLAB" by John H. Mathews and Kurtis D. Fink. They just show it 
%   for a simple derivative, and I expanded it to a Jacobian while looking 
%   at efficiency. This calculation can also be used to find simple
%   derivatives, btw.
%
%   syntax: 
%
%   [J,h,func_evals,E,S,D]=jacobianext(func,x);
%   [J,h,func_evals,E,S,D]=jacobianext(func,x,opt);
%   [J,h,func_evals,E,S,D]=jacobianext(func,x,name,value);
%   [J,h,func_evals,E,S,D]=jacobianext(func,x,opt,name,value);
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
%   E           matrix of estimates for errors
%   S           matrix of no of steps taken
%   D           extrapolation matrix
%
%   examples:
%
%   func=@(X)[sin(X(1));cos(X(2))]';
%   x=[0 0]';
%   J=jacobianext(func,x);
%   J2=jacobianext(func,x'); %yields identical results
%   opt.FinDiffRelStep=eps^(1/4);
%   J3=jacobianext(func,x,opt);
%   J4=jacobianext(func,x,'FinDiffRelStep',eps^(1/10));
%   J5=jacobianext(func,x,opt,'FinDiffRelStep',eps^(1/10));
%
%   The last calculation is done with FinDiffRelStep=eps^(1/10) as the 
%   name-value pair overrides the structure. 
%   You can pass on arguments. Say you use this routine for an optimization 
%   routine, and you find no better point given J. Then you can save D from 
%   your last Jacobian evaluation and recycle:
%
%   opt=[];
%   opt.MaxStepNo=1;
%   opt.FinDiffRelStep=eps^(1/10);
%   [J6,h6,func_evals6,E6,S6,D6]=jacobianext(func,x,opt);
%   opt.D=D6;
%   [J7,h7,func_evals7,E7,S7,D7]=jacobianext(func,x,opt);
%
%   this refines the jacobian while avoiding possibly expensive function
%   calls.
%
%   author: Alexander.Dentler at gmail.com, October 9th 2015, all errors
%   are yours to deal with. 
%
%   version 3
   
%%  SET OPTIONS
lb=-inf(size(x));           %   lower boundary of domain of X
ub=inf(size(x));            %   upper boundary of domain of X
FinDiffRelStep=eps^(1/3);   %   relative differentiation step size
TypicalX=1;                 %   sets lower limit of differentiation step 
%                               size, similar to the optimization routines 
%                               in matlab 
MaxStepNo=3;                %   max no of steps for the limiting process
toler=1e-6;                 %   absolute tolerances
%%  DYNAMIC READ IN
D=[];   %   old Romberg matrix from earlier evaluation that can be reused. Inofficial part.
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
m=1;                        %   no of output arguments
if isempty(D)
    D=NaN(MaxStepNo,MaxStepNo,n,m);
else
    siz=size(D);siz([1 2])=siz([1 2])+MaxStepNo;
    m=siz(4);
    D_new=NaN(siz);
    D_new(1:size(D,1),1:size(D,2),:,:)=D;
    D=D_new;
end
siz=size(D);
I=siz(1);
Q=inf(siz(1),n,m);
R=inf(siz(1),n,m);
J=NaN(m,n);
E=NaN(m,n);
S=ones(m,n);
%   size of increment before boundaries are considered
h=FinDiffRelStep.*max(abs(x),TypicalX);  
func_evals=0;
%%  LOOP
for k=1:n   
    %%  find final step size for given boundaries    
    if isfinite(lb(k)) || isfinite(ub(k))
        h(k)=min([h(k);abs(x(k)-lb(k))/2;abs(x(k)-ub(k))/2]);
        if h(k)==0
            error('evaluation on boundary not possible.')
        end
    end
    %%  prep for this loop
    current_J=sum(~isnan(D(:,1,k,1)));
    for j=2:current_J
        Q(j,k,:)=abs(D(j,j,k,:)-D(j-1,j-1,k,:));
        R(j,k,:)=2*squeeze(Q(j,k,:))./squeeze(abs(D(j,j,k,:))+abs(D(j-1,j-1,k,:))+eps);
    end
    %%  loop to increase depth
    flag_out=0;
    for j=(current_J+1):I
        if j==2
            if ~any(any(~(R(1:j-1,k,:)<toler),1),3)
                flag_out=1;
            end
        elseif j>2            
            if ~any(any(~(R(1:j-1,k,:)<toler),1) & all(Q(1:j-2,k,:)>=Q(2:j-1,k,:),1),3)
                flag_out=1;
            end
        end
        if flag_out==0
            %   only go on if one jacobian derivative element has not
            %   reduced its relative error below toler, and that element is
            %   still no suffering from numerical problems
            x_l=x;
            x_l(k)=x(k)-2^(-j)*h(k);
            x_r=x;
            x_r(k)=x(k)+2^(-j)*h(k);
            temp=(func(x_r)-func(x_l))/(2^(-j+1)*h(k));
            if m<numel(temp) && k==1 && j==1 % first evaluation, and we have more than one output argument
                m=numel(temp);
                Q(:,:,2:m)=inf;
                R(:,:,2:m)=inf;
            end
            D(j,1,k,1:m)=temp;%(func(x_r)-func(x_l))/(2^(-j+1)*h(k));
            func_evals=func_evals+2;
            for l=1:(j-1)
                D(j,l+1,k,:)=D(j,l,k,:)+(D(j,l,k,:)-D(j-1,l,k,:))/((4^l)-1);
            end
            if j>1
                Q(j,k,1:m)=abs(D(j,j,k,:)-D(j-1,j-1,k,:));
                R(j,k,1:m)=2*squeeze(Q(j,k,:))./squeeze(abs(D(j,j,k,:))+abs(D(j-1,j-1,k,:))+eps);
            end
        else
            break
        end
    end
    %%  choose best fit for jacobian
    temp=Q(:,k,:);
    temp(~cumprod(Q(1:end-1,k,:)>=Q(2:end,k,:),1))=NaN;
    [err,pos]=min(temp,[],1);
    for i1=1:m
        J(i1,k)= D(pos(i1),pos(i1),k,i1);
        E(i1,k)= err(i1);
        S(i1,k)= pos(i1);
    end
end
end

