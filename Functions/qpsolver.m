function [x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit)

%   Ax <= b           (inequality constraints)
%   Ax  = b           (equality   constraints)
%   lb <= x <= ub     (bound      constraints)
%
% Solving an inequality constrained QP of the form:
%   min  phi(x) = 0.5*x'*H*x + g*x     where     x is a member of the set R^n
%   such that lb  <= x  <= ub         (restriction 1)
%             bl  <= Ax <= bu,        (restriction 2)
% 
% Input parameters
%   H     : The Hessian matrix of the objective function, size n x n, must
%           be positive definite
%   g     : The linear term of the objective function, can be a vector or a
%           matrix, usually of size n x 1
%   A     : Matrix (limitations), usually of size k x n (A = [A  ; -A])
%           Should be m x n, m is number of inequalities and n = n in x0
%   b     : Vector (limitations), usually of size k x 1 (b = [bu ; bl])
%           Should be m x 1, where m is number of inequality constraints
%           (related to the A matrix)
%   Aeq   : Linear equality constraints, usually of size 1 x k
%           Should be Me x n matrix, Me is number of equalities, n = n in
%           x0
%   beq   : Vector, usually of size 1 x 1 I think
%           Linear equality constraints, should be Me x 1
%   bl    : Lower limits of general constraints
%   bu    : Upper limits of general constraints
%       lb    : Lower bounds, usually of size k x 1 (part of b)
%       lu    : Upper bounds, usually of size k x 1 (part of b)
%   x0    : Starting point (initial value), size is n x 1
% 
% Output parameters
%   x     : The optimal solution
%   info  : Performance information, vector with 3 elements:
%           info(1)   = Final values of the objective function
%           info(2)   = No. of iteration steps
%           info(3)   = Feasible solution found

AA   = [A ; -A];    bb   = [bu ; -bl];
Aeq = [];           beq = []; 

[x,fval] = quadprog(H,g,AA,bb,Aeq,beq,l,u,xinit);
info = fval;

end