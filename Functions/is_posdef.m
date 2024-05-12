function is_posdef(H)
% Is H positive definite?

[R,P] = chol(H);

if P > 0
error('H is not positive definite')
else
disp('H is positive definite')
end

end