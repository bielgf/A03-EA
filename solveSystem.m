function [u,r] = solveSystem(data,K,f,up,vp)
    vf = setdiff((1:data.ndof)',vp);
    u = zeros(data.ndof,1);
    u(vp) = up;
    LHS = K(vf,vf);
    RHS = f(vf) - K(vf,vp)*u(vp);
    
    method = 'Iterative';      % method = Direct or Iterative

    solver = Solver.create(LHS,RHS,method);
    solver.compute();

    u(vf) = solver.x;                               % free DOFs
    r = K(vp,:)*u-f(vp);                            % reaction loads
end