function [u,r] = solveSystem(ndof,K,f,up,vp)
    vf = setdiff((1:ndof)',vp);
    u = zeros(ndof,1);
    u(vp) = up;
    LHS = K(vf,vf);
    RHS = f(vf) - K(vf,vp)*u(vp);
    
    method = 'Direct';      % method = Direct or Iterative

    solver = Solver.create(LHS,RHS,method);
    solver.compute();

    u(vf) = solver.x;
    r = K(vp,:)*u-f(vp);
end