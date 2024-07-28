%% A03 Beam structure: Wingbox

% Change example number 1

% --------------------
clear;clc
close all
% --------------------

%% PREPROCESS

data.E = 210*10^(9);            % Pa
data.G = 80*10^(9);             % Pa

data.c = 2;                     % m
data.b = 16;                    % m
data.xi_p = 0.3*data.c;         % m
data.d = 0.3*data.c;            % m

data.h1 = 0.25*data.c;          % m
data.h2 = 0.15*data.c;          % m
data.t1 = 22*10^(-3);           % m
data.t2 = 15*10^(-3);           % m
data.t3 = 3.5*10^(-3);          % m

data.za = 0.25*data.c;          % m
data.zm = 0.48*data.c;          % m
data.ze = 0.3*data.c;           % m
data.be = 0.25*data.b;          % m
data.lambda = 140;              % Kg/m

data.Me = 2100;                 % Kg
data.Cl = 0.1;
data.rhoinf = 1.225;            % kg/m^3
data.vinf = 750*5/18;           % m/s

data.g = 9.81;                  % m/s^2


%% A) Cross Section Analysis

% A.0) Section's Geometric Discretization
data.N1 = 200;                  % Number of elements of the first section
data.N2 = 400;                  % Number of elements of the second section
data.N3 = 200;                  % Number of elements of the third section
data.open = 1;

if data.open == 0
    section = 'closed';
else
    section = 'open';
end

[x_prim,Tm,Tn] = GeometricDiscret(data);


% A.1) A.2) Open section

% Input values for analysis
M_x_prim = -1;
M_y_prim = 0;
M_z_prim = 1;
S_x_prim = 0;
S_y_prim = 1;

% Material properties matrix
mD1 = [ % Thickness
    data.t2
    data.t3
    data.t1
    ];


% SOLVER ------------------------------------------------------------------

% Function to get section properties
[x_0_prim,y_0_prim,x_s_prim,y_s_prim,A_tot,I_xx_prim,I_yy_prim,I_xy_prim,J,A_in] = sectionProperties(data,x_prim,Tn,mD1,Tm,data.open);

% Function to get the normal stress distribution
[sigma,s_norm] = normalStressDistribution(x_prim,Tn,x_0_prim,y_0_prim,I_xx_prim,I_yy_prim,I_xy_prim,M_x_prim,M_y_prim);

% Function to get the tangential stress distribution due to shear
[tau_s,s_shear] = tangentialStressDistributionShear(x_prim,Tn,mD1,Tm,x_0_prim,y_0_prim,I_xx_prim,I_yy_prim,I_xy_prim,S_x_prim,S_y_prim,x_s_prim,y_s_prim,A_in,data.open);

% Function to get the  tangential stress distribution due to torsion
[tau_t,s_tor] = tangentialStressDistributionTorsion(x_prim,Tn,mD1,Tm,M_z_prim,J,data.open,A_in);


% figure(1)
% plot(s_norm(1,:),sigma(1,:))
% xlabel('Node position, s [m]')
% ylabel('Normal stress, σ [Pa]')
% title(sprintf('Normal stress distribution on the %s section (Mx unitary)',section))
% grid on 
% 
% figure(2)
% plot(s_shear(1,:),tau_s(1,:))
% xlabel('Node position, s [m]')
% ylabel('Tangetial stress, τ [Pa]')
% title(sprintf('Tangential stress distribution due to shear on the %s section (Sy unitary)',section))
% grid on 
% 
% figure(3)
% plot(s_tor(1,:),tau_t(1,:))
% xlabel('Node position, s [m]')
% ylabel('Tangetial stress, τ [Pa]')
% title(sprintf('Tangential stress distribution due to torsion on the %s section (Mz unitary)',section))
% grid on

% plot2DBars(x_prim,Tn,sigma/10^(3),'kPa','Normal Stress Distribution',section)
% plot2DBars(x_prim,Tn,tau_s/10^(3),'kPa','Tangential Stress Distribution',section)
% plot2DBars(x_prim,Tn,tau_t/10^(3),'kPa','Torsional Stress Distribution',section)


%% B) Beam Analysis

% B.0) Beam's Geometric Discretization
data.nel = 16;                      % Number of elements
data.nne = 2;                       % Number of nodes in an element
data.ni = 3;                        % Number of degrees of freedom per node


% Element's force and moment computation
[xnod,fe,me] = GetForceMomentElement(data,x_s_prim);


% B.1)
% Nodal connectivities matrix
TnD2 = [1:data.nel ; 2:data.nel+1]';

% Material connectivities matrix
TmD2 = ones(data.nel,1);

% Degrees of freedom connectivities matrix
TdD2 = connectDOF(data,TnD2);

% Material properties matrix
mD2 = [ % Young's Modulus, Shear Modulus, Bending Inertia, Torsional Inertia 
       data.E   data.G      I_xx_prim       J 
       ];

% Fixed nodes matrix
pD2 = [ % node, direction, value
        1   1   0
        1   2   0
        1   3   0
       ];

FD2 = [
        find(xnod == data.be)       1       -data.Me*data.g
        find(xnod == data.be)       3       -data.Me*data.g*(((data.d + data.xi_p) - x_s_prim) - data.ze)
        ];

data.nnod = size(xnod,2);
data.ndof = data.nnod*data.ni;      % Total number of degrees of freedom


% SOLVER ------------------------------------------------------------------

[Kel] = stiffnessFunction(data,xnod',TnD2,mD2,TmD2);
[fel] = forceFunction(data,xnod',TnD2,fe,me);
[K,F] = assemblyFunction(data,TdD2,Kel,fel);
[up,vp] = applyBC(data,pD2);
[F] = pointLoads(data,F,FD2);
[u,r] = solveSystem(data,K,F,up,vp);
[xel,Sel,Mbel,Mtel] = internalforcesFunction(data,xnod',TnD2,TdD2,Kel,u);


% PLOTTING

% figure(6)
% plot(xnod,u(1:3:end-2))
% title(sprintf('Vertical Deflection for %s section',section))
% grid on
% xlabel('x (m)')
% ylabel('Vertical Deflection (m)')
% 
% figure(7)
% plot(xnod,u(2:3:end-1))
% title(sprintf('Section Bending Rotation for %s section',section))
% grid on
% xlabel('x (m)')
% ylabel('Section Bending Rotation (rad)')
% 
% figure(8)
% plot(xnod,u(3:3:end))
% title(sprintf('Section Torsion Rotation for %s section',section))
% grid on
% xlabel('x (m)')
% ylabel('Section Torsion Rotation (rad)')
%  
% figure(9)
% plot(xel,Sel)
% title(sprintf('Shear Force for %s section',section))
% grid on
% xlabel('x (m)')
% ylabel('Shear Force (N)')
%  
% figure(10)
% plot(xel,Mbel)
% title(sprintf('Bending Moment for %s section',section))
% grid on
% xlabel('x (m)')
% ylabel('Bending Moment (N*m)')
% 
% figure(11)
% plot(xel,Mtel)
% title(sprintf('Torsion Moment for %s section',section))
% grid on
% xlabel('x (m)')
% ylabel('Torsion Moment (N*m)')


% B.2) Study of convergence

el_conv = [4 8 16 32 64 128 256 512];
StudyOfConvergence(el_conv,data,mD2,x_s_prim,pD2,section);


%% C) Von Mises Criterion

[max_vms,mcp_sect,mcp_beam] = VonMises(data,mD1,xnod,Sel,Mbel,Mtel);



%% Unit testing

% Test 1
% load("resultsGood.mat",'KGood')
% error1 = norm(K(:) - KGood(:));
% if error1 < 1e-10
%     disp('Test 1 passed')
% else
%     disp('Test 1 failed')
% end

test1 = stiffnessMatrixTest(K);
run(test1)

% Test 2
% load("resultsGood.mat",'FGood')
% error2 = norm(F(:) - FGood(:));
% if error2 < 1e-10
%     disp('Test 2 passed')
% else
%     disp('Test 2 failed')
% end

test2 = forceTest(F);
run(test2)

% Test 3
% load("resultsGood.mat",'uGood')
% error3 = norm(u(:) - uGood(:));
% if error3 < 1e-10
%     disp('Test 3 passed')
% else
%     disp('Test 3 failed')
% end

test3 = displacementsTest(u);
run(test3)

% Test 4
% load("resultsGood.mat",'rGood')
% error4 = norm(r(:) - rGood(:))/norm(rGood(:));
% if error4 < 1e-10
%     disp('Test 4 passed')
% else
%     disp('Test 4 failed')
% end

test4 = reactionsTest(r);
run(test4)


% transformar MAIN en class fins abans tests -> pasar tests -> commit
