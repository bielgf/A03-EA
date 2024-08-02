%% A03 Beam structure: Wingbox

%% DATA -------------------------------------------------------------------

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


%% A) Cross Section Analysis ----------------------------------------------

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

% A.1) A.2) Open section

% Input values for analysis
data.M_x_prim = -1;
data.M_y_prim = 0;
data.M_z_prim = 1;
data.S_x_prim = 0;
data.S_y_prim = 1;

% Material properties matrix
data.mD1 = [ % Thickness
            data.t2
            data.t3
            data.t1
            ];

%% B) Beam Analysis -------------------------------------------------------

% B.0) Beam's Geometric Discretization
data.nel = 16;                      % Number of elements
data.nne = 2;                       % Number of nodes in an element
data.ni = 3;                        % Number of degrees of freedom per node


% B.1)
% Nodal connectivities matrix
data.TnD2 = [1:data.nel ; 2:data.nel+1]';

% Material connectivities matrix
data.TmD2 = ones(data.nel,1);

% Degrees of freedom connectivities matrix
data.TdD2 = connectDOF(data.nel,data.ni,data.TnD2);

% Fixed nodes matrix
data.pD2 = [ % node, direction, value
            1   1   0
            1   2   0
            1   3   0
           ];

%% OOP --------------------------------------------------------------------

m = FEMBeamComputer(data);
m.computeGeoDiscret();

% SECTION SOLVER ----------------------------------------------------------

m.computeSectionSolver();

% BEAM SOLVER -------------------------------------------------------------

m.computeBeamSolver();

%% PLOTTING ---------------------------------------------------------------

% figure(1)
% plot(m.s_norm(1,:),m.sigma(1,:))
% xlabel('Node position, s [m]')
% ylabel('Normal stress, σ [Pa]')
% title(sprintf('Normal stress distribution on the %s section (Mx unitary)',section))
% grid on 
% 
% figure(2)
% plot(m.s_shear(1,:),m.tau_s(1,:))
% xlabel('Node position, s [m]')
% ylabel('Tangetial stress, τ [Pa]')
% title(sprintf('Tangential stress distribution due to shear on the %s section (Sy unitary)',section))
% grid on 
% 
% figure(3)
% plot(m.s_tor(1,:),m.tau_t(1,:))
% xlabel('Node position, s [m]')
% ylabel('Tangetial stress, τ [Pa]')
% title(sprintf('Tangential stress distribution due to torsion on the %s section (Mz unitary)',section))
% grid on
%
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

% el_conv = [4 8 16 32 64 128 256 512];
% StudyOfConvergence(el_conv,data,mD2,x_s_prim,pD2,section);


%% C) Von Mises Criterion

% [max_vms,mcp_sect,mcp_beam] = VonMises(data,mD1,xnod,Sel,Mbel,Mtel);



%% Unit testing

test1 = stiffnessMatrixTest(m.K);
run(test1)

test2 = forceTest(m.F);
run(test2)

test3 = displacementsTest(m.u);
run(test3)

test4 = reactionsTest(m.r);
run(test4)