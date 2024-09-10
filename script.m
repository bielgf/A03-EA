%% A03 Beam structure: Wingbox

%% DATA -------------------------------------------------------------------

data.E = 210*10^(9);                          % Pa
data.G = 80*10^(9);                           % Pa

data.chord     = 2;                           % m
data.wingspan  = 16;                          % m
data.chiP      = 0.3*data.chord;              % m
data.beamWidth = 0.3*data.chord;              % m

data.h1 = 0.25*data.chord;                    % m
data.h2 = 0.15*data.chord;                    % m
data.t1 = 22*10^(-3);                         % m
data.t2 = 15*10^(-3);                         % m
data.t3 = 3.5*10^(-3);                        % m

data.aeroCenter   = 0.25*data.chord;          % m
data.gravCenter   = 0.48*data.chord;          % m
data.zEngine      = 0.3*data.chord;           % m
data.xEngine      = 0.25*data.wingspan;       % m
data.lambda       = 140;                      % Kg/m

data.engineMass = 2100;                       % Kg
data.Cl         = 0.1;
data.rhoInf     = 1.225;                      % kg/m^3
data.vInf       = 750*5/18;                   % m/s

data.g = 9.81;                                % m/s^2


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
data.xBendMoment = -1;
data.yBendMoment = 0;
data.zBendMoment = 1;
data.xShearForce = 0;
data.yShearForce = 1;

% Material properties matrix
data.materialProp = [ % Thickness
                    data.t2
                    data.t3
                    data.t1
                    ];

%% B) Beam Analysis -------------------------------------------------------

% B.0) Beam's Geometric Discretization
data.numElements   = 16;                      % Number of elements
data.numNodesElem  = 2;                       % Number of nodes in an element
data.numDOFsNode   = 3;                       % Number of degrees of freedom per node


% B.1)
% Nodal connectivities matrix
data.nodalConnec = [1:data.numElements ; 2:data.numElements+1]';

% Material connectivities matrix
data.materialConnec = ones(data.numElements,1);

% Degrees of freedom connectivities matrix
data.dofsConnec = connectDOF(data.numElements,data.numDOFsNode,data.nodalConnec);

% Fixed nodes matrix
data.fixedNodes = [ % node, direction, value
                    1   1   0
                    1   2   0
                    1   3   0
                   ];

%% OOP --------------------------------------------------------------------

s = data; % NO!!

s.geomParams.beamWidth  = data.beamWidth;
s.geomParams.chord      = data.chord;
s.geomParams.wingspan   = data.wingspan;
s.geomParams.chiP       = data.chiP;
s.geomParams.aeroCenter = data.aeroCenter;
s.geomParams.gravCenter = data.gravCenter;
s.geomParams.h1         = data.h1;
s.geomParams.h2         = data.h2;
s.geomParams.N1         = data.N1;
s.geomParams.N2         = data.N2;
s.geomParams.N3         = data.N3;
s.geomParams.open       = data.open;

s.connec.nodalConnec    = data.nodalConnec;
s.connec.materialConnec = data.materialConnec;
s.connec.dofsConnec     = data.dofsConnec;
s.connec.fixedNodes     = data.fixedNodes;

s.aeroParams.Cl     = data.Cl;
s.aeroParams.rhoInf = data.rhoInf;
s.aeroParams.vInf   = data.vInf;
s.aeroParams.g      = data.g;


m = FEMBeamComputer(s);
m.compute();

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

% stycnv = studyOfCnvg(data);
% stycnv.compute();



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

close all;