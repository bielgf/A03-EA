%% A03 Beam structure: Wingbox

% ToDo:
% - [BeamSolver] Crear funcions privades de forces, assembly... similar a StiffnessFunctionClass
% - computeForceMomentElem simplificar propietats
% - Fer [F11,F21,F13,F23] = ... private function
% - Unir beamParams i beamProp a BeamSolver i dependències
% - Revisar noms de les classes

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

% A.1) A.2) Open section

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

s.geomParams.E            = data.E;
s.geomParams.G            = data.G;
s.geomParams.beamWidth    = data.beamWidth;
s.geomParams.chord        = data.chord;
s.geomParams.wingspan     = data.wingspan;
s.geomParams.chiP         = data.chiP;
s.geomParams.aeroCenter   = data.aeroCenter;
s.geomParams.gravCenter   = data.gravCenter;
s.geomParams.h1           = data.h1;
s.geomParams.h2           = data.h2;
s.geomParams.N1           = data.N1;
s.geomParams.N2           = data.N2;
s.geomParams.N3           = data.N3;
s.geomParams.materialProp = data.materialProp;
s.geomParams.open         = data.open;

s.beamParams.numElements  = data.numElements;
s.beamParams.numNodesElem = data.numNodesElem;
s.beamParams.numDOFsNode  = data.numDOFsNode;
s.beamParams.lambda       = data.lambda;
s.beamParams.engineMass   = data.engineMass;
s.beamParams.xEngine      = data.xEngine;
s.beamParams.zEngine      = data.zEngine;

s.aeroParams.Cl     = data.Cl;
s.aeroParams.rhoInf = data.rhoInf;
s.aeroParams.vInf   = data.vInf;
s.aeroParams.g      = data.g;

s.connec.nodalConnec    = data.nodalConnec;
s.connec.materialConnec = data.materialConnec;
s.connec.dofsConnec     = data.dofsConnec;
s.connec.fixedNodes     = data.fixedNodes;


m = FEMBeamComputer(s);
m.compute();


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