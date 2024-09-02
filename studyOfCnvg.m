classdef studyOfCnvg < handle

    properties (Access = private)
        data
    end

    methods (Access = public)

        function obj = studyOfCnvg(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeStudy();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.data = cParams;
        end

        function computeStudy(obj)
            
            el_conv = [4 8 16 32 64 128 256 512];
            u_conv = zeros(3,size(el_conv,2));

            for nel = 1:size(el_conv,2)
                obj.data.numElements = el_conv(1,nel);
                obj.data.nodalConnec = [1:obj.data.numElements ; 2:obj.data.numElements+1]';
                obj.data.materialConnec = ones(obj.data.numElements,1);
                obj.data.dofsConnec = connectDOF(obj.data.numElements,obj.data.numDOFsNode,obj.data.nodalConnec);
                
                m = FEMBeamComputer(obj.data);
                m.computeGeoDiscret();
                m.computeSectionSolver();
                m.computemD2();
                m.computeForceMomentElem();
                m.computeBeamSolver();
                
                u_conv(:,nel) = m.u(end-2:end);
            end

            error_uy = zeros(1,size(u_conv,2)-1);
            error_rotx = zeros(1,size(u_conv,2)-1);
            error_rotz = zeros(1,size(u_conv,2)-1);
        
            for ii = 1:size(u_conv,2)-1
                error_uy(1,ii) = abs((u_conv(1,ii) - u_conv(1,end))/(u_conv(1,end)))*100;
                error_rotx(1,ii) = abs((u_conv(2,ii) - u_conv(2,end))/(u_conv(2,end)))*100;
                error_rotz(1,ii) = abs((u_conv(3,ii) - u_conv(3,end))/(u_conv(3,end)))*100;
            end
            
            figure(12)
            semilogx(el_conv(1:end-1),error_uy,el_conv(1:end-1),error_rotx,el_conv(1:end-1),error_rotz)
            legend('error in u_y','error in rot_x','error in rot_z')
            grid on
            xlabel('Beam''s number of elements' )
            ylabel('Error (%)')
            title(sprintf('Study of Convergence'))

        end

    end

end