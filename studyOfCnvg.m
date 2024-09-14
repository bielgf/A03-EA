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
            
            elConv = [4 8 16 32 64 128 256 512];
            uConv = zeros(3,size(elConv,2));

            for nel = 1:size(elConv,2)
                obj.data.numElements = elConv(1,nel);
                obj.data.nodalConnec = [1:obj.data.numElements ; 2:obj.data.numElements+1]';
                obj.data.materialConnec = ones(obj.data.numElements,1);
                obj.data.dofsConnec = connectDOF(obj.data.numElements,obj.data.numDOFsNode,obj.data.nodalConnec);
                
                m = FEMBeamComputer(obj.data);
                m.compute();
                
                uConv(:,nel) = m.u(end-2:end);
            end

            error_uy = zeros(1,size(uConv,2)-1);
            error_rotx = zeros(1,size(uConv,2)-1);
            error_rotz = zeros(1,size(uConv,2)-1);
        
            for ii = 1:size(uConv,2)-1
                error_uy(1,ii) = abs((uConv(1,ii) - uConv(1,end))/(uConv(1,end)))*100;
                error_rotx(1,ii) = abs((uConv(2,ii) - uConv(2,end))/(uConv(2,end)))*100;
                error_rotz(1,ii) = abs((uConv(3,ii) - uConv(3,end))/(uConv(3,end)))*100;
            end
            
            figure(12)
            semilogx(elConv(1:end-1),error_uy,elConv(1:end-1),error_rotx,elConv(1:end-1),error_rotz)
            legend('error in u_y','error in rot_x','error in rot_z')
            grid on
            xlabel('Beam''s number of elements' )
            ylabel('Error (%)')
            title(sprintf('Study of Convergence'))

        end

    end

end