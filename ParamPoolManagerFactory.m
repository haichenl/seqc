classdef ParamPoolManagerFactory < handle
    
    methods (Static)
        
        function paramPoolManager = GetManagerInstance(theoryType)
            if(theoryType == EnumTheory.CNDO2)
                paramPoolManager = ParamPoolManagerCndo2.GetInstance();
            elseif(theoryType == EnumTheory.INDO)
                paramPoolManager = ParamPoolManagerIndo.GetInstance();
            else
                throw(MException('ParamPoolManagerFactory:GetManagerInstance', 'Theory type not implemented yet.'));
            end
        end
        
    end
    
end