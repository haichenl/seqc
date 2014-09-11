classdef ParamPoolManagerFactory < handle
    
    methods (Static)
        
        function paramPoolManager = GetManagerInstance(theoryType)
            if(theoryType == EnumTheory.CNDO2)
                paramPoolManager = ParamPoolManagerCndo2.GetInstance();
            elseif(theoryType == EnumTheory.INDO)
                paramPoolManager = ParamPoolManagerIndo.GetInstance();
            elseif(theoryType == EnumTheory.MNDO)
                paramPoolManager = ParamPoolManagerMndo.GetInstance();
            elseif(theoryType == EnumTheory.AM1)
                paramPoolManager = ParamPoolManagerAm1.GetInstance();
            elseif(theoryType == EnumTheory.PM3)
                paramPoolManager = ParamPoolManagerPm3.GetInstance();
            else
                throw(MException('ParamPoolManagerFactory:GetManagerInstance', 'Theory type not implemented yet.'));
            end
        end
        
    end
    
end