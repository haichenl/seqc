classdef ParamPoolManagerFactory < handle
    
    methods (Static)
        
        function paramPoolManager = Create()
            TheoryType = Arguments.GetInstance().currentTheory;
            if(TheoryType == EnumTheory.CNDO2)
                paramPoolManager = ParamPoolManagerCndo2.GetInstance();
            else
                throw(MException('AtomFactory:Create', 'Atom type not implemented yet.'));
            end
        end
        
    end
    
end