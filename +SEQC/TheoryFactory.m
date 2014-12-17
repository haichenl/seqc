classdef TheoryFactory < handle
    
    methods (Static)
        
        function theory = Create(theoryType)
            import SEQC.EnumTheory ...
                SEQC.Cndo2 SEQC.Indo SEQC.Mndo SEQC.Am1 SEQC.Pm3;
            
            if(theoryType == EnumTheory.CNDO2)
                theory = Cndo2();
            elseif(theoryType == EnumTheory.INDO)
                theory = Indo();
            elseif(theoryType == EnumTheory.MNDO)
                theory = Mndo();
            elseif(theoryType == EnumTheory.AM1)
                theory = Am1();
            elseif(theoryType == EnumTheory.PM3)
                theory = Pm3();
            else
                throw(MException('TheoryFactory:GetInstance', 'Theory type not implemented yet.'));
            end
        end
        
    end
    
end