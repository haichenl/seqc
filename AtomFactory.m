classdef AtomFactory < handle
    
    methods (Static)
        
        function atom = Create(atomType, index, xyz, charge)
            if(atomType == AtomType.H)
                atom = HAtom(index);
            elseif(atomType == AtomType.C)
                atom = CAtom(index);
            elseif(atomType == AtomType.Cl)
                atom = ClAtom(index);
            else
                throw(MException('AtomFactory:Create', 'Atom type not implemented yet.'));
            end
            atom.SetXyz(xyz);
            if(nargin>3)
                atom.SetCoreCharge(charge);
            end
        end
        
    end
    
    methods (Access = private)
        
        function obj = AtomFactory()
        end
        
    end
    
end