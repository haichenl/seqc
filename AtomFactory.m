classdef AtomFactory < handle
    
    methods (Static)
        
        function atom = Create(atomType, index, xyz, charge)
            if(atomType == EnumAtom.H)
                atom = AtomH(index);
            elseif(atomType == EnumAtom.C)
                atom = AtomC(index);
            elseif(atomType == EnumAtom.N)
                atom = AtomN(index);
            elseif(atomType == EnumAtom.O)
                atom = AtomO(index);
            elseif(atomType == EnumAtom.F)
                atom = AtomF(index);
            elseif(atomType == EnumAtom.Cl)
                atom = AtomCl(index);
            else
                throw(MException('AtomFactory:Create', 'Atom type not implemented yet.'));
            end
            atom.SetXyz(xyz);
            if(nargin>3)
                atom.SetCoreCharge(charge);
            end
        end
        
    end
    
end