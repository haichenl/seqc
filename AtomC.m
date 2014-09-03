classdef AtomC < Atom
    
    methods
        
        function obj = AtomC(ind)
            obj@Atom(ind);
            obj.SetAtomicParameters();
        end
        
    end
    
    methods (Access = protected)
        
        function SetAtomicParameters(obj)
            obj.atomType = EnumAtom.C;
            obj.atomicMass = 12.0107*Arguments.GetInstance().GetGMolin2AU();
            obj.coreCharge = 4.0;
            obj.numberValenceElectrons = 4;
%             obj.valenceShellType = EnumShell.lShell;
            obj.valenceShellType = 2;
%             obj.valence{1} = EnumOrbital.s;
%             obj.valence{2} = EnumOrbital.py;
%             obj.valence{3} = EnumOrbital.pz;
%             obj.valence{4} = EnumOrbital.px;
            obj.valence(1) = 1;
            obj.valence(2) = 2;
            obj.valence(3) = 3;
            obj.valence(4) = 4;
            for i=1:length(obj.valence)
                obj.realSphericalHarmonicsIndices{i} = RealSphericalHarmonicsIndex(obj.valence(i));
            end
            obj.vdWCoefficient = 1.65*Arguments.GetInstance().GetJ2AU()...
                *power(Arguments.GetInstance().GetNm2AU(),6.0)...
                /Arguments.GetInstance().GetAvogadro();
            obj.vdWRadii = 1.610*Arguments.GetInstance().GetAngstrom2AU();
            
            obj.paramPool = ParamPoolC.GetInstance();
            
        end
        
    end
    
end