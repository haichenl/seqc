classdef AtomO < Atom
    
    methods
        
        function obj = AtomO(ind)
            obj@Atom(ind);
            obj.SetAtomicParameters();
        end
        
    end
    
    methods (Access = protected)
        
        function SetAtomicParameters(obj)
            obj.atomType = EnumAtom.O;
            obj.atomicMass = 15.9994*Arguments.GetInstance().GetGMolin2AU();
            obj.coreCharge = 6.0;
            obj.numberValenceElectrons = 6;
            obj.valenceShellType = 2; % lShell
            obj.valence(1) = 1; % s
            obj.valence(2) = 2; % py
            obj.valence(3) = 3; % pz
            obj.valence(4) = 4; % px
            for i=1:length(obj.valence)
                obj.realSphericalHarmonicsIndices{i} = RealSphericalHarmonicsIndex(obj.valence(i));
            end
            obj.vdWCoefficient = 0.70*Arguments.GetInstance().GetJ2AU()...
                *power(Arguments.GetInstance().GetNm2AU(),6.0)...
                /Arguments.GetInstance().GetAvogadro();
            obj.vdWRadii = 1.490*Arguments.GetInstance().GetAngstrom2AU();
            
            obj.paramPool = ParamPoolO.GetInstance();
            
        end
        
    end
    
end