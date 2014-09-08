classdef AtomF < Atom
    
    methods
        
        function obj = AtomF(ind)
            obj@Atom(ind);
            obj.SetAtomicParameters();
        end
        
    end
    
    methods (Access = protected)
        
        function SetAtomicParameters(obj)
            obj.atomType = uint8(EnumAtom.F);
            obj.atomicMass = 18.9984032*Arguments.GetInstance().GetGMolin2AU();
            obj.coreCharge = 7.0;
            obj.numberValenceElectrons = 7;
            obj.valenceShellType = 2; % lShell
            obj.valence(1) = 1; % s
            obj.valence(2) = 2; % py
            obj.valence(3) = 3; % pz
            obj.valence(4) = 4; % px
            for i=1:length(obj.valence)
                obj.realSphericalHarmonicsIndices{i} = RealSphericalHarmonicsIndex(obj.valence(i));
            end
            obj.vdWCoefficient = 0.57*Arguments.GetInstance().GetJ2AU()...
                *power(Arguments.GetInstance().GetNm2AU(),6.0)...
                /Arguments.GetInstance().GetAvogadro();
            obj.vdWRadii = 1.430*Arguments.GetInstance().GetAngstrom2AU();
            
            obj.paramPool = ParamPoolF.GetInstance();
            
        end
        
    end
    
end