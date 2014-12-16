classdef AtomF < SEQC.Atom
    
    methods
        
        function obj = AtomF(ind)
            obj@SEQC.Atom(ind);
            obj.SetAtomicParameters();
        end
        
    end
    
    methods (Access = protected)
        
        function SetAtomicParameters(obj)
            import SEQC.Arguments SEQC.RealSphericalHarmonicsIndex SEQC.ParamPoolF SEQC.EnumAtom;
            
            obj.atomType = uint8(EnumAtom.F);
            obj.atomicMass = 18.9984032*Arguments.GetInstance().GetGMolin2AU();
            obj.coreCharge = 7.0;
            obj.numberValenceElectrons = 7;
            obj.valenceShellType = 2; % lShell
            obj.nShell = 2;
            obj.valence = (1:4)';
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