classdef AtomO < SEQC.Atom
    
    methods
        
        function obj = AtomO(ind)
            obj@SEQC.Atom(ind);
            obj.SetAtomicParameters();
        end
        
    end
    
    methods (Access = protected)
        
        function SetAtomicParameters(obj)
            import SEQC.Arguments SEQC.RealSphericalHarmonicsIndex SEQC.ParamPoolC SEQC.EnumAtom;
            
            obj.atomType = uint8(EnumAtom.O);
            obj.atomicMass = 15.9994*Arguments.GetInstance().GetGMolin2AU();
            obj.coreCharge = 6.0;
            obj.numberValenceElectrons = 6;
            obj.valenceShellType = 2; % lShell
            obj.nShell = 2;
            obj.valence = (1:4)';
            obj.lVec = zeros(length(obj.valence), 1);
            obj.mVec = zeros(length(obj.valence), 1);
            for i=1:length(obj.valence)
                obj.realSphericalHarmonicsIndices{i} = RealSphericalHarmonicsIndex(obj.valence(i));
                obj.lVec(i) = obj.realSphericalHarmonicsIndices{i}.l;
                obj.mVec(i) = obj.realSphericalHarmonicsIndices{i}.m;
            end
            obj.vdWCoefficient = 0.70*Arguments.GetInstance().GetJ2AU()...
                *power(Arguments.GetInstance().GetNm2AU(),6.0)...
                /Arguments.GetInstance().GetAvogadro();
            obj.vdWRadii = 1.490*Arguments.GetInstance().GetAngstrom2AU();
            
            obj.paramPool = ParamPoolO.GetInstance();
            
        end
        
    end
    
end