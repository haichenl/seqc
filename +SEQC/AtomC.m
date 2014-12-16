classdef AtomC < SEQC.Atom
    
    methods
        
        function obj = AtomC(ind)
            obj@SEQC.Atom(ind);
            obj.SetAtomicParameters();
        end
        
    end
    
    methods (Access = protected)
        
        function SetAtomicParameters(obj)
            import SEQC.Arguments SEQC.RealSphericalHarmonicsIndex SEQC.ParamPoolC SEQC.EnumAtom;
            
            obj.atomType = uint8(EnumAtom.C);
            obj.atomicMass = 12.0107*Arguments.GetInstance().GetGMolin2AU();
            obj.coreCharge = 4.0;
            obj.numberValenceElectrons = 4;
%             obj.valenceShellType = EnumShell.lShell;
            obj.valenceShellType = 2;
            obj.nShell = 2;
            obj.valence = (1:4)';
            obj.lVec = zeros(length(obj.valence), 1);
            obj.mVec = zeros(length(obj.valence), 1);
            for i=1:length(obj.valence)
                obj.realSphericalHarmonicsIndices{i} = RealSphericalHarmonicsIndex(obj.valence(i));
                obj.lVec(i) = obj.realSphericalHarmonicsIndices{i}.l;
                obj.mVec(i) = obj.realSphericalHarmonicsIndices{i}.m;
            end
            obj.vdWCoefficient = 1.65*Arguments.GetInstance().GetJ2AU()...
                *power(Arguments.GetInstance().GetNm2AU(),6.0)...
                /Arguments.GetInstance().GetAvogadro();
            obj.vdWRadii = 1.610*Arguments.GetInstance().GetAngstrom2AU();
            
            obj.paramPool = ParamPoolC.GetInstance();
            
        end
        
    end
    
end