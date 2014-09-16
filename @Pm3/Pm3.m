classdef Pm3 < Am1
    
    methods (Access = public)
        
        function obj = Pm3()
            obj.theory = EnumTheory.PM3;
            obj.SetEnableAtomTypes();
        end
        
    end
    
    methods (Access = protected)
        
        function SetEnableAtomTypes(obj)
            obj.enableAtomTypes = {};
            obj.enableAtomTypes{end+1} = EnumAtom.H;
            obj.enableAtomTypes{end+1} = EnumAtom.C;
            obj.enableAtomTypes{end+1} = EnumAtom.N;
            obj.enableAtomTypes{end+1} = EnumAtom.O;
            obj.enableAtomTypes{end+1} = EnumAtom.F;
            obj.enableAtomTypes{end+1} = EnumAtom.S;
            obj.enableAtomTypes{end+1} = EnumAtom.Cl;
            obj.enableAtomTypes{end+1} = EnumAtom.Zn;
        end
        
        function value = AtomGetOrbitalExponent(~, atom, ~, orbitalType)
            if(orbitalType == 1)
                value = atom.paramPool.pm3OrbitalExponentS;
            elseif(orbitalType == 4 ||...
                    orbitalType == 2 ||...
                    orbitalType == 3)
                value = atom.paramPool.pm3OrbitalExponentP;
            else
                throw(MException('Pm3:AtomGetOrbitalExponent', 'Orbital type wrong.'));
            end
        end
        
        function value = AtomGetNddoAlpha(~, atom)
            value = atom.paramPool.pm3Alpha;
        end
        
        function value = AtomGetCoreIntegral(~, atom, orbital)
            if(orbital == 1)
                value = atom.paramPool.pm3CoreintegralS;
            elseif(orbital == 4 || orbital == 2 || orbital == 3)
                value = atom.paramPool.pm3CoreintegralP;
            else
                throw(MException('Pm3:AtomGetCoreIntegral', 'Orbital type wrong.'));
            end
        end
        
        function value = AtomGetBondingParameter(~, atom, orbital)
            if(orbital == 1)
                value = atom.paramPool.pm3BondingParameterS;
            elseif(( orbital == 4 ||...
                    orbital == 2 ||...
                    orbital == 3 ) )
                value = atom.paramPool.pm3BondingParameterP;
            else
                throw(MException('Pm3:AtomGetBondingParameter', 'Orbital type wrong.'));
            end
        end
        
        function valueVec = AtomGetCoreIntegralVec(~, atom)
            valueVec = [atom.paramPool.pm3CoreintegralS; ...
                atom.paramPool.pm3CoreintegralP];
        end
        
        function valueVec = AtomGetBondingParameterVec(~, atom)
            valueVec = [atom.paramPool.pm3BondingParameterS; ...
                atom.paramPool.pm3BondingParameterP];
        end
        
        function res = AtomGetNddoGss(~, atom)
            res = atom.paramPool.pm3Gss;
        end
        
        function res = AtomGetNddoGsp(~, atom)
            res = atom.paramPool.pm3Gsp;
        end
        
        function res = AtomGetNddoGpp(~, atom)
            res = atom.paramPool.pm3Gpp;
        end
        
        function res = AtomGetNddoGpp2(~, atom)
            res = atom.paramPool.pm3Gpp2;
        end
        
        function res = AtomGetNddoHsp(~, atom)
            res = atom.paramPool.pm3Hsp;
        end
        
        function res = AtomGetNddoHpp(~, atom)
            res =  0.5*(atom.paramPool.pm3Gpp - atom.paramPool.pm3Gpp2);
        end
        
        function res = AtomGetNddoDerivedParameterD(~, atom, multipole)
            switch(multipole)
                case 1
                    dIndex = 1;
                case {8, 9, 10}
                    dIndex = 2;
                case {2,3,4,5,6,7}
                    dIndex = 3;
                otherwise
                    throw(MException('Pm3:AtomGetNddoDerivedParameterD', 'Multipole type wrong.'));
            end
            res = atom.paramPool.pm3DerivedParameterD(dIndex);
        end
        
        function res = AtomGetNddoDerivedParameterRho(~, atom, multipole)
            switch(multipole)
                case 1
                    rhoIndex = 1;
                case {8, 9, 10}
                    rhoIndex = 2;
                case {2,3,4,5,6,7}
                    rhoIndex = 3;
                otherwise
                    throw(MException('Pm3:AtomGetNddoDerivedParameterRho', 'Multipole type wrong.'));
            end
            res = atom.paramPool.pm3DerivedParameterRho(rhoIndex);
        end
        
        function res = AtomGetNddoDerivedParameterDVec(~, atom)
            res = atom.paramPool.pm3DerivedParameterD;
        end
        
        function res = AtomGetNddoDerivedParameterRhoVec(~, atom)
            res = atom.paramPool.pm3DerivedParameterRho;
        end
        
        function res = AtomGetNddoParameterK(~, atom, kIndex)
            if(kIndex == 1 || kIndex == 2 || kIndex == 3 || kIndex == 4)
                res = atom.paramPool.pm3ParameterK(kIndex);
            else
                throw(MException('Pm3:AtomGetNddoParameterK', 'kIndex wrong.'));
            end
        end
        
        function res = AtomGetNddoParameterL(~, atom, lIndex)
            if(lIndex == 1 || lIndex == 2 || lIndex == 3 || lIndex == 4)
                res = atom.paramPool.pm3ParameterL(lIndex);
            else
                throw(MException('Pm3:AtomGetNddoParameterL', 'lIndex wrong.'));
            end
        end
        
        function res = AtomGetNddoParameterM(~, atom, mIndex)
            if(mIndex == 1 || mIndex == 2 || mIndex == 3 || mIndex == 4)
                res = atom.paramPool.pm3ParameterM(mIndex);
            else
                throw(MException('Pm3:AtomGetNddoParameterM', 'mIndex wrong.'));
            end
        end
        
    end
    
end