classdef Am1 < SEQC.Mndo
    
    methods (Access = public)
        
        function obj = Am1()
            obj.theory = SEQC.EnumTheory.AM1;
        end
        
    end
    
    methods (Access = protected)
        
        function SetEnableAtomTypes(obj)
            import SEQC.EnumAtom;
            
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
        
        function res = GetDiatomCoreRepulsionEnergy(obj, atomA, atomB)
            % MNDO term
            mndoTerm = GetDiatomCoreRepulsionEnergy@SEQC.Mndo(obj, atomA, atomB);
            
            % additional term, Eq. (4) in [S_1989].
            distance   = obj.molecule.GetDistanceAtoms(atomA, atomB);
            ang2AU     = SEQC.Arguments.GetInstance().GetAngstrom2AU();
            
            kA = obj.AtomGetNddoParameterKVec(atomA);
            lA = obj.AtomGetNddoParameterLVec(atomA);
            mA = obj.AtomGetNddoParameterMVec(atomA);
            temp = sum(obj.GetAdditionalDiatomCoreRepulsionTerm(kA, lA, mA, distance));
            
            kB = obj.AtomGetNddoParameterKVec(atomB);
            lB = obj.AtomGetNddoParameterLVec(atomB);
            mB = obj.AtomGetNddoParameterMVec(atomB);
            temp = temp + sum(obj.GetAdditionalDiatomCoreRepulsionTerm(kB, lB, mB, distance));

            additionalTerm = atomA.coreCharge*atomB.coreCharge*temp*ang2AU/distance;
            
            res = mndoTerm + additionalTerm;
        end
        
        %    virtual double GetDiatomCoreRepulsion1stDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                       const MolDS_base_atoms::Atom& atomB,
        %                                                       MolDS_base::CartesianType axisA) const;
        %    virtual double GetDiatomCoreRepulsion2ndDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                       const MolDS_base_atoms::Atom& atomB,
        %                                                       MolDS_base::CartesianType axisA1,
        %                                                       MolDS_base::CartesianType axisA2) const;
        %
        
        function value = AtomGetOrbitalExponent(~, atom, ~, orbitalType)
            if(orbitalType == 1)
                value = atom.paramPool.am1OrbitalExponentS;
            elseif(orbitalType == 4 ||...
                    orbitalType == 2 ||...
                    orbitalType == 3)
                value = atom.paramPool.am1OrbitalExponentP;
            else
                throw(MException('Am1:AtomGetOrbitalExponent', 'Orbital type wrong.'));
            end
        end
        
        function value = AtomGetNddoAlpha(~, atom)
            value = atom.paramPool.am1Alpha;
        end
        
        function value = AtomGetCoreIntegral(~, atom, orbital)
            if(orbital == 1)
                value = atom.paramPool.am1CoreintegralS;
            elseif(orbital == 4 || orbital == 2 || orbital == 3)
                value = atom.paramPool.am1CoreintegralP;
            else
                throw(MException('Am1:AtomGetCoreIntegral', 'Orbital type wrong.'));
            end
        end
        
        function value = AtomGetBondingParameter(~, atom, orbital)
            if(orbital == 1)
                value = atom.paramPool.am1BondingParameterS;
            elseif(( orbital == 4 ||...
                    orbital == 2 ||...
                    orbital == 3 ) )
                value = atom.paramPool.am1BondingParameterP;
            else
                throw(MException('Am1:AtomGetBondingParameter', 'Orbital type wrong.'));
            end
        end
        
        function valueVec = AtomGetCoreIntegralVec(~, atom)
            valueVec = [atom.paramPool.am1CoreintegralS; ...
                atom.paramPool.am1CoreintegralP];
        end
        
        function valueVec = AtomGetBondingParameterVec(~, atom)
            valueVec = [atom.paramPool.am1BondingParameterS; ...
                atom.paramPool.am1BondingParameterP];
        end
        
        function res = AtomGetNddoGss(~, atom)
            res = atom.paramPool.am1Gss;
        end
        
        function res = AtomGetNddoGsp(~, atom)
            res = atom.paramPool.am1Gsp;
        end
        
        function res = AtomGetNddoGpp(~, atom)
            res = atom.paramPool.am1Gpp;
        end
        
        function res = AtomGetNddoGpp2(~, atom)
            res = atom.paramPool.am1Gpp2;
        end
        
        function res = AtomGetNddoHsp(~, atom)
            res = atom.paramPool.am1Hsp;
        end
        
        function res = AtomGetNddoHpp(~, atom)
            res =  0.5*(atom.paramPool.am1Gpp - atom.paramPool.am1Gpp2);
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
                    throw(MException('Am1:AtomGetNddoDerivedParameterD', 'Multipole type wrong.'));
            end
            res = atom.paramPool.am1DerivedParameterD(dIndex);
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
                    throw(MException('Am1:AtomGetNddoDerivedParameterRho', 'Multipole type wrong.'));
            end
            res = atom.paramPool.am1DerivedParameterRho(rhoIndex);
        end
        
        function res = AtomGetNddoDerivedParameterDVec(~, atom)
            res = atom.paramPool.am1DerivedParameterD;
        end
        
        function res = AtomGetNddoDerivedParameterRhoVec(~, atom)
            res = atom.paramPool.am1DerivedParameterRho;
        end
        
        function res = AtomGetNddoParameterK(~, atom, kIndex)
            if(kIndex == 1 || kIndex == 2 || kIndex == 3 || kIndex == 4)
                res = atom.paramPool.am1ParameterK(kIndex);
            else
                throw(MException('Am1:AtomGetNddoParameterK', 'kIndex wrong.'));
            end
        end
        
        function res = AtomGetNddoParameterL(~, atom, lIndex)
            if(lIndex == 1 || lIndex == 2 || lIndex == 3 || lIndex == 4)
                res = atom.paramPool.am1ParameterL(lIndex);
            else
                throw(MException('Am1:AtomGetNddoParameterL', 'lIndex wrong.'));
            end
        end
        
        function res = AtomGetNddoParameterM(~, atom, mIndex)
            if(mIndex == 1 || mIndex == 2 || mIndex == 3 || mIndex == 4)
                res = atom.paramPool.am1ParameterM(mIndex);
            else
                throw(MException('Am1:AtomGetNddoParameterM', 'mIndex wrong.'));
            end
        end
        
        function res = AtomGetNddoParameterKVec(~, atom)
            res = atom.paramPool.am1ParameterK;
        end
        
        function res = AtomGetNddoParameterLVec(~, atom)
            res = atom.paramPool.am1ParameterL;
        end
        
        function res = AtomGetNddoParameterMVec(~, atom)
            res = atom.paramPool.am1ParameterM;
        end
        
    end
    
    methods (Access = private)
        
        function res = GetAdditionalDiatomCoreRepulsionTerm(~, k, l, m, distance)
            res =  k.*exp(-l.*power(distance-m,2.0));
        end
        
        function res =  GetAdditionalDiatomCoreRepulsionTerm1stDerivative(~, k, l, m, distance)
            res =  -2.0.*l*(distance-m).*k.*exp(-l.*power(distance-m,2.0));
        end
        
        function res = GetAdditionalDiatomCoreRepulsionTerm2ndDerivative(~, k, l, m, distance)
            res = (-2.0.*l + power(2.0.*l.*(distance-m),2.0)).*k.*exp(-l.*power(distance-m,2.0));
        end
        
    end
    
end