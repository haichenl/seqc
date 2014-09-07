classdef Indo < Cndo2
    
    methods (Access = public)
        
        function obj = Indo()
        end
        
    end
    
    methods (Access = protected)
        
        function SetEnableAtomTypes(obj)
            obj.enableAtomTypes = {};
            obj.enableAtomTypes{end+1} = EnumAtom.H;
            obj.enableAtomTypes{end+1} = EnumAtom.Li;
            % obj.enableAtomTypes{end+1} = EnumAtom.Be;
            % obj.enableAtomTypes{end+1} = EnumAtom.B;
            obj.enableAtomTypes{end+1} = EnumAtom.C;
            obj.enableAtomTypes{end+1} = EnumAtom.N;
            obj.enableAtomTypes{end+1} = EnumAtom.O;
            obj.enableAtomTypes{end+1} = EnumAtom.F;
        end
        
        function value = GetFockDiagElement(obj, atomA, mu, isGuess)
            indexAtomA = atomA.index;
            firstAOIndexA = atomA.GetFirstAOIndex();
            value = obj.AtomGetCoreIntegral(atomA, atomA.valence(mu-firstAOIndexA+1), ...
                obj.gammaAB(indexAtomA,indexAtomA), isGuess);
            
            if(~isGuess)
                temp = 0.0;
                orbitalMu = atomA.valence(mu-firstAOIndexA+1);
                for v = 1:atomA.GetValenceSize()
                    orbitalLam = atomA.valence(v);
                    coulomb  = obj.GetCoulombInt(orbitalMu, orbitalLam, obj.gammaAB(indexAtomA,indexAtomA), atomA);
                    exchange = obj.GetExchangeInt(orbitalMu, orbitalLam, obj.gammaAB(indexAtomA,indexAtomA), atomA);
                    lammda = firstAOIndexA + v - 1;
                    temp = temp + obj.orbitalElectronPopulation(lammda,lammda)*(coulomb - 0.5*exchange);
                end
                value = value + temp;
                
                temp = 0.0;
                for B = 1:length(obj.molecule.atomVect)
                    if(B ~= indexAtomA)
                        atomB = obj.molecule.atomVect{B};
                        temp = temp + ( obj.atomicElectronPopulation(B) - atomB.coreCharge  )...
                            *obj.gammaAB(indexAtomA,B);
                    end
                end
                value = value + temp;
            end
        end
        
        function value = GetFockOffDiagElement(obj, atomA, atomB, mu, nu, isGuess)
            indexAtomA = atomA.index;
            indexAtomB = atomB.index;
            K = obj.GetBondingAdjustParameterK(atomA.valenceShellType, atomB.valenceShellType);
            bondParameter = 0.5*K*(obj.AtomGetBondingParameter(atomA) + obj.AtomGetBondingParameter(atomB));
            
            if(isGuess)
                value = bondParameter*obj.overlapAOs(mu,nu);
            else
                if(indexAtomA == indexAtomB)
                    orbitalMu = atomA.valence(mu-atomA.GetFirstAOIndex()+1);
                    orbitalNu = atomA.valence(nu-atomA.GetFirstAOIndex()+1);
                    coulomb  = obj.GetCoulombInt(orbitalMu, orbitalNu, obj.gammaAB(indexAtomA,indexAtomA), atomA);
                    exchange = obj.GetExchangeInt(orbitalMu, orbitalNu, obj.gammaAB(indexAtomA,indexAtomA), atomA);
                    value = (1.5*exchange - 0.5*coulomb)*obj.orbitalElectronPopulation(mu,nu);
                else
                    value = bondParameter*obj.overlapAOs(mu,nu);
                    value = value - 0.5*obj.orbitalElectronPopulation(mu,nu)*obj.gammaAB(indexAtomA,indexAtomB);
                end
            end
        end

        function value = GetMolecularIntegralElement(obj, moI, moJ, moK, moL)
            
            % CNDO terms
            value = GetMolecularIntegralElement@Cndo2(moI, moJ, moK, moL);
            
            % Aditional terms for INDO, see Eq. (10) in (RZ_1973)
            for A = 1:length(obj.molecule.atomVect)
                atomA = obj.molecule.atomVect{A};
                firstAOIndexA = atomA.GetFirstAOIndex();
                numberAOsA = atomA.GetValenceSize();
                
                for mu = firstAOIndexA:firstAOIndexA+numberAOsA-1
                    orbitalMu = atomA.valence(mu-firstAOIndexA+1);
                    for nu = firstAOIndexA:firstAOIndexA+numberAOsA-1
                        orbitalNu = atomA.valence(nu-firstAOIndexA+1);
                        
                        if(mu~=nu)
                            exchange = obj.GetExchangeInt(orbitalMu, orbitalNu, obj.gammaAB(A,A), atomA);
                            value = value + exchange...
                                *obj.fockMatrix(moI,mu)...
                                *obj.fockMatrix(moJ,nu)...
                                *obj.fockMatrix(moK,nu)...
                                *obj.fockMatrix(moL,mu);
                            value = value + exchange...
                                *obj.fockMatrix(moI,mu)...
                                *obj.fockMatrix(moJ,nu)...
                                *obj.fockMatrix(moK,mu)...
                                *obj.fockMatrix(moL,nu);
                        end
                        
                        coulomb = obj.GetCoulombInt(orbitalMu, orbitalNu, gammaAB(A,A), atomA);
                        value = value + (coulomb-obj.gammaAB(A,A))...
                            *obj.fockMatrix(moI,mu)...
                            *obj.fockMatrix(moJ,mu)...
                            *obj.fockMatrix(moK,nu)...
                            *obj.fockMatrix(moL,nu);
                    end
                end
            end
        end
        
        function res = AtomGetBondingParameter(~, atom)
            res = atom.paramPool.bondingParameter;
        end
        
        function value = AtomGetCoreIntegral(~, atom, orbital, gamma, isGuess)
            if(orbital == 1)
                value = -1.0*atom.paramPool.imuAmuS;
                if(~isGuess)
                    value = value - (atom.paramPool.indoF0CoefficientS*gamma ...
                        +atom.paramPool.indoG1CoefficientS*atom.paramPool.indoG1...
                        +atom.paramPool.indoF2CoefficientS*atom.paramPool.indoF2);
                end
            elseif(orbital == 4 || orbital == 2 || orbital == 3)
                value = -1.0*atom.paramPool.imuAmuP;
                if(~isGuess)
                    value = value - (atom.paramPool.indoF0CoefficientP*gamma ...
                        +atom.paramPool.indoG1CoefficientP*atom.paramPool.indoG1...
                        +atom.paramPool.indoF2CoefficientP*atom.paramPool.indoF2);
                end
            else
                throw(MException('Indo:AtomGetIndoCoreIntegral', 'Orbital type wrong.'));
            end
        end

    end
    
    methods (Access = private)
        
        % (3.87) - (3.91) in J. A. Pople book.
        % Indo Coulomb Interaction
        function value = GetCoulombInt(~, orbital1, orbital2, gamma, atom)
            if( orbital1 == 1 && orbital2 == 1)
                value = gamma;
            elseif( orbital1 == 1 && ( orbital2 == 4 || orbital2 == 2 || orbital2 == 3 ))
                value = gamma;
            elseif( (orbital1 == 4 || orbital1 == 2 || orbital1 == 3 ) && orbital2 == 1)
                value = gamma;
            elseif( (orbital1 == orbital2) && ( orbital1 == 4 || orbital1 == 2 || orbital1 == 3 ))
                value = gamma + 4.0*atom.paramPool.indoF2/25.0;
            elseif( (orbital1 ~= orbital2) ...
                    && ( orbital1 == 4 || orbital1 == 2 || orbital1 == 3 )...
                    && ( orbital2 == 4 || orbital2 == 2 || orbital2 == 3 ) )
                value = gamma - 2.0*atom.paramPool.indoF2/25.0;
            else
                throw(MException('Indo:GetCoulombInt', 'Orbital type wrong.'));
            end
        end

        function value = GetExchangeInt(obj, orbital1, orbital2, gamma, atom)
            if( orbital1 == orbital2)
                value = obj.GetCoulombInt(orbital1, orbital2, gamma, atom);
            elseif( (orbital1 == 1) && (orbital2 == 4 || orbital2 == 2 || orbital2 == 3 ) )
                value = atom.paramPool.indoG1/3.0;
            elseif( (orbital1 == 4 || orbital1 == 2 || orbital1 == 3) && orbital2 == 1  )
                value = atom.paramPool.indoG1/3.0;
            elseif( (orbital1 ~= orbital2) ...
                    && ( orbital1 == 4 || orbital1 == 2 || orbital1 == 3 )...
                    && ( orbital2 == 4 || orbital2 == 2 || orbital2 == 3 ) )
                value = 3.0*atom.paramPool.indoF2/25.0;
            else
                throw(MException('Indo:GetExchangeInt', 'Orbital type wrong.'));
            end
        end


    end
    
end




