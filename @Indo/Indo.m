classdef Indo < Cndo2
    
    properties (SetAccess = private)
        
        indoG1VecAtom;
        indoF2VecAtom;
        indoG1VecShell;
        indoF2VecShell;
        indoG1VecBasis;
        indoF2VecBasis;
        indoF0CoefficientVecBasis;
        indoG1CoefficientVecBasis;
        indoF2CoefficientVecBasis;
        
        indoCoulombMat;
        indoExchangeMat;
        betaS;
        gammaEC;
        
    end
    
    methods (Access = public)
        
        function obj = Indo()
            obj.theory = EnumTheory.INDO;
            obj.SetEnableAtomTypes();
        end
        
        function SetMolecule(obj, mol)
            SetMolecule@Cndo2(obj, mol);
            obj.indoF0CoefficientVecBasis = zeros(obj.nbf, 1);
            obj.indoG1CoefficientVecBasis = zeros(obj.nbf, 1);
            obj.indoF2CoefficientVecBasis = zeros(obj.nbf, 1);
            for i = 1:obj.nbf
                angularType = obj.mapBasis2AngularType(i);
                atom = obj.molecule.atomVect{obj.mapBasis2Atom(i)};
                if(angularType == 1) % s
                    obj.indoF0CoefficientVecBasis(i) = atom.paramPool.indoF0CoefficientS;
                    obj.indoG1CoefficientVecBasis(i) = atom.paramPool.indoG1CoefficientS;
                    obj.indoF2CoefficientVecBasis(i) = atom.paramPool.indoF2CoefficientS;
                elseif(angularType == 2) % py pz px
                    obj.indoF0CoefficientVecBasis(i) = atom.paramPool.indoF0CoefficientP;
                    obj.indoG1CoefficientVecBasis(i) = atom.paramPool.indoG1CoefficientP;
                    obj.indoF2CoefficientVecBasis(i) = atom.paramPool.indoF2CoefficientP;
                else
                    throw(MException('Indo:SetMolecule', 'Angular type wrong.'));
                end
            end
        end
        
    end
    
    methods %(Access = protected)
        
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
        
        % generate protected vectorization stuffs
        function Preiterations(obj)
            Preiterations@Cndo2(obj);
            obj.indoF2VecShell = zeros(obj.nshell, 1);
            for i = 1:obj.nshell
                atom = obj.molecule.atomVect{obj.mapShell2Atom(i)};
                if(obj.mapShell2AngularType(i) == 1) % s
                    valueF2 = 0;
                else
                    valueF2 = atom.paramPool.indoF2;
                end
                obj.indoF2VecShell(i) = valueF2;
            end
            obj.indoG1VecAtom = zeros(obj.natom, 1);
            obj.indoF2VecAtom = zeros(obj.natom, 1);
            for i = 1:obj.natom
                atom = obj.molecule.atomVect{i};
                obj.indoG1VecAtom(i) = atom.paramPool.indoG1;
                obj.indoF2VecAtom(i) = atom.paramPool.indoF2;
            end
            obj.indoG1VecShell = obj.indoG1VecAtom(obj.mapShell2Atom);
            obj.indoG1VecBasis = obj.indoG1VecAtom(obj.mapBasis2Atom);
            obj.indoF2VecBasis = obj.indoF2VecAtom(obj.mapBasis2Atom);
            
            F2Mat = diag(obj.indoF2VecShell);
            F2Mat = F2Mat(obj.mapBasis2Shell, obj.mapBasis2Shell);
            F2MatDiag = diag(diag(F2Mat));
            F2MatOffDiag = F2Mat - F2MatDiag;
            G1Mat = diag(obj.indoG1VecAtom);
            G1Mat = G1Mat(obj.mapShell2Atom, obj.mapShell2Atom);
            G1Mat = G1Mat - diag(diag(G1Mat));
            G1Mat = G1Mat(obj.mapBasis2Shell, obj.mapBasis2Shell);
            gammaSemiDiag = diag(diag(obj.gammaAB));
            gammaSemiDiag = gammaSemiDiag(obj.mapBasis2Atom, obj.mapBasis2Atom);
            obj.indoCoulombMat = gammaSemiDiag + 4/25 .* F2MatDiag - 2/25 .* F2MatOffDiag;
            obj.indoExchangeMat = 1/3 .* G1Mat + 3/25 .* F2MatOffDiag + diag(diag(obj.indoCoulombMat));
            
            betaS_ = obj.bondParamMatAtom .* obj.bondParamKMat .* 0.5;
            betaS_ = betaS_ - diag(diag(betaS_));
            betaS_ = betaS_(obj.mapBasis2Atom, obj.mapBasis2Atom);
            betaS_ = betaS_ .* obj.overlapAOs;
            obj.betaS = betaS_;
            
            gammaOff = obj.gammaAB - diag(diag(obj.gammaAB));
            gammaOff = gammaOff(obj.mapBasis2Atom, obj.mapBasis2Atom);
            gammaEC_ = 3.*obj.indoExchangeMat - obj.indoCoulombMat - gammaOff;
            gammaEC_ = gammaEC_ ./ 2;
            obj.gammaEC = gammaEC_;
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
        
        function fockdiag = GetFockDiag(obj)
            fockdiag = - obj.imuAmuVecBasis;
            gammaijdiag = diag(obj.gammaij);
            fockdiag = fockdiag - (obj.indoF0CoefficientVecBasis.*gammaijdiag ...
                +obj.indoG1CoefficientVecBasis.*obj.indoG1VecBasis...
                +obj.indoF2CoefficientVecBasis.*obj.indoF2VecBasis);
            
            diagDens = diag(obj.orbitalElectronPopulation);
            fockdiag = fockdiag + (obj.indoCoulombMat - 0.5.*obj.indoExchangeMat) * diagDens;
            
            temp = obj.atomicElectronPopulation - obj.coreChargeVecAtom;
            temp = (obj.gammaAB - diag(diag(obj.gammaAB))) * temp;
            fockdiag = fockdiag + temp(obj.mapBasis2Atom);
        end
        
        function offdiagFock = GetFockOffDiag(obj)
            offdiagFock = obj.gammaEC .* obj.orbitalElectronPopulation + obj.betaS;
            offdiagFock = offdiagFock - diag(diag(offdiagFock));
        end
        
        function fullH1 = GetH1(obj)
            % diagonal part
            diagH1 = - obj.imuAmuVecBasis...
                - (obj.indoF0CoefficientVecBasis.*diag(obj.gammaij)...
                +obj.indoG1CoefficientVecBasis.*obj.indoG1VecBasis...
                +obj.indoF2CoefficientVecBasis.*obj.indoF2VecBasis);
            
            temp = (obj.gammaAB - diag(diag(obj.gammaAB))) * (-obj.coreChargeVecAtom);
            diagH1 = diagH1 + temp(obj.mapBasis2Atom);
            
            % off diagonal part
            fullH1 = obj.betaS; % betaS already has diagonal terms = 0
            
            % full H1
            fullH1 = fullH1 + diag(diagH1);
        end
        
        function fullG = GetG(obj)
            % diagonal part
            diagG = (obj.indoCoulombMat - 0.5.*obj.indoExchangeMat)...
                * diag(obj.orbitalElectronPopulation);
            temp = (obj.gammaAB - diag(diag(obj.gammaAB)))... % zero out gammAB's diagonal
                * obj.atomicElectronPopulation;
            diagG = diagG + temp(obj.mapBasis2Atom);
            
            % off diagonal part
            fullG = obj.gammaEC .* obj.orbitalElectronPopulation;
            fullG = fullG - diag(diag(fullG));
            
            % full G
            fullG = fullG + diag(diagG);
        end
        
        % GetGuessH1() same as cndo/2

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
    
    methods %(Access = private)
        
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




