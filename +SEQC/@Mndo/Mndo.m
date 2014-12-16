classdef Mndo < SEQC.Cndo2
    
    properties (Access = protected)
        
        nddoCoreIntVecBasis;
        bondParamSumMatBasis;
        nddoCoulombMat;
        nddoExchangeMat
        
        twoElecsCache1;
        twoElecsCache2;
        
    end
    
    properties (Access = private)
        
        twoElecsTwoAtomCoresMpiBuff;
        twoElecsAtomEpcCoresMpiBuff;
        heatsFormation = 0;
        
    end
    
    methods (Access = public)
        
        function obj = Mndo()
            % protedted variables and methods
            obj.theory = SEQC.EnumTheory.MNDO;
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
        
        % generate protected vectorization stuffs
        function Preiterations(obj)
            nddoCoreIntVecShell = zeros(obj.nshell, 1);
            bondParamVecShell = zeros(obj.nshell, 1);
            for i = 1:obj.natom
                atom = obj.molecule.atomVect{i};
                coreIntVec = obj.AtomGetCoreIntegralVec(atom);
                dimShell = atom.GetFirstShellIndex():atom.GetLastShellIndex();
                nddoCoreIntVecShell(dimShell) = coreIntVec(1:1+atom.nShell-1);
                
                bondParamVec = obj.AtomGetBondingParameterVec(atom);
                bondParamVecShell(dimShell) = bondParamVec(1:1+atom.nShell-1);
                
            end
            obj.nddoCoreIntVecBasis = nddoCoreIntVecShell(obj.mapBasis2Shell);
            bondParamMatShell = bondParamVecShell(:,ones(1,obj.nshell));
            bondParamMatShell = bondParamMatShell + bondParamMatShell';
            obj.bondParamSumMatBasis = bondParamMatShell(obj.mapBasis2Shell,obj.mapBasis2Shell);
            obj.bondParamSumMatBasis = obj.bondParamSumMatBasis - diag(diag(obj.bondParamSumMatBasis));
            
            for i = 1:obj.natom
                atom = obj.molecule.atomVect{i};
                dimBlock = atom.GetFirstAOIndex():atom.GetLastAOIndex();
                obj.nddoCoulombMat(dimBlock,dimBlock) = obj.GetCoulombIntBlock(atom);
                obj.nddoExchangeMat(dimBlock,dimBlock) = obj.GetExchangeIntBlock(atom);
            end
            
            obj.twoElecsCache1 = cell(obj.natom);
            obj.twoElecsCache2 = cell(obj.natom);
            for A = 1:obj.natom
                for B = 1:obj.natom
                    atomA = obj.molecule.atomVect{A};
                    valLengthA = atomA.valence(end);
                    indexScaleA = 1:valLengthA;
                    twoElecsCache_ = reshape(obj.twoElecsTwoAtomCores(A,B,:,:,:,:), 4,4,4,4);
                    twoElecsCache_ = permute(twoElecsCache_, [3 2 4 1]);
                    obj.twoElecsCache1{A,B} = reshape(twoElecsCache_, 16, 16);
                    obj.twoElecsCache2{A,B} = reshape(obj.twoElecsTwoAtomCores(A,B,indexScaleA,indexScaleA,:,:),valLengthA^2,[]);
                end
            end
            
        end
        
        function gammaAB = CalcGammaAB(~)
            gammaAB = [];
        end
        
        %    virtual void CalcSCFProperties();
        %    virtual void CalcNormalModes(double** normalModes, double* normalForceConstants, const MolDS_base::Molecule& molecule) const;
        %    virtual void CalcForce(const std::vector<int>& elecStates);
        
        function res = GetAtomCoreEpcCoulombEnergy(obj, atom, epc)
            distance = obj.molecule.GetDistanceAtomEpc(atom.index, epc.index);
            res = atom.coreCharge()*epc.coreCharge()/distance;
        end
        
        function res = GetDiatomCoreRepulsionEnergy(obj, atomA, atomB)
            tmp = obj.GetAuxiliaryDiatomCoreRepulsionEnergy(atomA, atomB, ...
                obj.molecule.GetDistanceAtoms(atomA, atomB));
            res = atomA.coreCharge...
                *atomB.coreCharge...
                *obj.twoElecsTwoAtomCores(atomA.index,atomB.index,1,1,1,1)...
                *tmp;
        end
        
        %    virtual double GetDiatomCoreRepulsion1stDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                       const MolDS_base_atoms::Atom& atomB,
        %                                                       MolDS_base::CartesianType axisA) const;
        %    virtual double GetDiatomCoreRepulsion2ndDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                       const MolDS_base_atoms::Atom& atomB,
        %                                                       MolDS_base::CartesianType axisA1,
        %                                                       MolDS_base::CartesianType axisA2) const;
        
        function guessH1 = GetGuessH1(obj)
            fockDiag = obj.nddoCoreIntVecBasis;
            guessH1 = 0.5 .* obj.bondParamSumMatBasis .* obj.overlapAOs;
            guessH1 = guessH1 - diag(diag(guessH1)) + diag(fockDiag);
        end
        
        function fullH1 = GetH1(obj)
            marginizer = ones(obj.natom);
            marginizer = marginizer - diag(diag(marginizer));
            marginizer = marginizer(obj.mapBasis2Atom,obj.mapBasis2Atom);
            marginH1 = 0.5 .* marginizer .* obj.bondParamSumMatBasis .* obj.overlapAOs;
            marginH1 = triu(marginH1, +1);
            marginH1 = marginH1 + marginH1';
            
            diagH1 = obj.nddoCoreIntVecBasis;
            centerH1 = diag(diagH1);
            
            for A = 1:obj.natom
                atomA = obj.molecule.atomVect{A};
                valLengthA = atomA.valence(end);
                indexScaleA = 1:valLengthA;
                temp = zeros(valLengthA^2, 1);
                for B = 1:obj.natom
                    atomB = obj.molecule.atomVect{B};
                    coreAttr = -atomB.coreCharge.*reshape(obj.twoElecsTwoAtomCores(A,B,indexScaleA,indexScaleA,1,1), [], 1);
                    if(B ~= A)
                        temp = temp + coreAttr;
                    end
                end
                ind = indexScaleA+atomA.firstAOIndex-1;
                centerH1(ind, ind) = centerH1(ind, ind) + reshape(temp,valLengthA,valLengthA);
            end
            
            
            fullH1 = marginH1 + centerH1;
        end
        
        function fullG = GetG(obj)
            marginG = zeros(obj.nbf);
            for A = 1:obj.natom
                valLengthA = obj.atomValLength(A);
                indScaleA = obj.atomAOinds{A};
                for B = 1:A-1
                    valLengthB = obj.atomValLength(B);
                    indScaleB = obj.atomAOinds{B};
                    twoElecsCache_ = obj.twoElecsCache1{A,B};
                    linedDens = reshape(obj.orbitalElectronPopulation(indScaleA,indScaleB)',[],1);
                    twoElecsCache_ = twoElecsCache_(1:length(linedDens), 1:length(linedDens));
                    marginG(indScaleB, indScaleA) = marginG(indScaleB, indScaleA) - 0.5 .* reshape(twoElecsCache_ * linedDens, valLengthB, valLengthA);
                end
            end
            marginG = triu(marginG, +1);
            marginG = marginG + marginG';
            
            diagDens = diag(obj.orbitalElectronPopulation);
            diagG = (obj.nddoCoulombMat - 0.5.*obj.nddoExchangeMat) * diagDens;
            centerG = (1.5.*obj.nddoExchangeMat - 0.5.*obj.nddoCoulombMat) .* obj.orbitalElectronPopulation;
            centerG = centerG - diag(diag(centerG)) + diag(diagG);
            
            for A = 1:obj.natom
                valLengthA = obj.atomValLength(A);
                temp = zeros(valLengthA^2, 1);
                for B = 1:obj.natom
                    indScaleB = obj.atomAOinds{B};
                    if(B ~= A)
                        linedDens = reshape(obj.orbitalElectronPopulation(indScaleB,indScaleB),[],1);
                        twoElecsMat = obj.twoElecsCache2{A,B};
                        twoElecsMat = twoElecsMat(:,1:length(linedDens));
                        temp = temp + twoElecsMat * linedDens;
                    end
                end
                ind = obj.atomAOinds{A};
                centerG(ind, ind) = centerG(ind, ind) + reshape(temp,valLengthA,valLengthA);
            end
            fullG = marginG + centerG;
        end
        
        function diatomicOverlapAOs = CalcDiatomicOverlapAOsInDiatomicFrame(obj, atomA, atomB)
            diatomicOverlapAOs = CalcDiatomicOverlapAOsInDiatomicFrame@SEQC.Cndo2(obj, atomA, atomB);
        end
        
        %    virtual void CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(double** diatomicOverlapAOsDeri,
        %                                                                    const MolDS_base_atoms::Atom& atomA,
        %                                                                    const MolDS_base_atoms::Atom& atomB) const;
        %    virtual void CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(double** diatomicOverlapAOs2ndDeri,
        %                                                                    const MolDS_base_atoms::Atom& atomA,
        %                                                                    const MolDS_base_atoms::Atom& atomB) const;
        
        % NDDO Coulomb interaction
        function value = GetCoulombInt(obj, orbital1, orbital2, atom)
            if( orbital1 == 1 && orbital2 == 1)
                value = obj.AtomGetNddoGss(atom);
            elseif( orbital1 == 1 && ( orbital2 == 4 || orbital2 == 2 || orbital2 == 3 ))
                value = obj.AtomGetNddoGsp(atom);
            elseif( orbital2 == 1 && ( orbital1 == 4 || orbital1 == 2 || orbital1 == 3 ))
                value = obj.GetCoulombInt(orbital2, orbital1, atom);
            elseif( (orbital1 == orbital2) && ( orbital1 == 4 || orbital1 == 2 || orbital1 == 3 ))
                value = obj.AtomGetNddoGpp(atom);
            elseif( (orbital1 ~= orbital2) ...
                    && ( orbital1 == 4 || orbital1 == 2 || orbital1 == 3 )...
                    && ( orbital2 == 4 || orbital2 == 2 || orbital2 == 3 ) )
                value = obj.AtomGetNddoGpp2(atom);
            else
                throw(MException('Mndo:GetCoulombInt', 'Orbital type wrong.'));
            end
        end
        
        % NDDO Exchange Interaction
        function value = GetExchangeInt(obj, orbital1, orbital2, atom)
            if( orbital1 == orbital2)
                value = obj.GetCoulombInt(orbital1, orbital2, atom);
            elseif( orbital1 == 1 && (orbital2 == 4 || orbital2 == 2 || orbital2 == 3 ) )
                value = obj.AtomGetNddoHsp(atom);
            elseif( orbital2 == 1 && (orbital1 == 4 || orbital1 == 2 || orbital1 == 3 ) )
                value = obj.GetExchangeInt(orbital2, orbital1, atom);
            elseif( (orbital1 ~= orbital2)...
                    && ( orbital1 == 4 || orbital1 == 2 || orbital1 == 3 )...
                    && ( orbital2 == 4 || orbital2 == 2 || orbital2 == 3 ) )
                value = obj.AtomGetNddoHpp(atom);
            else
                throw(MException('Mndo:GetExchangeInt', 'Orbital type wrong.'));
            end
        end
        
        function couBlock = GetCoulombIntBlock(obj, atom)
            couBlock = zeros(length(atom.valence));
            Gss = obj.AtomGetNddoGss(atom);
            couBlock(1,1) = Gss;
            if(length(atom.valence) == 1)
                return;
            end
            % deal with only p not d for now
            Gsp = obj.AtomGetNddoGsp(atom);
            couBlock(2:4,1) = Gsp;
            couBlock(1,2:4) = Gsp;
            couBlock(2:4,2:4) = obj.AtomGetNddoGpp2(atom);
            Gpp = obj.AtomGetNddoGpp(atom);
            diagBlock = [Gss;Gpp;Gpp;Gpp];
            couBlock = couBlock - diag(diag(couBlock)) + diag(diagBlock);
        end
        
        function excBlock = GetExchangeIntBlock(obj, atom)
            excBlock = zeros(length(atom.valence));
            Gss = obj.AtomGetNddoGss(atom);
            excBlock(1,1) = Gss;
            if(length(atom.valence) == 1)
                return;
            end
            % deal with only p not d for now
            Hsp = obj.AtomGetNddoHsp(atom);
            excBlock(2:4,1) = Hsp;
            excBlock(1,2:4) = Hsp;
            excBlock(2:4,2:4) = obj.AtomGetNddoHpp(atom);
            Gpp = obj.AtomGetNddoGpp(atom);
            diagblock = [Gss;Gpp;Gpp;Gpp];
            excBlock = excBlock - diag(diag(excBlock)) + diag(diagblock);
        end
        
        function CalcTwoElecsTwoCores(obj)
            obj.CalcTwoElecsTwoAtomCores();
            obj.CalcTwoElecsAtomEpcCores();
        end
        
        %    virtual double GetMolecularIntegralElement(int moI,
        %                                               int moJ,
        %                                               int moK,
        %                                               int moL,
        %                                               const MolDS_base::Molecule& molecule,
        %                                               double const* const* fockMatrix,
        %                                               double const* const* gammaAB) const;
        %    virtual void CalcCISMatrix(double** matrixCIS) const;
        %    virtual double GetSmallQElement(int moI,
        %                                    int moP,
        %                                    double const* const* xiOcc,
        %                                    double const* const* xiVir,
        %                                    double const* const* eta) const;
        %    virtual double GetAuxiliaryKNRKRElement(int moI, int moJ, int moK, int moL) const;
        
        function value = AtomGetOrbitalExponent(~, atom, ~, orbitalType)
            if(orbitalType == 1)
                value = atom.paramPool.mndoOrbitalExponentS;
            elseif(orbitalType == 4 ||...
                    orbitalType == 2 ||...
                    orbitalType == 3)
                value = atom.paramPool.mndoOrbitalExponentP;
            else
                throw(MException('Mndo:AtomGetOrbitalExponent', 'Orbital type wrong.'));
            end
        end
        
        function value = AtomGetNddoAlpha(~, atom)
            value = atom.paramPool.mndoAlpha;
        end
        
        function value = AtomGetCoreIntegral(~, atom, orbital)
            if(orbital == 1)
                value = atom.paramPool.mndoCoreintegralS;
            elseif(orbital == 4 || orbital == 2 || orbital == 3)
                value = atom.paramPool.mndoCoreintegralP;
            else
                throw(MException('Mndo:AtomGetCoreIntegral', 'Orbital type wrong.'));
            end
        end
        
        function value = AtomGetBondingParameter(~, atom, orbital)
            if(orbital == 1)
                value = atom.paramPool.mndoBondingParameterS;
            elseif(( orbital == 4 ||...
                    orbital == 2 ||...
                    orbital == 3 ) )
                value = atom.paramPool.mndoBondingParameterP;
            else
                throw(MException('Mndo:AtomGetBondingParameter', 'Orbital type wrong.'));
            end
        end
        
        function valueVec = AtomGetCoreIntegralVec(~, atom)
            valueVec = [atom.paramPool.mndoCoreintegralS; ...
                atom.paramPool.mndoCoreintegralP];
        end
        
        function valueVec = AtomGetBondingParameterVec(~, atom)
            valueVec = [atom.paramPool.mndoBondingParameterS; ...
                atom.paramPool.mndoBondingParameterP];
        end
        
        function res = AtomGetNddoGss(~, atom)
            res = atom.paramPool.mndoGss;
        end
        
        function res = AtomGetNddoGsp(~, atom)
            res = atom.paramPool.mndoGsp;
        end
        
        function res = AtomGetNddoGpp(~, atom)
            res = atom.paramPool.mndoGpp;
        end
        
        function res = AtomGetNddoGpp2(~, atom)
            res = atom.paramPool.mndoGpp2;
        end
        
        function res = AtomGetNddoHsp(~, atom)
            res = atom.paramPool.mndoHsp;
        end
        
        function res = AtomGetNddoHpp(~, atom)
            res =  0.5*(atom.paramPool.mndoGpp - atom.paramPool.mndoGpp2);
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
                    throw(MException('Mndo:AtomGetNddoDerivedParameterD', 'Multipole type wrong.'));
            end
            res = atom.paramPool.mndoDerivedParameterD(dIndex);
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
                    throw(MException('Mndo:AtomGetNddoDerivedParameterRho', 'Multipole type wrong.'));
            end
            res = atom.paramPool.mndoDerivedParameterRho(rhoIndex);
        end
        
        function res = AtomGetNddoDerivedParameterDVec(~, atom)
            res = atom.paramPool.mndoDerivedParameterD;
        end
        
        function res = AtomGetNddoDerivedParameterRhoVec(~, atom)
            res = atom.paramPool.mndoDerivedParameterRho;
        end
        
    end
    
    methods (Access = private)
        
%         function CalcTwoElecsTwoAtomCores(obj)
%             dxy = 4;
%             totalNumberAtoms = length(obj.molecule.atomVect);
%             obj.twoElecsTwoAtomCores = zeros(totalNumberAtoms,totalNumberAtoms,dxy,dxy,dxy,dxy);
%             
%             for a = 1:totalNumberAtoms
%                 % note that terms with condition a==b are not needed to calculate.
%                 for b = a+1:totalNumberAtoms
%                     diatomicTwoElecsTwoCores = ...
%                         obj.CalcDiatomicTwoElecsTwoCores( ...
%                         obj.molecule.atomVect{a}, obj.molecule.atomVect{b});
%                     i=1;
%                     for mu = 1:dxy
%                         for nu = mu:dxy
%                             j=1;
%                             for lambda = 1:dxy
%                                 for sigma = lambda:dxy
%                                     obj.twoElecsTwoAtomCoresMpiBuff(a,b,i,j) ...
%                                         = diatomicTwoElecsTwoCores(mu,nu,lambda,sigma);
%                                     j = j + 1;
%                                 end
%                             end
%                             i = i + 1;
%                         end
%                     end
%                 end
%             end
%             
%             for a = 1:totalNumberAtoms
%                 for b = a+1:totalNumberAtoms
%                     i=1;
%                     for mu = 1:dxy
%                         for nu = mu:dxy
%                             j=1;
%                             for lambda = 1:dxy
%                                 for sigma = lambda:dxy
%                                     value = obj.twoElecsTwoAtomCoresMpiBuff(a,b,i,j);
%                                     obj.twoElecsTwoAtomCores(a,b,mu,nu,lambda,sigma) = value;
%                                     obj.twoElecsTwoAtomCores(a,b,mu,nu,sigma,lambda) = value;
%                                     obj.twoElecsTwoAtomCores(a,b,nu,mu,lambda,sigma) = value;
%                                     obj.twoElecsTwoAtomCores(a,b,nu,mu,sigma,lambda) = value;
%                                     obj.twoElecsTwoAtomCores(b,a,lambda,sigma,mu,nu) = value;
%                                     obj.twoElecsTwoAtomCores(b,a,lambda,sigma,nu,mu) = value;
%                                     obj.twoElecsTwoAtomCores(b,a,sigma,lambda,mu,nu) = value;
%                                     obj.twoElecsTwoAtomCores(b,a,sigma,lambda,nu,mu) = value;
%                                     j = j + 1;
%                                 end
%                             end
%                             i = i + 1;
%                         end
%                     end
%                 end
%             end
%         end
        
        function CalcTwoElecsTwoAtomCores(obj)
            dxy = 4;
            totalNumberAtoms = length(obj.molecule.atomVect);
            obj.twoElecsTwoAtomCores = zeros(totalNumberAtoms,totalNumberAtoms,dxy,dxy,dxy,dxy);
            
            for a = 1:totalNumberAtoms
                % note that terms with condition a==b are not needed to calculate.
                for b = a+1:totalNumberAtoms
                    diatomicTwoElecsTwoCores = ...
                        obj.CalcDiatomicTwoElecsTwoCores( ...
                        obj.molecule.atomVect{a}, obj.molecule.atomVect{b});
                    obj.twoElecsTwoAtomCores(a,b,:,:,:,:) = diatomicTwoElecsTwoCores;
                    obj.twoElecsTwoAtomCores(b,a,:,:,:,:) = permute(diatomicTwoElecsTwoCores, [3 4 1 2]);
                end
            end
            
        end
        
        function CalcTwoElecsAtomEpcCores(obj)
            if(isempty(obj.molecule.epcVect))
                return;
            end
            dxy = 4;
            totalNumberAtoms = length(obj.molecule.atomVect);
            totalNumberEpcs  = length(obj.molecule.epcVect);
            obj.twoElecsAtomEpcCores = zeros(totalNumberAtoms, totalNumberEpcs,...
                dxy, dxy, dxy, dxy);
            
            for a = 1:totalNumberAtoms
                atom = obj.molecule.atomVect{a};
                % note that terms with condition a==b are not needed to calculate. spring: ?
                for b = 1:totalNumberEpcs
                    epc  = obj.molecule.epcVect{b};
                    diatomicTwoElecsTwoCores = obj.CalcDiatomicTwoElecsTwoCores(atom, epc);
                    i=1;
                    for mu = 1:dxy
                        for nu = mu:dxy
                            j=1;
                            for lambda = 1:dxy
                                for sigma = lambda:dxy
                                    obj.twoElecsAtomEpcCoresMpiBuff(a,b,i,j) ...
                                        = diatomicTwoElecsTwoCores(mu,nu,lambda,sigma);
                                    j = j + 1;
                                end
                            end
                            i = i + 1;
                        end
                    end
                end
            end
            
            for a = 1:totalNumberAtoms
                for b = 1:totalNumberEpcs
                    i=1;
                    for mu = 1:dxy
                        for nu = mu:dxy
                            j=1;
                            for lambda = 1:dxy
                                for sigma = lambda:dxy
                                    value = obj.twoElecsAtomEpcCoresMpiBuff(a,b,i,j);
                                    obj.twoElecsAtomEpcCores(a,b,mu,nu,lambda,sigma) = value;
                                    obj.twoElecsAtomEpcCores(a,b,mu,nu,sigma,lambda) = value;
                                    obj.twoElecsAtomEpcCores(a,b,nu,mu,lambda,sigma) = value;
                                    obj.twoElecsAtomEpcCores(a,b,nu,mu,sigma,lambda) = value;
                                    j = j + 1;
                                end
                            end
                            i = i + 1;
                        end
                    end
                end
            end
        end
        
        function value = GetAuxiliaryDiatomCoreRepulsionEnergy(obj, atomA, atomB, distanceAB)
            alphaA = obj.AtomGetNddoAlpha(atomA);
            alphaB = obj.AtomGetNddoAlpha(atomB);
            ang2AU = SEQC.Arguments.GetInstance().GetAngstrom2AU();
            if(atomA.atomType == 1 && (atomB.atomType == 7 || atomB.atomType == 8) ) % H N O
                value = 1.0 + (distanceAB/ang2AU)*exp(-alphaB*distanceAB) + exp(-alphaA*distanceAB);
            elseif(atomB.atomType == 1 && (atomA.atomType == 7 || atomA.atomType == 8) ) % H N O
                value = 1.0 + (distanceAB/ang2AU)*exp(-alphaA*distanceAB) + exp(-alphaB*distanceAB);
            else
                value = 1.0 + exp(-alphaA*distanceAB) + exp(-alphaB*distanceAB);
            end
        end
        
        %    double GetAuxiliaryDiatomCoreRepulsionEnergy1stDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                              const MolDS_base_atoms::Atom& atomB,
        %                                                              double distanceAB,
        %                                                              MolDS_base::CartesianType axisA) const;
        %    double GetAuxiliaryDiatomCoreRepulsionEnergy2ndDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                              const MolDS_base_atoms::Atom& atomB,
        %                                                              double distanceAB,
        %                                                              MolDS_base::CartesianType axisA1,
        %                                                              MolDS_base::CartesianType axisA2) const;
        %    double GetCISCoefficientMOEnergy(int k, int l, int r, int numberActiveVir) const;
        %    double GetCISCoefficientTwoElecIntegral(int k, int l, int p, int q, int r, int s, int numberActiveVir) const;
        %    void CalcHessianSCF(double** hessianSCF, bool isMassWeighted) const;
        %    double GetHessianElementSameAtomsSCF(int indexAtomA,
        %                                         MolDS_base::CartesianType axisA1,
        %                                         MolDS_base::CartesianType axisA2,
        %                                         double const* const*               orbitalElectronPopulation,
        %                                         double const* const* const* const* orbitalElectronPopulation1stDerivs,
        %                                         double const* const* const* const*                      diatomicOverlapAOs1stDerivs,
        %                                         double const* const* const* const* const*               diatomicOverlapAOs2ndDerivs,
        %                                         double const* const* const* const* const* const*        diatomicTwoElecsTwoCores1stDerivs,
        %                                         double const* const* const* const* const* const* const* diatomicTwoElecsTwoCores2ndDerivs) const;
        %    double GetHessianElementDifferentAtomsSCF(int indexAtomA,
        %                                              int indexAtomB,
        %                                              MolDS_base::CartesianType axisA,
        %                                              MolDS_base::CartesianType axisB,
        %                                              double const* const*               orbitalElectronPopulation,
        %                                              double const* const* const* const* orbitalElectronPopulation1stDerivs,
        %                                              double const* const* const* const*                      diatomicOverlapAOs1stDerivs,
        %                                              double const* const* const* const* const*               diatomicOverlapAOs2ndDerivs,
        %                                              double const* const* const* const* const* const*        diatomicTwoElecsTwoCores1stDerivs,
        %                                              double const* const* const* const* const* const* const* diatomicTwoElecsTwoCores2ndDerivs) const;
        %    void MallocTempMatricesEachThreadCalcHessianSCF(double*****    diatomicOverlapAOs1stDerivs,
        %                                                    double******   diatomicOverlapAOs2ndDerivs,
        %                                                    double*******  diatomicTwoElecsTwoCores1stDerivs,
        %                                                    double******** diatomicTwoElecsTwoCores2ndDerivs,
        %                                                    double***      tmpRotMat,
        %                                                    double***      tmpRotMat1stDeriv,
        %                                                    double****     tmpRotMat1stDerivs,
        %                                                    double*****    tmpRotMat2ndDerivs,
        %                                                    double*****    tmpDiatomicTwoElecsTwoCores,
        %                                                    double******   tmpDiatomicTwoElecsTwoCores1stDerivs,
        %                                                    double***      tmpDiaOverlapAOsInDiaFrame,
        %                                                    double***      tmpDiaOverlapAOs1stDerivInDiaFrame,
        %                                                    double***      tmpDiaOverlapAOs2ndDerivInDiaFrame,
        %                                                    double****     tmpDiaOverlapAOs1stDerivs,
        %                                                    double*****    tmpDiaOverlapAOs2ndDerivs,
        %                                                    double***      tmpRotatedDiatomicOverlap,
        %                                                    double**       tmpRotatedDiatomicOverlapVec,
        %                                                    double***      tmpMatrixBC,
        %                                                    double**      tmpVectorBC) const;
        %    void FreeTempMatricesEachThreadCalcHessianSCF(double*****    diatomicOverlapAOs1stDerivs,
        %                                                  double******   diatomicOverlapAOs2ndDerivs,
        %                                                  double*******  diatomicTwoElecsTwoCores1stDerivs,
        %                                                  double******** diatomicTwoElecsTwoCores2ndDerivs,
        %                                                  double***      tmpRotMat,
        %                                                  double***      tmpRotMat1stDeriv,
        %                                                  double****     tmpRotMat1stDerivs,
        %                                                  double*****    tmpRotMat2ndDerivs,
        %                                                  double*****    tmpDiatomicTwoElecsTwoCores,
        %                                                  double******   tmpDiatomicTwoElecsTwoCores1stDerivs,
        %                                                  double***      tmpDiaOverlapAOsInDiaFrame,
        %                                                  double***      tmpDiaOverlapAOs1stDerivInDiaFrame,
        %                                                  double***      tmpDiaOverlapAOs2ndDerivInDiaFrame,
        %                                                  double****     tmpDiaOverlapAOs1stDerivs,
        %                                                  double*****    tmpDiaOverlapAOs2ndDerivs,
        %                                                  double***      tmpRotatedDiatomicOverlap,
        %                                                  double**       tmpRotatedDiatomicOverlapVec,
        %                                                  double***      tmpMatrixBC,
        %                                                  double**       tmpVectorBC) const;
        %    double GetAuxiliaryHessianElement1(int mu,
        %                                       int nu,
        %                                       int indexAtomA,
        %                                       int indexAtomC,
        %                                       MolDS_base::CartesianType axisA1,
        %                                       MolDS_base::CartesianType axisA2,
        %                                       double const* const* orbitalElectronPopulation,
        %                                       double const* const* const* const* const* const* diatomicTwoElecsTwoCores2ndDerivs) const;
        %    double GetAuxiliaryHessianElement2(int mu,
        %                                       int nu,
        %                                       int indexAtomA,
        %                                       int indexAtomB,
        %                                       int indexAtomC,
        %                                       MolDS_base::CartesianType axisA,
        %                                       MolDS_base::CartesianType axisB,
        %                                       double const* const* const* const* orbitalElectronPopulation1stDerivs,
        %                                       double const* const* const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
        %    double GetAuxiliaryHessianElement3(int lambda,
        %                                       int sigma,
        %                                       int indexAtomA,
        %                                       int indexAtomC,
        %                                       MolDS_base::CartesianType axisA1,
        %                                       MolDS_base::CartesianType axisA2,
        %                                       double const* const* orbitalElectronPopulation,
        %                                       double const* const* const* const* const* const* diatomicTwoElecsTwoCores2ndDerivs) const;
        %    double GetAuxiliaryHessianElement4(int lambda,
        %                                       int sigma,
        %                                       int indexAtomA,
        %                                       int indexAtomB,
        %                                       int indexAtomC,
        %                                       MolDS_base::CartesianType axisA,
        %                                       MolDS_base::CartesianType axisB,
        %                                       double const* const* const* const* orbitalElectronPopulation1stDerivs,
        %                                       double const* const* const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
        %    double GetAuxiliaryHessianElement5(int mu,
        %                                       int lambda,
        %                                       int indexAtomA,
        %                                       int indexAtomC,
        %                                       MolDS_base::CartesianType axisA1,
        %                                       MolDS_base::CartesianType axisA2,
        %                                       double const* const* orbitalElectronPopulation,
        %                                       double const* const* const* const* diatomicOverlapAOs2ndDerivs) const;
        %    double GetAuxiliaryHessianElement6(int mu,
        %                                       int lambda,
        %                                       int indexAtomA,
        %                                       int indexAtomB,
        %                                       int indexAtomC,
        %                                       MolDS_base::CartesianType axisA,
        %                                       MolDS_base::CartesianType axisB,
        %                                       double const* const* const* const* orbitalElectronPopulation1stDerivs,
        %                                       double const* const* const* diatomicOverlapAOs1stDerivs) const;
        %    double GetAuxiliaryHessianElement7(int mu,
        %                                       int nu,
        %                                       int lambda,
        %                                       int sigma,
        %                                       int indexAtomA,
        %                                       int indexAtomC,
        %                                       MolDS_base::CartesianType axisA1,
        %                                       MolDS_base::CartesianType axisA2,
        %                                       double const* const* orbitalElectronPopulation,
        %                                       double const* const* const* const* const* const* diatomicTwoElecsTwoCores2ndDerivs) const;
        %    double GetAuxiliaryHessianElement8(int mu,
        %                                       int nu,
        %                                       int lambda,
        %                                       int sigma,
        %                                       int indexAtomA,
        %                                       int indexAtomB,
        %                                       int indexAtomC,
        %                                       MolDS_base::CartesianType axisA,
        %                                       MolDS_base::CartesianType axisB,
        %                                       double const* const* orbitalElectronPopulation,
        %                                       double const* const* const* const* orbitalElectronPopulation1stDerivs,
        %                                       double const* const* const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
        %    void CalcOrbitalElectronPopulation1stDerivatives(double**** orbitalElectronPopulation1stDerivatives) const;
        %    void SolveCPHF(double** solutionsCPHF,
        %                   const std::vector<MoIndexPair>& nonRedundantQIndeces,
        %                   const std::vector<MoIndexPair>& redundantQIndeces) const;
        %    void CalcStaticFirstOrderFocks(double** staticFirstOrderFocks,
        %                                   const std::vector<MoIndexPair>& nonRedundantQIndeces,
        %                                   const std::vector<MoIndexPair>& redundantQIndeces) const;
        %    void CalcStaticFirstOrderFock(double* staticFirstOrderFock,
        %                                  const std::vector<MoIndexPair>& nonRedundantQIndeces,
        %                                  const std::vector<MoIndexPair>& redundantQIndeces,
        %                                  int indexAtomA,
        %                                  MolDS_base::CartesianType axisA) const;
        %    void MallocTempMatricesStaticFirstOrderFock(double****** diatomicTwoElecsTwoCores1stDeriv,
        %                                                double****   diatomicOverlapAOs1stDeriv,
        %                                                double***    tmpRotMat,
        %                                                double****   tmpRotMat1stDerivs,
        %                                                double*****  tmpDiatomicTwoElecTwo) const;
        %    void FreeTempMatricesStaticFirstOrderFock(double****** diatomicTwoElecsTwoCores1stDeriv,
        %                                              double****   diatomicOverlapAOs1stDeriv,
        %                                              double***    tmpRotMat,
        %                                              double****   tmpRotMat1stDerivs,
        %                                              double*****  tmpDiatomicTwoElecTwo) const;
        %    void CalcMatrixCPHF(double** matrixCPHF,
        %                        const std::vector<MoIndexPair>& nonRedundantQIndeces,
        %                        const std::vector<MoIndexPair>& redundantQIndeces) const;
        %    void MallocTempMatricesSolveCPHF(double*** matrixCPHF,
        %                                     int dimensionCPHF) const;
        %    void FreeTempMatricesSolveCPHF(double*** matrixCPHF,
        %                                   int dimensionCPHF) const;
        %    void CalcHeatsFormation(double* heatsFormation,
        %                            const MolDS_base::Molecule& molecule) const;
        
        function res = GetElectronCoreAttraction(obj, indexAtomA, indexAtomB, mu, nu)
            atomB = obj.molecule.atomVect{indexAtomB};
            res = -atomB.coreCharge*reshape(obj.twoElecsTwoAtomCores(indexAtomA,indexAtomB,mu,nu,1,1), [], 1);
        end
        
        %    double GetElectronCoreAttraction1stDerivative(int indexAtomA,
        %                                                  int indexAtomB,
        %                                                  int mu,
        %                                                  int nu,
        %                                                  double const* const* const* const* const* diatomicTwoElecsTwoCores1stDerivatives,
        %                                                  MolDS_base::CartesianType axisA) const;
        
        function matrix = CalcDiatomicTwoElecsTwoCores(obj, atomA, atomB)
            dxy = 4;
            if(atomA.atomType ~= 255 && atomB.atomType ~= 255) % uint8(EnumAtom.EPC) = 255
                if(atomA.index == atomB.index)
                    throw(MException('Mndo:CalcDiatomicTwoElecsTwoCores', 'Atom indices cannot be equal.'));
                end
            end
            if(atomA.atomType == 255 && atomB.atomType == 255)
                if(atomA.index == atomB.index)
                    throw(MException('Mndo:CalcDiatomicTwoElecsTwoCores', 'EPC atom indices cannot be equal.'));
                end
            end
            
            matrix = zeros(dxy, dxy, dxy, dxy);
            
            % calclation in diatomic frame
            for mu = 1:atomA.GetValenceSize()
                for nu = mu:atomA.GetValenceSize()
                    for lambda = 1:atomB.GetValenceSize()
                        for sigma = lambda:atomB.GetValenceSize()
%                             value = obj.GetNddoRepulsionIntegral(...
%                                 atomA, ...
%                                 atomA.valence(mu),...
%                                 atomA.valence(nu),...
%                                 atomB, ...
%                                 atomB.valence(lambda),...
%                                 atomB.valence(sigma));...
%                                 matrix(mu,nu,lambda,sigma) = value;
                            value = obj.Cpp_GetNddoRepulsionIntegral(...
                                atomA, ...
                                atomA.valence(mu),...
                                atomA.valence(nu),...
                                atomB, ...
                                atomB.valence(lambda),...
                                atomB.valence(sigma));...
                                matrix(mu,nu,lambda,sigma) = value;
                            matrix(mu,nu,sigma,lambda) = value;
                            matrix(nu,mu,lambda,sigma) = value;
                            matrix(nu,mu,sigma,lambda) = value;
                        end
                    end
                end
            end
            % rotate matirix into the space frame
            tmpRotMat = obj.CalcRotatingMatrix(atomA, atomB);
            tmpRotMat = tmpRotMat(1:dxy, 1:dxy);
            lined_tmpRotMat = reshape(tmpRotMat',dxy*dxy,1);
            expand_tmpRotMat = reshape(permute(reshape(lined_tmpRotMat*lined_tmpRotMat',dxy,dxy,dxy,dxy),[1 3 2 4]),dxy*dxy,dxy*dxy);
            matrix = reshape(expand_tmpRotMat'*reshape(matrix,16,16)*expand_tmpRotMat,4,4,4,4);
            
        end
        
        %    void CalcDiatomicTwoElecsTwoCores1stDerivatives(double***** matrix,
        %                                                    double**    tmpRotMat,
        %                                                    double***   tmpRotMat1stDerivs,
        %                                                    double****  tmpDiatomicTwoElecsTwoCores,
        %                                                    int indexAtomA,
        %                                                    int indexAtomB) const;
        %    void CalcDiatomicTwoElecsTwoCores2ndDerivatives(double****** matrix,
        %                                                    double**     tmpRotMat,
        %                                                    double***    tmpRotMat1stDerivs,
        %                                                    double****   tmpRotMat2ndDerivs,
        %                                                    double****   tmpDiatomicTwoElecsTwoCores,
        %                                                    double*****  tmpDiatomicTwoElecsTwoCores1stDerivs,
        %                                                    int indexAtomA,
        %                                                    int indexAtomB) const;
        
        
        
        %    void RotateDiatomicTwoElecsTwoCores1stDerivativesToSpaceFrame(double***** matrix,
        %                                                                  double const* const* const* const* diatomicTwoElecsTwoCores,
        %                                                                  double const* const* rotatingMatrix,
        %                                                                  double const* const* const* rotMat1stDerivatives) const;
        %    void RotateDiatomicTwoElecsTwoCores2ndDerivativesToSpaceFrame(double****** matrix,
        %                                                                  double const* const* const* const*        diatomicTwoElecsTwoCores,
        %                                                                  double const* const* const* const* const* diatomicTwoElecsTwoCores1stDerivatives,
        %                                                                  double const* const* rotatingMatrix,
        %                                                                  double const* const* const* rotMat1stDerivatives,
        %                                                                  double const* const* const* const* rotMat2ndDerivatives) const;
        %    void MallocTempMatricesRotateDiatomicTwoElecsTwoCores1stDerivs(double*** twiceRotatingMatrix,
        %                                                                   double*** twiceRotatingMatrixDerivA,
        %                                                                   double*** twiceRotatingMatrixDerivB,
        %                                                                   double*** oldMatrix,
        %                                                                   double*** rotatedMatrix,
        %                                                                   double**  tmpRotatedVec,
        %                                                                   double*** tmpMatrix,
        %                                                                   double**  tmpVector,
        %                                                                   double*** ptrDiatomic) const;
        %    void FreeTempMatricesRotateDiatomicTwoElecsTwoCores1stDerivs(double*** twiceRotatingMatrix,
        %                                                                 double*** twiceRotatingMatrixDerivA,
        %                                                                 double*** twiceRotatingMatrixDerivB,
        %                                                                 double*** oldMatrix,
        %                                                                 double*** rotatedMatrix,
        %                                                                 double**  tmpRotatedVec,
        %                                                                 double*** tmpMatrix,
        %                                                                 double**  tmpVector,
        %                                                                 double*** ptrDiatomic) const;
        
        value = GetNddoRepulsionIntegral(obj, atomA, mu, nu, atomB, lambda, sigma);
        
        value = Cpp_GetNddoRepulsionIntegral(obj, atomA, mu, nu, atomB, lambda, sigma);
        
        %    double GetNddoRepulsionIntegral1stDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                 MolDS_base::OrbitalType mu,
        %                                                 MolDS_base::OrbitalType nu,
        %                                                 const MolDS_base_atoms::Atom& atomB,
        %                                                 MolDS_base::OrbitalType lambda,
        %                                                 MolDS_base::OrbitalType sigma,
        %                                                 MolDS_base::CartesianType axisA) const;
        %    double GetNddoRepulsionIntegral2ndDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                 MolDS_base::OrbitalType mu,
        %                                                 MolDS_base::OrbitalType nu,
        %                                                 const MolDS_base_atoms::Atom& atomB,
        %                                                 MolDS_base::OrbitalType lambda,
        %                                                 MolDS_base::OrbitalType sigma,
        %                                                 MolDS_base::CartesianType axisA1,
        %                                                 MolDS_base::CartesianType axisA2) const;
        
        value = GetSemiEmpiricalMultipoleInteraction(obj, atomA, atomB, multipoleA, multipoleB, rAB);
        
        %    double GetSemiEmpiricalMultipoleInteraction1stDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                             const MolDS_base_atoms::Atom& atomB,
        %                                                             MolDS_base::MultipoleType multipoleA,
        %                                                             MolDS_base::MultipoleType multipoleB,
        %                                                             double rAB) const;
        %    double GetSemiEmpiricalMultipoleInteraction2ndDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                             const MolDS_base_atoms::Atom& atomB,
        %                                                             MolDS_base::MultipoleType multipoleA,
        %                                                             MolDS_base::MultipoleType multipoleB,
        %                                                             double rAB) const;
        %    void MallocTempMatricesCalcForce(double****   diatomicOverlapAOs1stDerivs,
        %                                     double****** diatomicTwoElecsTwoCores1stDerivs,
        %                                     double***    tmpDiaOverlapAOsInDiaFrame,
        %                                     double***    tmpDiaOverlapAOs1stDerivInDiaFrame,
        %                                     double***    tmpRotMat,
        %                                     double***    tmpRotMat1stDeriv,
        %                                     double****   tmpRotMat1stDerivs,
        %                                     double***    tmpRotatedDiatomicOverlap,
        %                                     double**     tmpRotatedDiatomicOverlapVec,
        %                                     double***    tmpMatrixBC,
        %                                     double**     tmpVectorBC,
        %                                     double*****  tmpDiatomicTwoElecsTwoCores) const;
        %    void FreeTempMatricesCalcForce(double****   diatomicOverlapAOs1stDerivs,
        %                                   double****** diatomicTwoElecsTwoCores1stDerivs,
        %                                   double***    tmpDiaOverlapAOsInDiaFrame,
        %                                   double***    tmpDiaOverlapAOs1stDerivInDiaFrame,
        %                                   double***    tmpRotMat,
        %                                   double***    tmpRotMat1stDeriv,
        %                                   double****   tmpRotMat1stDerivs,
        %                                   double***    tmpRotatedDiatomicOverlap,
        %                                   double**     tmpRotatedDiatomicOverlapVec,
        %                                   double***    tmpMatrixBC,
        %                                   double**     tmpVectorBC,
        %                                   double*****  tmpDiatomicTwoElecsTwoCores) const;
        %    void CalcForceSCFElecCoreAttractionPart(double* force,
        %                                            int indexAtomA,
        %                                            int indexAtomB,
        %                                            double const* const* const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
        %    void CalcForceSCFOverlapAOsPart(double* force,
        %                                    int indexAtomA,
        %                                    int indexAtomB,
        %                                    double const* const* const* diatomicOverlapAOs1stDerivs) const;
        %    void CalcForceSCFTwoElecPart(double* force,
        %                                 int indexAtomA,
        %                                 int indexAtomB,
        %                                 double const* const* const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
        %    void CalcForceExcitedStaticPart(double* force,
        %                                    int elecStateIndex,
        %                                    int indexAtomA,
        %                                    int indexAtomB,
        %                                    double const* const* const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
        %    void CalcForceExcitedElecCoreAttractionPart(double* force,
        %                                                int elecStateIndex,
        %                                                int indexAtomA,
        %                                                int indexAtomB,
        %                                                double const* const* const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
        %    void CalcForceExcitedTwoElecPart(double* force,
        %                                     int elecStateIndex,
        %                                     int indexAtomA,
        %                                     int indexAtomB,
        %                                     double const* const* const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
        
    end
    
    methods (Static, Access = protected)
        
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, lambda, sigma, rAB);
        
    end
    
end