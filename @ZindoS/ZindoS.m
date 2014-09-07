classdef ZindoS < Cndo2
    
    properties (SetAccess = protected)
        
        zMatrixForce;
        etaMatrixForce;
        MoIndexPair; % MoIndexPair{int moI; int moJ; bool isMoICIMO; bool isMoJCIMO;};
        
    end
    
    properties (SetAccess = private)
        
        nishimotoMatagaMatrix;
        matrixForceElecStatesNum;
        nishimotoMatagaParamA;
        nishimotoMatagaParamB;
        overlapAOsCorrectionSigma;
        overlapAOsCorrectionPi;
        zMatrixForceElecStatesNum;
        etaMatrixForceElecStatesNum;
        
    end
    
    methods (Access = public)
        
        function obj = ZindoS()
            % protected variables and methods
            obj.theory = EnumTheory.ZINDOS;
            obj.SetEnableAtomTypes();
            
            obj.zMatrixForceElecStatesNum = 0;
            obj.etaMatrixForceElecStatesNum = 0;
            obj.zMatrixForce = [];
            obj.etaMatrixForce = [];
            
            % private variables
            obj.nishimotoMatagaMatrix = [];
            obj.matrixForceElecStatesNum = 0;
            obj.nishimotoMatagaParamA = 1.2;
            obj.nishimotoMatagaParamB = 2.4;
            obj.overlapAOsCorrectionSigma = 1.267;
            obj.overlapAOsCorrectionPi = 0.585;
        end
        
        function SetMolecule(obj, mol)
            SetMolecule@Cndo2(obj, mol);
            obj.nishimotoMatagaMatrix = zeros(length(obj.molecule.atomVect), 9, ...
                length(obj.molecule.atomVect), 9);
        end
        
        %        void DoCIS();
        %        void OutputCISResults() const;
        %        void CalcOverlapSingletSDsWithAnotherElectronicStructure(double** overlapSingletSDs,
        %        double const* const* overlapMOs) const;
        %        void CalcOverlapESsWithAnotherElectronicStructure(double** overlapESs,
        %        double const* const* overlapSingletSDs,
        %        const MolDS_base::ElectronicStructure& lhsElectronicStructure) const;
        
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
        end
        
        %    virtual void CalcCISProperties();
        %    virtual void CalcElectronicTransitionDipoleMoment(double* transitionDipoleMoment,
        %                                                      int to, int from,
        %                                                      double const* const* fockMatrix,
        %                                                      double const* const* matrixCIS,
        %                                                      double const* const* const* cartesianMatrix,
        %                                                      const MolDS_base::Molecule& molecule,
        %                                                      double const* const* orbitalElectronPopulation,
        %                                                      double const* const* overlapAOs,
        %                                                      double const* groundStateDipole) const;
        
        function gammaAB = CalcGammaAB(~)
            gammaAB = [];
        end
        
        function value = GetFockDiagElement(obj, atomA, mu, isGuess)
            indexAtomA = atomA.index;
            firstAOIndexA = atomA.firstAOIndex;
            value = obj.AtomGetCoreIntegral(atomA, atomA.valence(mu-firstAOIndexA+1));
            if(~isGuess)
                
                orbitalElectronPopulationDiagPart = diag(obj.orbitalElectronPopulation);
                
                orbitalMu = atomA.valence(mu-firstAOIndexA+1);
                atomANumberValence = atomA.GetValenceSize();
                temp = 0.0;
                for v = 1:atomANumberValence
                    orbitalLam = atomA.valence(v);
                    coulomb  = obj.GetCoulombInt(orbitalMu, orbitalLam, atomA);
                    exchange = obj.GetExchangeInt(orbitalMu, orbitalLam, atomA);
                    lammda = v + firstAOIndexA - 1;
                    temp = temp + orbitalElectronPopulationDiagPart(lammda)*(coulomb - 0.5*exchange);
                end
                value = value + temp;
                
                temp = 0.0;
                totalNumberAtoms = length(obj.molecule.atomVect);
                for B = 1:totalNumberAtoms
                    if(B ~= indexAtomA)
                        atomB = obj.molecule.atomVect{B};
                        atomBNumberValence = atomB.GetValenceSize();
                        for i = 1:atomBNumberValence
                            sigma = i + atomB.firstAOIndex - 1;
                            orbitalSigma = atomB.valence(i);
                            temp = temp + orbitalElectronPopulationDiagPart(sigma)...
                                *obj.nishimotoMatagaMatrix(indexAtomA,orbitalMu,B,orbitalSigma);
                        end
                        temp = temp - atomB.coreCharge ...
                            *obj.nishimotoMatagaMatrix(indexAtomA,1,B,1);
                    end
                end
                value = value + temp;
            end
        end
        
        function value = GetFockOffDiagElement(obj, atomA, atomB, mu, nu, isGuess)
            indexAtomA = atomA.index;
            indexAtomB = atomB.index;
            orbitalMu = atomA.valence(mu-atomA.firstAOIndex+1);
            orbitalNu = atomB.valence(nu-atomB.firstAOIndex+1);
            bondParameter = 0.5*(obj.AtomGetBondingParameter(atomA, orbitalMu) ...
                +obj.AtomGetBondingParameter(atomB, orbitalNu));
            if(isGuess)
%                 disp([mu,nu,bondParameter]);
                value = bondParameter*obj.overlapAOs(mu,nu);
            else
                if(indexAtomA == indexAtomB)
                    coulomb  = obj.GetCoulombInt(orbitalMu, orbitalNu, atomA);
                    exchange = obj.GetExchangeInt(orbitalMu, orbitalNu, atomA);
                    value = (1.5*exchange - 0.5*coulomb)*obj.orbitalElectronPopulation(mu,nu);
                else
                    value = bondParameter*obj.overlapAOs(mu,nu);
                    value = value - 0.5*obj.orbitalElectronPopulation(mu,nu)...
                        *obj.nishimotoMatagaMatrix(indexAtomA,orbitalMu,indexAtomB,orbitalNu);
                end
            end
        end
        
        function diatomicOverlapAOs = CalcDiatomicOverlapAOsInDiatomicFrame(obj, atomA, atomB)
            diatomicOverlapAOs = CalcDiatomicOverlapAOsInDiatomicFrame@Cndo2(obj, atomA, atomB);
            
            % see (4f) in [AEZ_1986]
            diatomicOverlapAOs(3,3) = diatomicOverlapAOs(3,3) * obj.overlapAOsCorrectionSigma;
            diatomicOverlapAOs(2,2) = diatomicOverlapAOs(2,2) * obj.overlapAOsCorrectionPi;
            diatomicOverlapAOs(4,4) = diatomicOverlapAOs(4,4) * obj.overlapAOsCorrectionPi;
            
        end
        
        %    virtual void CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(double** diatomicOverlapAOsDeri,
        %                                                                    const MolDS_base_atoms::Atom& atomA,
        %                                                                    const MolDS_base_atoms::Atom& atomB) const;
        %    virtual void CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(double** diatomicOverlapAOs2ndDeri,
        %                                                                    const MolDS_base_atoms::Atom& atomA,
        %                                                                    const MolDS_base_atoms::Atom& atomB) const;
        
        function value = GetCoulombInt(~, orbital1, orbital2, atom)
            if( orbital1 == 1 && orbital2 == 1)
                value = atom.paramPool.GetZindoF0ssLower();
            elseif( orbital1 == 1 && ( orbital2 == 4 || orbital2 == 2 || orbital2 == 3 ))
                value = atom.paramPool.GetZindoF0ssLower();
            elseif( orbital2 == 1 && ( orbital1 == 4 || orbital1 == 2 || orbital1 == 3 ))
                value = atom.paramPool.GetZindoF0ssLower();
            elseif( (orbital1 == orbital2) && ( orbital1 == 4 || orbital1 == 2 || orbital1 == 3 ))
                value = atom.paramPool.GetZindoF0ssLower()...
                    +atom.paramPool.GetZindoF2ppLower()*4.0;
            elseif( (orbital1 ~= orbital2) ...
                    && ( orbital1 == 4 || orbital1 == 2 || orbital1 == 3 )...
                    && ( orbital2 == 4 || orbital2 == 2 || orbital2 == 3 ) )
                value = atom.paramPool.GetZindoF0ssLower()...
                    -atom.paramPool.GetZindoF2ppLower()*2.0;
            else
                throw(MException('ZindoS:GetCoulombInt', 'Orbital type wrong.'));
            end
        end
        
        function value = GetExchangeInt(obj, orbital1, orbital2, atom)
            if( orbital1 == orbital2)
                value = obj.GetCoulombInt(orbital1, orbital2, atom);
            elseif( orbital1 == 1 && (orbital2 == 4 || orbital2 == 2 || orbital2 == 3 ) )
                value = atom.paramPool.GetZindoG1spLower();
            elseif( orbital2 == 1 && (orbital1 == 4 || orbital1 == 2 || orbital1 == 3 ) )
                value = atom.paramPool.GetZindoG1spLower();
            elseif( (orbital1 ~= orbital2)...
                    && ( orbital1 == 4 || orbital1 == 2 || orbital1 == 3 )...
                    && ( orbital2 == 4 || orbital2 == 2 || orbital2 == 3 ) )
                value = atom.paramPool.GetZindoF2ppLower()*3.0;
            else
                throw(MException('ZindoS:GetExchangeInt', 'Orbital type wrong.'));
            end
        end
        
        function CalcTwoElecsTwoCores(obj)
            obj.CalcNishimotoMatagaMatrix();
        end
        
        %    virtual double GetMolecularIntegralElement(int moI,
        %                                               int moJ,
        %                                               int moK,
        %                                               int moL,
        %                                               const MolDS_base::Molecule& molecule,
        %                                               double const* const* fockMatrix,
        %                                               double const* const* gammaAB) const;
        %    virtual void CalcCISMatrix(double** matrixCIS) const;
        %    double GetCISDiagElement(double const* energiesMO,
        %                             double const* const* const* const* nishimotoMatagaMatrix,
        %                             const MolDS_base::Molecule& molecule,
        %                             double const* const* fockMatrix,
        %                             int moI,
        %                             int moA) const;
        %    double GetCISOffDiagElement(double const* const* const* const* nishimotoMatagaMatrix,
        %                                const MolDS_base::Molecule& molecule,
        %                                double const* const* fockMatrix,
        %                                int moI,
        %                                int moA,
        %                                int moJ,
        %                                int moB) const;
        %    bool RequiresExcitedStatesForce(const std::vector<int>& elecStates) const;
        %    virtual void CalcForce(const std::vector<int>& elecStates);
        %    int GetSlaterDeterminantIndex(int activeOccIndex, int activeVirIndex) const;
        %    int GetActiveOccIndex(const MolDS_base::Molecule& molecule, int matrixCISIndex) const;
        %    int GetActiveVirIndex(const MolDS_base::Molecule& molecule, int matrixCISIndex) const;
        %    void CheckMatrixForce(const std::vector<int>& elecStates);
        %    void CalcEtaMatrixForce(const std::vector<int>& elecStates);
        %    void CalcZMatrixForce(const std::vector<int>& elecStates);
        %    void CalcActiveSetVariablesQ(std::vector<MoIndexPair>* nonRedundantQIndeces,
        %                                 std::vector<MoIndexPair>* redundantQIndeces,
        %                                 int numberActiveOcc,
        %                                 int numberActiveVir) const;
        %    virtual double GetSmallQElement(int moI,
        %                                    int moP,
        %                                    double const* const* xiOcc,
        %                                    double const* const* xiVir,
        %                                    double const* const* eta) const;
        %    double GetGammaNRElement(int moI, int moJ, int moK, int moL) const;
        %    double GetGammaRElement (int moI, int moJ, int moK, int moL) const;
        %    double GetNNRElement    (int moI, int moJ, int moK, int moL) const;
        %    double GetNRElement     (int moI, int moJ, int moK, int moL) const;
        %    double GetKNRElement    (int moI, int moJ, int moK, int moL) const;
        %    double GetKRElement     (int moI, int moJ, int moK, int moL) const;
        %    virtual double GetAuxiliaryKNRKRElement(int moI, int moJ, int moK, int moL) const;
        %    void CalcForceExcitedOverlapAOsPart(double* force,
        %                                        int elecStateIndex,
        %                                        int indexAtomA,
        %                                        int indexAtomB,
        %                                        double const* const* const* diatomicOverlapAOs1stDerivs) const;
        
        function value = AtomGetCoreIntegral(~, atom, orbital)
            if(orbital == 1)
                value = -1.0*atom.paramPool.zindoIonPotS...
                    - atom.paramPool.GetZindoJss()*(atom.paramPool.zindoL-1)...
                    - atom.paramPool.GetZindoJsp()*(atom.paramPool.zindoM)...
                    - atom.paramPool.GetZindoJsd()*(atom.paramPool.zindoN);
            elseif(orbital == 4 || orbital == 2 || orbital == 3)
                value = -1.0*atom.paramPool.zindoIonPotP...
                    - atom.paramPool.GetZindoJpp()*(atom.paramPool.zindoM-1)...
                    - atom.paramPool.GetZindoJsp()*(atom.paramPool.zindoL)...
                    - atom.paramPool.GetZindoJpd()*(atom.paramPool.zindoN);
            elseif(orbital == 5 || orbital == 6 || orbital == 7 || orbital == 8 || orbital == 9 )
                value = -1.0*atom.paramPool.zindoIonPotD...
                    - atom.paramPool.GetZindoJdd()*(atom.paramPool.zindoN-1)...
                    - atom.paramPool.GetZindoJsd()*(atom.paramPool.zindoL)...
                    - atom.paramPool.GetZindoJpd()*(atom.paramPool.zindoM);
            else
                throw(MException('ZindoS:AtomGetCoreIntegral', 'Orbital type wrong.'));
            end
        end
        
        function value = AtomGetBondingParameter(~, atom, orbital)
            if(orbital == 1 || orbital == 4 || orbital == 2 || orbital == 3 )
                value = atom.paramPool.zindoBondingParameterS;
            elseif( orbital == 5 ||...
                    orbital == 6 ||...
                    orbital == 7 ||...
                    orbital == 8 ||...
                    orbital == 9 )
                value = atom.paramPool.zindoBondingParameterD;
            else
                throw(MException('ZindoS:AtomGetBondingParameter', 'Orbital type wrong.'));
            end
        end
        
    end
    
    methods (Access = private)
        
        %    void DoCISDirect();
        %    void DoCISDavidson();
        %    void OutputCISDipole() const;
        %    void OutputCISTransitionDipole() const;
        %    void OutputCISMulliken() const;
        %    void OutputCISUnpairedPop() const;
        %    void CalcFreeExcitonEnergies(double** freeExcitonEnergiesCIS,
        %                                 const MolDS_base::Molecule& molecule,
        %                                 double const* energiesMO,
        %                                 double const* const* matrixCIS,
        %                                 int matrixCISdimension) const;
        %    void CalcOrbitalElectronPopulationCIS(double**** orbitalElectronPopulationCIS,
        %                                          double const* const* orbitalElectronPopulation,
        %                                          const MolDS_base::Molecule& molecule,
        %                                          double const* const* fockMatrix,
        %                                          double const* const* matrixCIS) const;
        %    void CalcAtomicElectronPopulationCIS(double*** atomicElectronPopulationCIS,
        %                                         double const* const* const* orbitalElectronPopulationCIS,
        %                                         const MolDS_base::Molecule& molecule) const;
        %    void CalcAtomicUnpairedPopulationCIS(double*** atomicUnpairedPopulationCIS,
        %                                         double const* const* const* orbitalElectronPopulationCIS,
        %                                         const MolDS_base::Molecule& molecule) const;
        %    void CalcElectronicDipoleMomentsExcitedStates(double*** electronicTransitionDipoleMoments,
        %                                                  double const* const* fockMatrix,
        %                                                  double const* const* matrixCIS,
        %                                                  double const* const* const* cartesianMatrix,
        %                                                  const MolDS_base::Molecule& molecule,
        %                                                  double const* const* orbitalElectronPopulation,
        %                                                  double const* const* overlapAOs) const;
        %    void CalcElectronicTransitionDipoleMoments(double*** electronicTransitionDipoleMoments,
        %                                               double const* const* fockMatrix,
        %                                               double const* const* matrixCIS,
        %                                               double const* const* const* cartesianMatrix,
        %                                               const MolDS_base::Molecule& molecule,
        %                                               double const* const* orbitalElectronPopulation,
        %                                               double const* const* overlapAOs) const;
        
        function res = GetNishimotoMatagaTwoEleInt(obj, atomA, orbitalA, atomB, orbitalB, rAB)
            
            if(orbitalA == 1 || ...
                    orbitalA == 4 ||...
                    orbitalA == 2 ||...
                    orbitalA == 3 )
                gammaAA = atomA.paramPool.zindoF0ss;
            else
                throw(MException('ZindoS:GetNishimotoMatagaTwoEleInt', 'Orbital type wrong.'));
            end
            
            if(orbitalB == 1 || ...
                    orbitalB == 4 ||...
                    orbitalB == 2 ||...
                    orbitalB == 3 )
                gammaBB = atomB.paramPool.zindoF0ss;
            else
                throw(MException('ZindoS:GetNishimotoMatagaTwoEleInt', 'Orbital type wrong.'));
            end
            
            gamma=gammaAA+gammaBB;
            res = obj.nishimotoMatagaParamA/( rAB+obj.nishimotoMatagaParamB/gamma );
            
        end
        
        %    double GetNishimotoMatagaTwoEleInt1stDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                    MolDS_base::OrbitalType orbitalA,
        %                                                    const MolDS_base_atoms::Atom& atomB,
        %                                                    MolDS_base::OrbitalType orbitalB,
        %                                                    MolDS_base::CartesianType axisA) const;// ref. [MN_1957] and (5a) in [AEZ_1986]
        %    double GetNishimotoMatagaTwoEleInt1stDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                    MolDS_base::OrbitalType orbitalA,
        %                                                    const MolDS_base_atoms::Atom& atomB,
        %                                                    MolDS_base::OrbitalType orbitalB,
        %                                                    const double rAB,
        %                                                    MolDS_base::CartesianType axisA) const;// ref. [MN_1957] and (5a) in [AEZ_1986]
        
        function CalcNishimotoMatagaMatrix(obj)
            totalNumberAtoms = length(obj.molecule.atomVect);
            for A = 1:totalNumberAtoms
                atomA = obj.molecule.atomVect{A};
                firstAOIndexA = atomA.firstAOIndex;
                lastAOIndexA  = atomA.GetLastAOIndex();
                for mu = firstAOIndexA:lastAOIndexA
                    orbitalMu = atomA.valence(mu-firstAOIndexA+1);
                    for B = A:totalNumberAtoms
                        atomB = obj.molecule.atomVect{B};
                        firstAOIndexB = atomB.firstAOIndex;
                        lastAOIndexB  = atomB.GetLastAOIndex();
                        rAB = obj.molecule.GetDistanceAtoms(atomA, atomB);
                        for nu = firstAOIndexB:lastAOIndexB
                            orbitalNu = atomB.valence(nu-firstAOIndexB+1);
                            obj.nishimotoMatagaMatrix(A,orbitalMu,B,orbitalNu) = obj.GetNishimotoMatagaTwoEleInt(atomA, ...
                                orbitalMu, ...
                                atomB, ...
                                orbitalNu,...
                                rAB);
                            if(A~=B)
                                obj.nishimotoMatagaMatrix(B,orbitalNu,A,orbitalMu) = obj.nishimotoMatagaMatrix(A,orbitalMu,B,orbitalNu);
                            end
                        end
                    end
                end
            end
        end
        
        %    void CalcRitzVector(double* ritzVector,
        %                        double const* const* expansionVectors,
        %                        double const* const* interactionMatrix,
        %                        int interactionMatrixDimension,
        %                        int ritzVectorIndex) const;
        %    void CalcResidualVectorAndNorm(double* residualVector,
        %                                   double* norm,
        %                                   double const* ritzVector,
        %                                   double const* interactionEigenEnergies,
        %                                   int residualVectorIndex) const;
        %    void SortCISEigenVectorCoefficients(std::vector<CISEigenVectorCoefficient>* cisEigenVectorCoefficients,
        %                                        double* cisEigenVector) const;
        %    void SortSingleExcitationSlaterDeterminants(std::vector<MoEnergyGap>* moEnergyGaps) const;
        %    void UpdateExpansionVectors(double** expansionVectors,
        %                                int* notConvergedStates,
        %                                double const* interactionEigenEnergies,
        %                                double const* residualVector,
        %                                int interactionMatrixDimension,
        %                                int residualVectorIndex) const;
        %    void CalcInteractionMatrix(double** interactionMatrix,
        %                               double const* const* expansionVectors,
        %                               int interactionMatrixDimension) const;
        %    void FreeDavidsonCISTemporaryMtrices(double*** expansionVectors,
        %                                         double** residualVector,
        %                                         double** ritzVector) const;
        %    void FreeDavidsonRoopCISTemporaryMtrices(double*** interactionMatrix,
        %                                             int interactionMatrixDimension,
        %                                             double** interactionEigenEnergies) const;
        %    void CalcDiatomicTwoElecsTwoCores1stDerivatives(double*** matrix,
        %                                                    int indexAtomA,
        %                                                    int indexAtomB) const;
        %    void MallocTempMatricesCalcForce(double**** diatomicOverlapAOs1stDerivs,
        %                                     double**** diatomicTwoElecsTwoCores1stDerivs,
        %                                     double***  tmpDiaOverlapAOsInDiaFrame,
        %                                     double***  tmpDiaOverlapAOs1stDerivInDiaFrame,
        %                                     double***  tmpRotMat,
        %                                     double***  tmpRotMat1stDeriv,
        %                                     double**** tmpRotMat1stDerivs,
        %                                     double***  tmpRotatedDiatomicOverlap,
        %                                     double**   tmpRotatedDiatomicOverlapVec,
        %                                     double***  tmpMatrixBC,
        %                                     double**   tmpVectorBC) const;
        %    void FreeTempMatricesCalcForce(double**** diatomicOverlapAOs1stDerivs,
        %                                   double**** diatomicTwoElecsTwoCores1stDerivs,
        %                                   double***  tmpDiaOverlapAOsInDiaFrame,
        %                                   double***  tmpDiaOverlapAOs1stDerivInDiaFrame,
        %                                   double***  tmpRotMat,
        %                                   double***  tmpRotMat1stDeriv,
        %                                   double**** tmpRotMat1stDerivs,
        %                                   double***  tmpRotatedDiatomicOverlap,
        %                                   double**   tmpRotatedDiatomicOverlapVec,
        %                                   double***  tmpMatrixBC,
        %                                   double**   tmpVectorBC) const;
        %    void CalcForceExcitedStaticPart(double* force,
        %                                    int elecStateIndex,
        %                                    int indexAtomA,
        %                                    int indexAtomB,
        %                                    double const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
        %    void CalcForceExcitedElecCoreAttractionPart(double* force,
        %                                                int elecStateIndex,
        %                                                int indexAtomA,
        %                                                int indexAtomB,
        %                                                double const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
        %    void CalcForceExcitedTwoElecPart(double* force,
        %                                     int elecStateIndex,
        %                                     int indexAtomA,
        %                                     int indexAtomB,
        %                                     double const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
        %    void CheckZMatrixForce(const std::vector<int>& elecStates);
        %    void CheckEtaMatrixForce(const std::vector<int>& elecStates);
        %    double GetZMatrixForceElement(double const* y,
        %                                  double const* q,
        %                                  double const* const* transposedFockMatrix,
        %                                  const std::vector<MoIndexPair>& nonRedundantQIndeces,
        %                                  const std::vector<MoIndexPair>& redundantQIndeces,
        %                                  int mu,
        %                                  int nu) const;
        %    void MallocTempMatrixForZMatrix(double** delta,
        %                                    double** q,
        %                                    double*** gammaNRMinusKNR,
        %                                    double*** kRDag,
        %                                    double** y,
        %                                    double*** transposedFockMatrix,
        %                                    double*** xiOcc,
        %                                    double*** xiVir,
        %                                    int sizeQNR,
        %                                    int sizeQR) const;
        %    void FreeTempMatrixForZMatrix(double** delta,
        %                                  double** q,
        %                                  double*** gammaNRMinusKNR,
        %                                  double*** kRDag,
        %                                  double** y,
        %                                  double*** transposedFockMatrix,
        %                                  double*** xiOcc,
        %                                  double*** xiVir,
        %                                  int sizeQNR,
        %                                  int sizeQR) const;
        %    void CalcDeltaVector(double* delta, int exciteState) const;
        %    void CalcQVector(double* q,
        %                     double const* delta,
        %                     double const* const* xiOcc,
        %                     double const* const* xiVir,
        %                     double const* const* eta,
        %                     const std::vector<MoIndexPair>& nonRedundantQIndeces,
        %                     const std::vector<MoIndexPair>& redundantQIndeces) const;
        %    void CalcXiMatrices(double** xiOcc,
        %                        double** xiVir,
        %                        int exciteState,
        %                        double const* const* transposedFockMatrix) const;
        %    void CalcAuxiliaryVector(double* y,
        %                             double const* q,
        %                             double const* const* kRDagerGammaRInv,
        %                             const std::vector<MoIndexPair>& nonRedundantQIndeces,
        %                             const std::vector<MoIndexPair>& redundantQIndeces) const;
        %    double GetKRDagerElement(int moI, int moJ, int moK, int moL) const;
        %    void CalcGammaNRMinusKNRMatrix(double** gammaNRMinusKNR,
        %                                   const std::vector<MoIndexPair>& nonRedundantQIndeces) const;
        %    void CalcKRDagerGammaRInvMatrix(double** kRDagerGammaRInv,
        %                                    const std::vector<MoIndexPair>& nonRedundantQIndeces,
        %                                    const std::vector<MoIndexPair>& redundantQIndeces) const;
        
    end
    
end




