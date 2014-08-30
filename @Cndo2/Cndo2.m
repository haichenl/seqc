classdef Cndo2 < handle
    
    properties (SetAccess = protected)
        
        enableAtomTypes = {};
        molecule;
        theory;

        coreRepulsionEnergy;
        coreEpcCoulombEnergy;
        vdWCorrectionEnergy;
        matrixCISdimension;
        fockMatrix;
        energiesMO;
        orbitalElectronPopulation; %P_{\mu\nu} of (2.50) in J. A. Pople book.
        orbitalElectronPopulationCIS;
        atomicElectronPopulation; %P_{AB} of (3.21) in J. A. Pople book.
        atomicElectronPopulationCIS;
        atomicUnpairedPopulationCIS;
        overlapAOs; % overlap integral between AOs
        twoElecsTwoAtomCores;
        twoElecsAtomEpcCores;
        cartesianMatrix; % cartesian matrix represented by AOs
        
        % Diagnonal terms are electronic dipole moments of each eigenstates
        % (i.e. electronicDipole[0,0,XAxis] is the x-component of the electronic
        % dipole moment of the ground state. electronicDipole[10,10,XAxis] is
        % the x-component of the electronic dipole moment of the 10-th excited state).
        % Off-diagonal terms are transition dipole moments between eigenstates
        % (i.e. electronicDipole[10,0,XAxis] is the x-component of the transition
        % dipole moment from the ground state to 10-th excited state.).
        electronicTransitionDipoleMoments;
        
        coreDipoleMoment; % dipole moment of configuration.
        normalForceConstants; % force constants of normal modes
        normalModes; % in mass-weighted coordinates
        matrixCIS;
        excitedEnergies;
        freeExcitonEnergiesCIS;
        matrixForce;

    end
    
    properties (SetAccess = private)
        
        elecSCFEnergy;
        bondingAdjustParameterK = zeros(1,2); %see (3.79) in J. A. Pople book
        gammaAB;
        enableAtomTypesVdW = {};
        
        ReducedOverlapAOsParameters_Z;
        ReducedOverlapAOsParameters_Y;
   
    end
    
    methods (Access = public)
        
        function obj = Cndo2()
            %protected variables
            obj.molecule = [];
            obj.theory = TheoryType.CNDO2;
            obj.coreRepulsionEnergy  = 0.0;
            obj.coreEpcCoulombEnergy = 0.0;
            obj.vdWCorrectionEnergy  = 0.0;
            obj.matrixCISdimension   = 0;
            obj.fockMatrix = [];
            obj.energiesMO = [];
            obj.orbitalElectronPopulation    = [];
            obj.orbitalElectronPopulationCIS = [];
            obj.atomicElectronPopulation     = [];
            obj.atomicElectronPopulationCIS  = [];
            obj.atomicUnpairedPopulationCIS  = [];
            obj.overlapAOs                   = [];
            obj.twoElecsTwoAtomCores         = [];
            obj.twoElecsAtomEpcCores         = [];
            obj.cartesianMatrix              = [];
            obj.electronicTransitionDipoleMoments = [];
            obj.coreDipoleMoment             = [];
            obj.normalForceConstants         = [];
            obj.normalModes                  = [];
            obj.matrixCIS                    = [];
            obj.excitedEnergies              = [];
            obj.freeExcitonEnergiesCIS       = [];
            obj.matrixForce                  = [];
            
            %protected methods
            obj.SetEnableAtomTypes();
            obj.SetEnableAtomTypesVdW();
            
            %private variables
            obj.elecSCFEnergy = 0.0;
            obj.bondingAdjustParameterK(1) = 1.000; %see (3.79) in J. A. Pople book
            obj.bondingAdjustParameterK(2) = 0.750; %see (3.79) in J. A. Pople book
            obj.gammaAB = [];
            
            % ReducedOverlapAOsParameters nested class stuff
            zcontainer = load('cndo2_z.mat', 'Z');
            ycontainer = load('cndo2_y.mat', 'Y');
            obj.ReducedOverlapAOsParameters_Z = zcontainer.Z;
            obj.ReducedOverlapAOsParameters_Y = ycontainer.Y;
            
        end
        
        function SetMolecule(obj, mol)
            obj.molecule = mol;
            obj.CheckNumberValenceElectrons(mol);
            obj.CheckEnableAtomType(mol);
            obj.CheckEnableAtomTypeVdW(mol);
            nbf = obj.molecule.GetTotalNumberAOs();
            obj.fockMatrix = zeros(nbf);
            obj.energiesMO = zeros(1, nbf);
            obj.orbitalElectronPopulation = zeros(nbf);
            obj.atomicElectronPopulation = zeros(1, length(obj.molecule.GetAtomVect()));
            obj.overlapAOs = zeros(nbf);
            obj.cartesianMatrix = zeros(3, nbf, nbf);
            electronicTransitionDipoleMomentsDim = 1;
%             if(Parameters::GetInstance()->RequiresCIS()){
%                 electronicTransitionDipoleMomentsDim += Parameters::GetInstance()->GetNumberExcitedStatesCIS();
%                 }
            obj.electronicTransitionDipoleMoments = zeros(electronicTransitionDipoleMomentsDim);
            obj.coreDipoleMoment = zeros(1, 3);
            if(obj.theory == TheoryType.CNDO2 || obj.theory == TheoryType.INDO)
                obj.gammaAB = zeros(length(obj.molecule.GetAtomVect()));
            end
        end
        
        function DoSCF(obj)
            requiresGuess = true;
            
            % Cndo2::MallocSCFTemporaryMatrices
            diisNumErrorVect = Parameters.GetInstance().GetDiisNumErrorVectSCF();
            oldOrbitalElectronPopulation = zeros(obj.molecule.GetTotalNumberAOs());
            if(0<diisNumErrorVect)
                diisStoredDensityMatrix = zeros(diisNumErrorVect, obj.molecule.GetTotalNumberAOs(), obj.molecule.GetTotalNumberAOs());
                diisStoredErrorVect = zeros(diisNumErrorVect, obj.molecule.GetTotalNumberAOs(), obj.molecule.GetTotalNumberAOs());
                diisErrorProducts = zeros(diisNumErrorVect+1);
                tmpDiisErrorProducts = zeros(diisNumErrorVect+1);
                diisErrorCoefficients = zeros(1, diisNumErrorVect+1);
            end
            
            % calculate electron integral
            obj.gammaAB = obj.CalcGammaAB(obj.molecule);
            obj.overlapAOs = obj.CalcOverlapAOs(obj.molecule);
            obj.cartesianMatrix = obj.CalcCartesianMatrixByGTOExpansion(obj.molecule, STOnGType.STO6G);
            [obj.twoElecsTwoAtomCores, obj.twoElecsAtomEpcCores] = obj.CalcTwoElecsTwoCores(obj.molecule);
            
            % SCF
            maxIterationsSCF = Parameters.GetInstance().GetMaxIterationsSCF();
            hasAppliedDIIS=false;
            hasAppliedDamping=false;
            diisError=0.0;
            for iterationStep = 1:maxIterationsSCF
                obj.atomicElectronPopulation = obj.CalcAtomicElectronPopulation(...
                    obj.orbitalElectronPopulation, ...
                    obj.molecule);
                oldOrbitalElectronPopulation = obj.orbitalElectronPopulation;
                
                isGuess = (iterationStep==1 && requiresGuess);
                % continue here
                obj.fockMatrix = obj.CalcFockMatrix(obj.molecule, ...
                    obj.overlapAOs, ...
                    obj.gammaAB,...
                    obj.orbitalElectronPopulation, ...
                    obj.atomicElectronPopulation,...
                    obj.twoElecsTwoAtomCores,...
                    isGuess);
                
                % diagonalization of the Fock matrix
                [orbital, eigenvalue] = eig(obj.fockMatrix);
                [obj.energiesMO, order] = sort(diag(eigenvalue));
                orbital = orbital(:, order);
                numberTotalValenceElectrons = obj.molecule.GetTotalNumberValenceElectrons();
                obj.orbitalElectronPopulation = 2 .* orbital(:, 1:numberTotalValenceElectrons/2) ...
                    *orbital(:, 1:numberTotalValenceElectrons/2)';
                
                % check convergence
                hasConverged = obj.SatisfyConvergenceCriterion(oldOrbitalElectronPopulation, ...
                    obj.orbitalElectronPopulation,...
                    obj.molecule.GetTotalNumberAOs());
                if(hasConverged)
                    obj.CalcSCFProperties();
                    break;
                else
                    if(~isGuess)
                        %                obj.DoDIIS(obj.orbitalElectronPopulation,
                        %                             oldOrbitalElectronPopulation,
                        %                             diisStoredDensityMatrix,
                        %                             diisStoredErrorVect,
                        %                             diisErrorProducts,
                        %                             tmpDiisErrorProducts,
                        %                             diisErrorCoefficients,
                        %                             diisError,
                        %                             hasAppliedDIIS,
                        %                             Parameters::GetInstance()->GetDiisNumErrorVectSCF(),
                        %                             *obj.molecule,
                        %                             iterationStep);
                        %                obj.DoDamp(rmsDensity,
                        %                             hasAppliedDamping,
                        %                             obj.orbitalElectronPopulation,
                        %                             oldOrbitalElectronPopulation,
                        %                             *obj.molecule);
                    end
                end
                
                % SCF fails
                if(iterationStep==maxIterationsSCF-1)
                    throw(MException('Cndo2:DoSCF', 'SCF not converged.'));
                end
            end
        end
        
        %         GetFockMatrix()
        %         GetEnergiesMO()
        %
        %         % todo: cis, force
        %
        %         GetElectronicEnergy(int elecState)
        %         GetCoreRepulsionEnergy()
        %         GetVdWCorrectionEnergy()
        %         CalcOverlapAOsWithAnotherConfiguration(double** overlapAOs,const MolDS_base::Molecule& lhsMolecule)
        %         CalcOverlapMOsWithAnotherElectronicStructure(double** overlapMOs,
        %                                                      double const* const* overlapAOs,
        %                                                      const MolDS_base::ElectronicStructure& lhsElectronicStructure)
        %         CalcOverlapSingletSDsWithAnotherElectronicStructure(double** overlapSingletSDs,
        %                                                                     double const* const* overlapMOs) const;
        %         CalcOverlapESsWithAnotherElectronicStructure(double** overlapESs,
        %                                                              double const* const* overlapSingletSDs,
        %                                                              const MolDS_base::ElectronicStructure& lhsElectronicStructure)
        
        function res = GetTheoryType(obj)
            res = obj.theory;
        end
        
    end
    
    methods (Access = protected)
        
        function SetEnableAtomTypes(obj)
            obj.enableAtomTypes = {};
            obj.enableAtomTypes{end+1} = AtomType.H;
            obj.enableAtomTypes{end+1} = AtomType.Li;
            %obj.enableAtomTypes{end+1} = AtomType.Be;
            %obj.enableAtomTypes{end+1} = AtomType.B;
            obj.enableAtomTypes{end+1} = AtomType.C;
            obj.enableAtomTypes{end+1} = AtomType.N;
            obj.enableAtomTypes{end+1} = AtomType.O;
            obj.enableAtomTypes{end+1} = AtomType.F;
            %obj.enableAtomTypes{end+1} = AtomType.Na;
            %obj.enableAtomTypes{end+1} = AtomType.Mg;
            %obj.enableAtomTypes{end+1} = AtomType.Al;
            %obj.enableAtomTypes{end+1} = AtomType.Si;
            %obj.enableAtomTypes{end+1} = AtomType.P;
            obj.enableAtomTypes{end+1} = AtomType.S;
            obj.enableAtomTypes{end+1} = AtomType.Cl;
        end
        
        function CalcSCFProperties(obj)
            obj.atomicElectronPopulation = obj.CalcAtomicElectronPopulation( ...
                obj.orbitalElectronPopulation, ...
                obj.molecule);
            obj.CalcCoreRepulsionEnergy();
            if(Parameters.GetInstance().RequiresVdWSCF())
                obj.CalcVdWCorrectionEnergy();
            end
            obj.elecSCFEnergy = obj.CalcElecSCFEnergy(obj.molecule, ...
                obj.energiesMO, ...
                obj.fockMatrix, ...
                obj.gammaAB, ...
                obj.coreRepulsionEnergy, ...
                obj.coreEpcCoulombEnergy, ...
                obj.vdWCorrectionEnergy);
            obj.coreDipoleMoment = obj.CalcCoreDipoleMoment(obj.molecule);
            obj.electronicTransitionDipoleMoments = obj.CalcElectronicDipoleMomentGroundState( ...
                obj.cartesianMatrix, ...
                obj.molecule, ...
                obj.orbitalElectronPopulation, ...
                obj.overlapAOs);
            %    groundState = 0;
            %    if(Parameters::GetInstance()->RequiresFrequencies() &&
            %       Parameters::GetInstance()->GetElectronicStateIndexFrequencies() == groundState){
            %       this->CalcNormalModes(this->normalModes, this->normalForceConstants, *this->molecule);
            %    }
        end
        
        %         CalcNormalModes(double** normalModes, double* normalForceConstants, const MolDS_base::Molecule& molecule) const;
        
        function transitionDipoleMoment = CalcElectronicTransitionDipoleMoment(~, ...
                to, from, ~, ~, cartesianMatrix, molecule, ...
                orbitalElectronPopulation, overlapAOs, ~)
            groundState = 1;
            if(from == groundState && to == groundState)
                dipoleCenter = molecule.GetXyzDipoleCenter();
                transitionDipoleMoment = zeros(1, 3);
                transitionDipoleMoment(1) = transitionDipoleMoment(1) - sum(sum(orbitalElectronPopulation.*reshape(cartesianMatrix(1,:,:), size(cartesianMatrix(1,:,:),2), [])));
                transitionDipoleMoment(2) = transitionDipoleMoment(2) - sum(sum(orbitalElectronPopulation.*reshape(cartesianMatrix(2,:,:), size(cartesianMatrix(2,:,:),2), [])));
                transitionDipoleMoment(3) = transitionDipoleMoment(3) - sum(sum(orbitalElectronPopulation.*reshape(cartesianMatrix(3,:,:), size(cartesianMatrix(3,:,:),2), [])));
                
                % set orign of dipole
                temp = sum(sum(orbitalElectronPopulation.*overlapAOs));
                transitionDipoleMoment = transitionDipoleMoment + dipoleCenter.*temp;
            else
                throw(MException('Cndo2:CalcElectronicTransitionDipoleMoment', 'Excited transition moment not enabled for this theory.'));
            end
        end
        
        function value = GetBondingAdjustParameterK(obj, shellA, shellB)
             value=1.0;
             if(shellA >= ShellType.mShell || shellB >= ShellType.mShell)
                 value = obj.bondingAdjustParameterK(2); % notice c++->matlab index conversion
             end
        end
        
        %         GetAtomCoreEpcCoulombEnergy (const MolDS_base_atoms::Atom& atom,
        %                                                const MolDS_base_atoms::Atom& epc) const;
        
        function res = GetDiatomCoreRepulsionEnergy(obj, atomA, atomB)
            distance = obj.molecule.GetDistanceAtoms(atomA, atomB);
            res = atomA.GetCoreCharge()*atomB.GetCoreCharge()/distance;
        end
        
        %    virtual double GetDiatomCoreRepulsion1stDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                       const MolDS_base_atoms::Atom& atomB,
        %                                                       MolDS_base::CartesianType axisA) const;
        %    virtual double GetDiatomCoreRepulsion2ndDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                       const MolDS_base_atoms::Atom& atomB,
        %                                                       MolDS_base::CartesianType axisA1,
        %                                                       MolDS_base::CartesianType axisA2) const;
        
        function res = GetDiatomVdWCorrectionEnergy(obj, atomA, atomB)
            distance = obj.molecule.GetDistanceAtoms(atomA, atomB);
            vdWDistance = atomA.GetVdWRadii() + atomB.GetVdWRadii();
            tmpSum = atomA.GetVdWCoefficient()+atomB.GetVdWCoefficient();
            if(tmpSum<=0)
                res = 0;
                return;
            end
            vdWCoefficients = 2.0*atomA.GetVdWCoefficient()*atomB.GetVdWCoefficient()/tmpSum;
            damping = obj.GetVdWDampingValue(vdWDistance, distance);
            scalingFactor = Parameters.GetInstance().GetVdWScalingFactorSCF();
            res =  -1.0*scalingFactor*vdWCoefficients*damping ...
                /(distance*distance*distance*distance*distance*distance);
        end
        
        %    virtual double GetDiatomVdWCorrection1stDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                       const MolDS_base_atoms::Atom& atomB,
        %                                                       MolDS_base::CartesianType axisA) const;
        %    virtual double GetDiatomVdWCorrection2ndDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                       const MolDS_base_atoms::Atom& atomB,
        %                                                       MolDS_base::CartesianType axisA1,
        %                                                       MolDS_base::CartesianType axisA2) const;
        
        function value = GetReducedOverlapAOs(obj, arg1, arg2, arg3, arg4, arg5, arg6, arg7)
            if(nargin == 8)
                na = arg1;
                la = arg2;
                m = arg3;
                nb = arg4;
                lb = arg5;
                alpha = arg6;
                beta = arg7;
                value = 0.0;
                I = 2*double(ShellType.ShellType_end)-1;
                J = 2*double(ShellType.ShellType_end)-1;
                for i = 0:I-1
                    for j = 0:J-1
                        temp = obj.ReducedOverlapAOsParameters_Y(na+1,nb+1,la+1,lb+1,m+1,i+1,j+1);
                        if(0<abs(temp))
                            temp = temp * obj.GetAuxiliaryA(i, 0.5*(alpha+beta));
                            temp = temp * obj.GetAuxiliaryB(j, 0.5*(alpha-beta));
                            value = value + temp;
                        end
                    end
                end
                value = value * obj.GetAuxiliaryD(la, lb, m);
            elseif(nargin == 5)
                na = arg1;
                nb = arg2;
                alpha = arg3;
                beta = arg4;
                value = 0.0;
                for k = 0:na+nb
                    temp = obj.ReducedOverlapAOsParameters_Z(na+1,nb+1,k+1);
                    if(0<abs(temp))
                        temp = temp * obj.GetAuxiliaryA(k, 0.5*(alpha+beta));
                        temp = temp * obj.GetAuxiliaryB(na+nb-k, 0.5*(alpha-beta));
                        value = value + temp;
                    end
                end
                value = value * 0.5;
            else
                throw(MException('Cndo2:GetReducedOverlapAOs', 'Input argument number wrong.'));
            end
        end
        
        %    double GetReducedOverlapAOs1stDerivativeAlpha    (int na, int la, int m, int nb, int lb, double alpha, double beta) const;
        %    double GetReducedOverlapAOs1stDerivativeBeta     (int na, int la, int m, int nb, int lb, double alpha, double beta) const;
        %    double GetReducedOverlapAOs2ndDerivativeAlpha    (int na, int la, int m, int nb, int lb, double alpha, double beta) const;
        %    double GetReducedOverlapAOs2ndDerivativeBeta     (int na, int la, int m, int nb, int lb, double alpha, double beta) const;
        %    double GetReducedOverlapAOs2ndDerivativeAlphaBeta(int na, int la, int m, int nb, int lb, double alpha, double beta) const;
        %    double GetOverlapAOsElement1stDerivativeByGTOExpansion(const MolDS_base_atoms::Atom& atomA,
        %                                                           int valenceIndexA,
        %                                                           const MolDS_base_atoms::Atom& atomB,
        %                                                           int valenceIndexB,
        %                                                           MolDS_base::STOnGType stonG,
        %                                                           MolDS_base::CartesianType axisA) const; // See [DY_1977].
        
        % see J. Mol. Struc. (Theochem), 419, 19 (1997) (ref. [BFB_1997])
        % we set gamma=0 always.
        function rotatingMatrix = CalcRotatingMatrix(~, atomA, atomB)
            rotatingMatrix = zeros(double(OrbitalType.OrbitalType_end));
            eulerAngle = EulerAngle(atomB.GetXyz() - atomA.GetXyz());
            alpha = eulerAngle.GetAlpha();
            beta  = eulerAngle.GetBeta();
            
            s = double(OrbitalType.s);
            py = double(OrbitalType.py);
            pz = double(OrbitalType.pz);
            px = double(OrbitalType.px);
            dxy = double(OrbitalType.dxy);
            dyz = double(OrbitalType.dyz);
            dzz = double(OrbitalType.dzz);
            dzx = double(OrbitalType.dzx);
            dxxyy = double(OrbitalType.dxxyy);
            % rotating matrix for s-function
            rotatingMatrix(s,s) = 1.0;
            
            % rotating matrix for p-function
            % dMatrix is (53) with gamma=0 in J. Mol. Strct. 419, 19(1997) (ref. [BFB_1997])
            rotatingMatrix(py,py) = cos(alpha);
            rotatingMatrix(py,pz) = sin(alpha)*sin(beta);
            rotatingMatrix(py,px) = sin(alpha)*cos(beta);
            
            rotatingMatrix(pz,py) = 0.0;
            rotatingMatrix(pz,pz) = cos(beta);
            rotatingMatrix(pz,px) = -1.0*sin(beta);
            
            rotatingMatrix(px,py) = -1.0*sin(alpha);
            rotatingMatrix(px,pz) = cos(alpha)*sin(beta);
            rotatingMatrix(px,px) = cos(alpha)*cos(beta);
            
            % rotating matrix for d-function
            % dMatrix is (37) in J. Mol. Strct. 419, 19(1997) (ref. [BFB_1997])
            dMatrix = zeros(double(OrbitalType.OrbitalType_end));
            dMatrix(dzz,dzz) = 0.5*(3.0*(cos(beta)*cos(beta)) - 1.0);
            dMatrix(dxxyy,dxxyy) = cos(0.5*beta)*cos(0.5*beta)*cos(0.5*beta)*cos(0.5*beta);
            dMatrix(dzx,dzx) = (2.0*cos(beta)-1.0)*cos(0.5*beta)*cos(0.5*beta);
            dMatrix(dxxyy,dzx) = -2.0*sin(0.5*beta)*cos(0.5*beta)*cos(0.5*beta)*cos(0.5*beta);
            dMatrix(dxxyy,dzz) = sqrt(6.0)*(sin(0.5*beta)*sin(0.5*beta))*((cos(0.5*beta))*(cos(0.5*beta)));
            dMatrix(dxxyy,dyz) = -2.0*sin(0.5*beta)*sin(0.5*beta)*sin(0.5*beta)*cos(0.5*beta);
            dMatrix(dxxyy,dxy) = sin(0.5*beta)*sin(0.5*beta)*sin(0.5*beta)*sin(0.5*beta);
            dMatrix(dzx,dzz) = -sqrt(6.0)*cos(beta)*cos(0.5*beta)*sin(0.5*beta);
            dMatrix(dzx,dyz) = (2.0*cos(beta)+1.0)*(sin(0.5*beta)*sin(0.5*beta));
            
            rotatingMatrix(dxy,dxy) = cos(2.0*alpha)*            (dMatrix(dxxyy,dxxyy) - dMatrix(dxxyy,dxy));
            rotatingMatrix(dxy,dyz) = cos(2.0*alpha)*            (-1.0*dMatrix(dxxyy,dzx) - dMatrix(dxxyy,dyz));
            rotatingMatrix(dxy,dzz) = sqrt(2.0)*sin(2.0*alpha)*  dMatrix(dxxyy,dzz);
            rotatingMatrix(dxy,dzx) = sin(2.0*alpha)*            (-1.0*dMatrix(dxxyy,dzx) + dMatrix(dxxyy,dyz));
            rotatingMatrix(dxy,dxxyy) = sin(2.0*alpha)*          (dMatrix(dxxyy,dxxyy) + dMatrix(dxxyy,dxy));
            
            rotatingMatrix(dyz,dxy) = cos(alpha)*                (dMatrix(dxxyy,dzx) + dMatrix(dxxyy,dyz));
            rotatingMatrix(dyz,dyz) = cos(alpha)*                (dMatrix(dzx,dzx) + dMatrix(dzx,dyz));
            rotatingMatrix(dyz,dzz) = -1.0*sqrt(2.0)*sin(alpha)* dMatrix(dzx,dzz);
            rotatingMatrix(dyz,dzx) = sin(alpha)*                (dMatrix(dzx,dzx) - dMatrix(dzx,dyz));
            rotatingMatrix(dyz,dxxyy) = sin(alpha)*              (dMatrix(dxxyy,dzx) - dMatrix(dxxyy,dyz));
            
            rotatingMatrix(dzz,dxy) = 0.0;
            rotatingMatrix(dzz,dyz) = 0.0;
            rotatingMatrix(dzz,dzz) = dMatrix(dzz,dzz);
            rotatingMatrix(dzz,dzx) = sqrt(2.0)*dMatrix(dzx,dzz);
            rotatingMatrix(dzz,dxxyy) = sqrt(2.0)*dMatrix(dxxyy,dzz);
            
            rotatingMatrix(dzx,dxy) = -1.0*sin(alpha)*           (dMatrix(dxxyy,dzx) + dMatrix(dxxyy,dyz));
            rotatingMatrix(dzx,dyz) = -1.0*sin(alpha)*           (dMatrix(dzx,dzx) + dMatrix(dzx,dyz));
            rotatingMatrix(dzx,dzz) = -1.0*sqrt(2.0)*cos(alpha)* dMatrix(dzx,dzz);
            rotatingMatrix(dzx,dzx) = cos(alpha)*                (dMatrix(dzx,dzx) - dMatrix(dzx,dyz));
            rotatingMatrix(dzx,dxxyy) = cos(alpha)*              (dMatrix(dxxyy,dzx) - dMatrix(dxxyy,dyz));
            
            rotatingMatrix(dxxyy,dxy) = -1.0*sin(2.0*alpha)*     (dMatrix(dxxyy,dxxyy) - dMatrix(dxxyy,dxy));
            rotatingMatrix(dxxyy,dyz) = -1.0*sin(2.0*alpha)*     (-1.0*dMatrix(dxxyy,dzx) - dMatrix(dxxyy,dyz));
            rotatingMatrix(dxxyy,dzz) = sqrt(2.0)*cos(2.0*alpha)*dMatrix(dxxyy,dzz);
            rotatingMatrix(dxxyy,dzx) = cos(2.0*alpha)*          (-1.0*dMatrix(dxxyy,dzx) + dMatrix(dxxyy,dyz));
            rotatingMatrix(dxxyy,dxxyy) = cos(2.0*alpha)*        (dMatrix(dxxyy,dxxyy) + dMatrix(dxxyy,dxy));
        end
        
        function gammaAB = CalcGammaAB(obj, molecule)
            totalAtomNumber = length(molecule.GetAtomVect());
            gammaAB = zeros(totalAtomNumber);
            atomvect = molecule.GetAtomVect();
            for A = 1:totalAtomNumber
                atomA = atomvect{A};
                na = double(atomA.GetValenceShellType());
                orbitalExponentA = atomA.GetOrbitalExponent(atomA.GetValenceShellType(), OrbitalType.s, obj.theory);
                for B = A:totalAtomNumber
                    atomB = atomvect{B};
                    nb = double(atomB.GetValenceShellType());
                    orbitalExponentB = atomB.GetOrbitalExponent(atomB.GetValenceShellType(), OrbitalType.s, obj.theory);
                    
                    R = molecule.GetDistanceAtoms(A, B);
                    if(R>0.0)
                        % (B.56)
                        value = power(0.5*R, 2.0*na);
                        value = value * obj.GetReducedOverlapAOs(2*na-1, 0, 2.0*orbitalExponentA*R, 0);
                        
                        for l = 1:2*nb
                            temp = l;
                            temp = temp * power(2.0*orbitalExponentB, 2*nb-l);
                            temp = temp / (factorial(2*nb-l)*2.0*nb);
                            temp = temp * power(0.5*R, 2.0*nb-l+2.0*na);
                            temp = temp * obj.GetReducedOverlapAOs(2*na-1, 2*nb-l, 2.0*orbitalExponentA*R, 2.0*orbitalExponentB*R);
                            value = value - temp;
                        end
                        
                        value = value * power(2.0*orbitalExponentA, 2.0*na+1.0);
                        value = value / factorial(2*na);
                    else
                        % (B.62)
                        value =  factorial(2*na-1);
                        value = value / power(2.0*orbitalExponentA, 2.0*na);
                        
                        for l = 1:2*nb
                            temp = l;
                            temp = temp * power(2.0*orbitalExponentB, 2*nb-l);
                            temp = temp * factorial(2*na+2*nb-l-1);
                            temp = temp / factorial(2*nb-l);
                            temp = temp / (2.0*nb);
                            temp = temp / power( 2.0*orbitalExponentA + 2.0*orbitalExponentB, 2.0*(na+nb)-l );
                            value = value - temp;
                        end
                        value = value * power(2.0*orbitalExponentA, 2.0*na+1);
                        value = value / factorial(2*na);
                    end
                    gammaAB(A, B) = value;
                end
            end
            gammaAB = gammaAB + gammaAB' - diag(diag(gammaAB));
        end
        
        function value = GetFockDiagElement(obj, atomA, indexAtomA, mu, molecule, ...
                gammaAB, orbitalElectronPopulation, ...
                atomicElectronPopulation, ~, isGuess)
            firstAOIndexA = atomA.GetFirstAOIndex();
            value = atomA.GetCoreIntegral(atomA.GetValence(mu-firstAOIndexA+1), ...
                gammaAB(indexAtomA, indexAtomA), ...
                isGuess, obj.theory);
            if(~isGuess)
                atomvect = molecule.GetAtomVect();
                temp = atomicElectronPopulation(indexAtomA) ...
                    -0.5*orbitalElectronPopulation(mu, mu);
                value = value + temp*gammaAB(indexAtomA, indexAtomA);
                
                temp = 0.0;
                for BB = 1:length(molecule.GetAtomVect())
                    if(BB ~= indexAtomA)
                        atomBB = atomvect{BB};
                        temp = temp + ( atomicElectronPopulation(BB) - atomBB.GetCoreCharge()  )...
                            *gammaAB(indexAtomA, BB);
                    end
                end
                value = value + temp;
            end
        end
        
        function value = GetFockOffDiagElement(obj, atomA, atomB, ...
                indexAtomA, indexAtomB, mu, nu, ~, gammaAB, ... 
                overlapAOs, orbitalElectronPopulation, ~, isGuess)
            K = obj.GetBondingAdjustParameterK(atomA.GetValenceShellType(), atomB.GetValenceShellType());
            bondParameter = 0.5*K*(atomA.GetBondingParameter() + atomB.GetBondingParameter());
            value =  bondParameter*overlapAOs(mu, nu);
            if(~isGuess)
                value = value - 0.5*orbitalElectronPopulation(mu, nu)*gammaAB(indexAtomA, indexAtomB);
            end
        end
        
        %    void TransposeFockMatrixMatrix(double** transposedFockMatrix) const;
        
        function diatomicOverlapAOs = CalcDiatomicOverlapAOsInDiatomicFrame(obj, atomA, atomB)
            na = double(atomA.GetValenceShellType());
            nb = double(atomB.GetValenceShellType());
            
            diatomicOverlapAOs = zeros(double(OrbitalType.OrbitalType_end));
            rAB = obj.molecule.GetDistanceAtoms(atomA, atomB); % Inter nuclear distance between aton A and B.
            
            for a = 1:atomA.GetValenceSize()
                valenceOrbitalA = atomA.GetValence(a);
                realShpericalHarmonicsA = atomA.GetRealSphericalHarmonicsIndex(a);
                orbitalExponentA = atomA.GetOrbitalExponent(atomA.GetValenceShellType(), valenceOrbitalA, obj.theory);

                for b = 1:atomB.GetValenceSize()
                    valenceOrbitalB = atomB.GetValence(b);
                    realShpericalHarmonicsB = atomB.GetRealSphericalHarmonicsIndex(b);
                    orbitalExponentB = atomB.GetOrbitalExponent(atomB.GetValenceShellType(), valenceOrbitalB, obj.theory);

                    if(realShpericalHarmonicsA.GetM() == realShpericalHarmonicsB.GetM())
                        m = abs(realShpericalHarmonicsA.GetM());
                        alpha = orbitalExponentA * rAB;
                        beta =  orbitalExponentB * rAB;

                        reducedOverlapAOs = obj.GetReducedOverlapAOs(na, realShpericalHarmonicsA.GetL(), m,...
                            nb, realShpericalHarmonicsB.GetL(), alpha, beta);


                        pre =  power(2.0*orbitalExponentA, na+0.5);
                        pre = pre * power(2.0*orbitalExponentB, nb+0.5);
                        factorials = factorial(2*na)*factorial(2*nb);
                        pre = pre / sqrt(factorials);
                        pre = pre * power(rAB/2.0, na+nb+1.0);

                        diatomicOverlapAOs(valenceOrbitalA, valenceOrbitalB) = pre*reducedOverlapAOs;
                    end
                end
            end
        end
        
        %    virtual void CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(double** diatomicOverlapAOsDeri,
        %                                                                    const MolDS_base_atoms::Atom& atomA,
        %                                                                    const MolDS_base_atoms::Atom& atomB) const;
        %    virtual void CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(double** diatomicOverlapAOs2ndDeri,
        %                                                                    const MolDS_base_atoms::Atom& atomA,
        %                                                                    const MolDS_base_atoms::Atom& atomB) const;
        %    void CalcDiatomicOverlapAOs1stDerivatives(double*** diatomicOverlapAOs1stDerivs,
        %                                              double**  tmpDiaOverlapAOsInDiaFrame,
        %                                              double**  tmpDiaOverlapAOs1stDerivInDiaFrame,
        %                                              double**  tmpRotMat,
        %                                              double**  tmpRotMat1stDeriv,
        %                                              double*** tmpRotMat1stDerivs,
        %                                              double**  tmpRotatedDiatomicOverlap,
        %                                              double*   tmpRotatedDiatomicOverlapVec,
        %                                              double**  tmpMatrixBC,
        %                                              double*   tmpVectorBC,
        %                                              const MolDS_base_atoms::Atom& atomA,
        %                                              const MolDS_base_atoms::Atom& atomB) const;
        %    void CalcDiatomicOverlapAOs1stDerivatives(double*** diatomicOverlapAOs1stDerivs,
        %                                              double**  tmpDiaOverlapAOsInDiaFrame,
        %                                              double**  tmpDiaOverlapAOs1stDerivInDiaFrame,
        %                                              double**  tmpRotMat,
        %                                              double**  tmpRotMat1stDeriv,
        %                                              double*** tmpRotMat1stDerivs,
        %                                              double**  tmpRotatedDiatomicOverlap,
        %                                              double*   tmpRotatedDiatomicOverlapVec,
        %                                              double**  tmpMatrixBC,
        %                                              double*   tmpVectorBC,
        %                                              int indexAtomA,
        %                                              int indexAtomB) const;
        %    void CalcDiatomicOverlapAOs2ndDerivatives(double**** overlapAOs2ndDeri,
        %                                              double**   tmpDiaOverlapAOsInDiaFrame,
        %                                              double**   tmpDiaOverlapAOs1stDerivInDiaFrame,
        %                                              double**   tmpDiaOverlapAOs2ndDerivInDiaFrame,
        %                                              double***  tmpDiaOverlapAOs1stDerivs,
        %                                              double**** tmpDiaOverlapAOs2ndDerivs,
        %                                              double**   tmpRotMat,
        %                                              double***  tmpRotMat1stDerivs,
        %                                              double**** tmpRotMat2ndDerivs,
        %                                              const MolDS_base_atoms::Atom& atomA,
        %                                              const MolDS_base_atoms::Atom& atomB) const;
        %    void CalcDiatomicOverlapAOs2ndDerivatives(double**** overlapAOs2ndDeri,
        %                                              double**   tmpDiaOverlapAOsInDiaFrame,
        %                                              double**   tmpDiaOverlapAOs1stDerivInDiaFrame,
        %                                              double**   tmpDiaOverlapAOs2ndDerivInDiaFrame,
        %                                              double***  tmpDiaOverlapAOs1stDerivs,
        %                                              double**** tmpDiaOverlapAOs2ndDerivs,
        %                                              double**   tmpRotMat,
        %                                              double***  tmpRotMat1stDerivs,
        %                                              double**** tmpRotMat2ndDerivs,
        %                                              int indexAtomA,
        %                                              int indexAtomB) const;
        %    double Get2ndDerivativeElementFromDistanceDerivatives(double firstDistanceDeri,
        %                                                          double secondDistanceDeri,
        %                                                          MolDS_base::CartesianType axisA1,
        %                                                          MolDS_base::CartesianType axisA2,
        %                                                          double* cartesian,
        %                                                          double rAB) const;
        %    virtual double GetMolecularIntegralElement(int moI, int moJ, int moK, int moL,
        %                                               const MolDS_base::Molecule& molecule,
        %                                               double const* const* fockMatrix,
        %                                               double const* const* gammaAB) const;
        
        function [twoElecsTwoAtomCores, twoElecsAtomEpcCores] = CalcTwoElecsTwoCores(~, ~)
            % do nothing for CNDO, INDO, and ZINDO/S.
            % two electron two core integrals are not needed for CNDO, INDO, and ZINDO/S.
            twoElecsTwoAtomCores = [];
            twoElecsAtomEpcCores = [];
        end
        
        %    virtual void CalcForce(const std::vector<int>& elecStates);
        %    void CalcRotatingMatrix1stDerivatives(double*** rotMat1stDerivatives,
        %                                          const MolDS_base_atoms::Atom& atomA,
        %                                          const MolDS_base_atoms::Atom& atomB) const;
        %    void CalcRotatingMatrix2ndDerivatives(double**** rotMat2ndDerivatives,
        %                                          const MolDS_base_atoms::Atom& atomA,
        %                                          const MolDS_base_atoms::Atom& atomB) const;
        
    end
    
    methods (Access = private)
        
        function CalcCoreRepulsionEnergy(obj)
            energy = 0.0;
            atomvect = obj.molecule.GetAtomVect();
            for i = 1:length(atomvect)
                for j = i+1:length(atomvect)
                    energy = energy + obj.GetDiatomCoreRepulsionEnergy(atomvect{i}, atomvect{j});
                end
            end
            obj.coreRepulsionEnergy = energy;
        end
        
        function SetEnableAtomTypesVdW(obj)
            obj.enableAtomTypesVdW = {};
            obj.enableAtomTypesVdW{end+1} = AtomType.H;
            obj.enableAtomTypesVdW{end+1} = AtomType.C;
            obj.enableAtomTypesVdW{end+1} = AtomType.N;
            obj.enableAtomTypesVdW{end+1} = AtomType.O;
            obj.enableAtomTypesVdW{end+1} = AtomType.F;
            obj.enableAtomTypesVdW{end+1} = AtomType.S;
            obj.enableAtomTypesVdW{end+1} = AtomType.Cl;
        end
        
        function CheckEnableAtomTypeVdW(obj, molecule)
            if(Parameters.GetInstance().RequiresVdWSCF())
                if(obj.theory == TheoryType.AM1D || obj.theory == TheoryType.PM3D)
                    return;
                end
            else
                return;
            end
            for i = 1:length(molecule.GetAtomVect())
                atomvect = molecule.GetAtomVect();
                atomType = atomvect{i}.GetAtomType();
                enable = false;
                for j = 1:length(obj.enableAtomTypesVdW)
                    if(atomType == obj.enableAtomTypesVdW{j})
                        enable = true;
                        break;
                    end
                end
                if(~enable)
                    throw(MException('Cndo2:CheckEnableAtomTypeVdW', 'Atom type not enabled for this theory.'));
                end
            end
        end
        
        function CalcVdWCorrectionEnergy(obj)
            value = 0.0;
            atomvect = obj.molecule.GetAtomVect();
            for i = 1:length(atomvect)
                atomI = atomvect{i};
                for j = i+1:length(atomvect)
                    atomJ = atomvect{j};
                    value = value + obj.GetDiatomVdWCorrectionEnergy(atomI, atomJ);
                end
            end
            obj.vdWCorrectionEnergy = value;
        end
        
        function res = GetVdWDampingValue(~, vdWDistance, distance)
            dampingFactor = Parameters.GetInstance().GetVdWDampingFactorSCF();
            res = 1.0/(1.0+exp(-1.0*dampingFactor*(distance/vdWDistance - 1.0)));
        end
        
        %    double GetVdWDampingValue1stDerivative(double vdWDistance, double distance) const;
        %    double GetVdWDampingValue2ndDerivative(double vdWDistance, double distance) const;
        
        function electronicTransitionDipoleMoments = CalcElectronicDipoleMomentGroundState(obj, ...
                cartesianMatrix, molecule, orbitalElectronPopulation, overlapAOs)
            groundState = 1;
            electronicTransitionDipoleMoments(groundState, groundState, :) = obj.CalcElectronicTransitionDipoleMoment(...
                groundState,...
                groundState,...
                0,...
                0,...
                cartesianMatrix,...
                molecule,...
                orbitalElectronPopulation,...
                overlapAOs,...
                0);
        end
        
        function satisfy = SatisfyConvergenceCriterion(~, oldOrbitalElectronPopulation, ...
                orbitalElectronPopulation, numberAOs)
            satisfy = false;
            change = 0.0;
            for i = 1:numberAOs
                for j = 1:numberAOs
                    change = change + (oldOrbitalElectronPopulation(i,j) - orbitalElectronPopulation(i,j))...
                        *(oldOrbitalElectronPopulation(i,j) - orbitalElectronPopulation(i,j));
                end
            end
            
            change = change / numberAOs*numberAOs;
            rmsDensity = sqrt(change);
            
            if(rmsDensity < Parameters.GetInstance().GetThresholdSCF())
                satisfy = true;
            end
        end
        
        %    void CalcOrbitalElectronPopulation(double** orbitalElectronPopulation,
        %                                       const MolDS_base::Molecule& molecule,
        %                                       double const* const* fockMatrix) const;
        
        function atomicElectronPopulation = CalcAtomicElectronPopulation(~, orbitalElectronPopulation, molecule)
            atomvect = molecule.GetAtomVect();
            totalNumberAtoms = length(atomvect);
            atomicElectronPopulation = zeros(1, totalNumberAtoms);
            for A = 1:totalNumberAtoms
                firstAOIndex = atomvect{A}.GetFirstAOIndex();
                numberAOs    = atomvect{A}.GetValenceSize();
                for i = firstAOIndex:firstAOIndex+numberAOs-1
                    atomicElectronPopulation(A) = atomicElectronPopulation(A) + orbitalElectronPopulation(i, i);
                end
            end
        end
        
        function coreDipoleMoment = CalcCoreDipoleMoment(~, molecule)
             dipoleCenter = molecule.GetXyzDipoleCenter();
             atomvect = molecule.GetAtomVect();
             coreDipoleMoment = zeros(1, 3);
             for A = 1:length(atomvect)
                 coreDipoleMoment = coreDipoleMoment + atomvect{A}.GetCoreCharge()...
                     .*(atomvect{A}.GetXyz() - dipoleCenter);
             end
        end
        
        function cartesianMatrix = CalcCartesianMatrixByGTOExpansion(obj, molecule, stonG)
            atomvect = molecule.GetAtomVect();
            totalAtomNumber = length(atomvect);
            for A = 1:totalAtomNumber
                atomA  = atomvect{A};
                firstAOIndexA  = atomA.GetFirstAOIndex();
                numValenceAOsA = atomA.GetValenceSize();
                for a = 1:numValenceAOsA
                    mu = firstAOIndexA + a - 1;
                    for B = 1:totalAtomNumber
                        atomB  = atomvect{B};
                        firstAOIndexB  = atomB.GetFirstAOIndex();
                        numValenceAOsB = atomB.GetValenceSize();
                        for b = 1:numValenceAOsB
                            nu = firstAOIndexB + b - 1;
                            cartesianMatrix(:, mu, nu) = obj.CalcCartesianMatrixElementsByGTOExpansion(atomA, a, atomB, b, stonG);
                        end
                    end 
                end
            end
        
        end
        
        function xyz = CalcCartesianMatrixElementsByGTOExpansion(obj, atomA, valenceIndexA, atomB, valenceIndexB, stonG)
            xyz = zeros(1,3);
            shellTypeA = atomA.GetValenceShellType();
            shellTypeB = atomB.GetValenceShellType();
            valenceOrbitalA = atomA.GetValence(valenceIndexA);
            valenceOrbitalB = atomB.GetValence(valenceIndexB);
            orbitalExponentA = atomA.GetOrbitalExponent(atomA.GetValenceShellType(), valenceOrbitalA, obj.theory);
            orbitalExponentB = atomB.GetOrbitalExponent(atomB.GetValenceShellType(), valenceOrbitalB, obj.theory);
            dxyz = atomA.GetXyz() - atomB.GetXyz();
            rAB = norm(dxyz);
            for i = 1:double(stonG)
                for j = 1:double(stonG)
                    temp = GTOExpansionSTO.GetInstance().GetCoefficient(stonG, ...
                        shellTypeA, ...
                        valenceOrbitalA, ...
                        i);
                    temp = temp * GTOExpansionSTO.GetInstance().GetCoefficient(stonG, ...
                        shellTypeB, ...
                        valenceOrbitalB, ...
                        j);
                    gaussianExponentA = (orbitalExponentA*orbitalExponentA) *...
                        GTOExpansionSTO.GetInstance().GetExponent(stonG, ...
                        shellTypeA, ...
                        valenceOrbitalA, ...
                        i);
                    gaussianExponentB = (orbitalExponentB*orbitalExponentB) *...
                        GTOExpansionSTO.GetInstance().GetExponent(stonG, ...
                        shellTypeB, ...
                        valenceOrbitalB, ...
                        j);
                    overlapSASB = obj.GetGaussianOverlapAOsSASB(gaussianExponentA, gaussianExponentB, rAB);
                    tempX = obj.GetGaussianCartesianMatrix(atomA.GetAtomType(), valenceOrbitalA, gaussianExponentA, atomA.GetXyz(),...
                        atomB.GetAtomType(), valenceOrbitalB, gaussianExponentB, atomB.GetXyz(), ...
                        rAB, overlapSASB, 1);
                    tempY = obj.GetGaussianCartesianMatrix(atomA.GetAtomType(), valenceOrbitalA, gaussianExponentA, atomA.GetXyz(),...
                        atomB.GetAtomType(), valenceOrbitalB, gaussianExponentB, atomB.GetXyz(), ...
                        rAB, overlapSASB, 2);
                    tempZ = obj.GetGaussianCartesianMatrix(atomA.GetAtomType(), valenceOrbitalA, gaussianExponentA, atomA.GetXyz(),...
                        atomB.GetAtomType(), valenceOrbitalB, gaussianExponentB, atomB.GetXyz(), ...
                        rAB, overlapSASB, 3);
                    xyz = xyz + temp.*[tempX, tempY, tempZ];
                end
            end
        end
        
        value = GetGaussianCartesianMatrix(obj, ...
                atomTypeA, valenceOrbitalA, gaussianExponentA, xyzA,...
                atomTypeB, valenceOrbitalB, gaussianExponentB, xyzB,...
                rAB, overlapSASB, axis)
        
        function overlapAOs = CalcOverlapAOs(obj, molecule)
            totalAONumber = molecule.GetTotalNumberAOs();
            totalAtomNumber = length(molecule.GetAtomVect());
            overlapAOs = eye(totalAONumber);
            atomvect = molecule.GetAtomVect();
            for A = 1:totalAtomNumber
                atomA = atomvect{A};
                symmetrize = false;
                for B = A+1:totalAtomNumber
                    atomB = atomvect{B};
                    diatomicOverlapAOs = obj.CalcDiatomicOverlapAOsInDiatomicFrame(atomA, atomB);
                    rotatingMatrix = obj.CalcRotatingMatrix(atomA, atomB);
                    diatomicOverlapAOs = obj.RotateDiatmicOverlapAOsToSpaceFrame(diatomicOverlapAOs, rotatingMatrix);
                    overlapAOs = obj.SetOverlapAOsElement(overlapAOs, diatomicOverlapAOs, atomA, atomB, symmetrize);
                end
            end
            overlapAOs = overlapAOs + overlapAOs' - diag(diag(overlapAOs));
        end
        
        %    void CalcOverlapAOsByGTOExpansion(double** overlapAOs,
        %                                      const MolDS_base::Molecule& molecule,
        %                                      MolDS_base::STOnGType stonG) const; //See [DY_1977]
        %    double GetOverlapAOsElementByGTOExpansion(const MolDS_base_atoms::Atom& atomA,
        %                                              int valenceIndexA,
        %                                              const MolDS_base_atoms::Atom& atomB,
        %                                              int valenceIndexB,
        %                                              MolDS_base::STOnGType stonG) const; // see [DY_1977]
        
        function value = GetGaussianOverlapAOsSASB(~, gaussianExponentA, gaussianExponentB, rAB)
            gauPlusAB = gaussianExponentA+gaussianExponentB;
            gauMultAB = gaussianExponentA*gaussianExponentB;
            temp1 = 2.0*sqrt(gauMultAB)/gauPlusAB;
            temp2 = -1.0* gauMultAB/gauPlusAB;
            value = temp1*sqrt(temp1)*exp(temp2*rAB*rAB);
        end
        
        value = GetGaussianOverlapAOs(~, ~, valenceOrbitalA, gaussianExponentA,...
                ~, valenceOrbitalB,gaussianExponentB,...
                dxyz, ~, overlapSASB)
        
        %    double GetGaussianOverlapAOs1stDerivative(MolDS_base::AtomType atomTypeA,
        %                                              MolDS_base::OrbitalType valenceOrbitalA,
        %                                              double gaussianExponentA,
        %                                              MolDS_base::AtomType atomTypeB,
        %                                              MolDS_base::OrbitalType valenceOrbitalB,
        %                                              double gaussianExponentB,
        %                                              double dx,
        %                                              double dy,
        %                                              double dz,
        %                                              double rAB,
        %                                              MolDS_base::CartesianType axisA) const;// see [DY_1977]
        
        function fockMatrix = CalcFockMatrix(obj, molecule, overlapAOs, ...
                gammaAB, orbitalElectronPopulation, atomicElectronPopulation, ...
                twoElecsTwoAtomCores, isGuess)
            totalNumberAOs   = molecule.GetTotalNumberAOs();
            atomvect = molecule.GetAtomVect();
            totalNumberAtoms = length(atomvect);
            fockMatrix = zeros(totalNumberAOs);
            
            for A = 1:totalNumberAtoms
                atomA = atomvect{A};
                firstAOIndexA = atomA.GetFirstAOIndex();
                lastAOIndexA  = atomA.GetLastAOIndex();
                for mu = firstAOIndexA:lastAOIndexA
                    for B = A:totalNumberAtoms
                        atomB = atomvect{B};
                        firstAOIndexB = atomB.GetFirstAOIndex();
                        lastAOIndexB  = atomB.GetLastAOIndex();
                        for nu = firstAOIndexB:lastAOIndexB
                            if(mu == nu)
                                % diagonal part
                                fockMatrix(mu, mu) = obj.GetFockDiagElement(atomA, ...
                                    A, ...
                                    mu, ...
                                    molecule, ...
                                    gammaAB,...
                                    orbitalElectronPopulation, ...
                                    atomicElectronPopulation,...
                                    twoElecsTwoAtomCores,...
                                    isGuess);
                            elseif(mu < nu)
                                % upper right part
                                fockMatrix(mu, nu) = obj.GetFockOffDiagElement(atomA, ...
                                    atomB,...
                                    A, ...
                                    B, ...
                                    mu, ...
                                    nu, ...
                                    molecule, ...
                                    gammaAB,...
                                    overlapAOs,...
                                    orbitalElectronPopulation, ...
                                    twoElecsTwoAtomCores,...
                                    isGuess);
                                % whether to have this?
                                fockMatrix(nu, mu) = fockMatrix(mu, nu);
                            else
                                % lower left part (not calculated)
                            end
                        end
                    end
                end % end of loop mu
            end % end of loop A
        end
        
        function res = RotateDiatmicOverlapAOsToSpaceFrame(~, diatomicOverlapAOs, rotatingMatrix)
            res = rotatingMatrix * diatomicOverlapAOs * rotatingMatrix;
        end
        
        function overlapAOs = SetOverlapAOsElement(~, overlapAOs, diatomicOverlapAOs, atomA, atomB, isSymmetricOverlapAOs)
            if(nargin < 5)
                isSymmetricOverlapAOs = true;
            end
            firstAOIndexAtomA = atomA.GetFirstAOIndex();
            firstAOIndexAtomB = atomB.GetFirstAOIndex();
            
            for i = 1:atomA.GetValenceSize()
                orbitalA = double(atomA.GetValence(i));
                for j = 1:atomB.GetValenceSize()
                    orbitalB = double(atomB.GetValence(j));
                    mu = firstAOIndexAtomA + i - 1;
                    nu = firstAOIndexAtomB + j - 1;
                    overlapAOs(mu,nu) = diatomicOverlapAOs(orbitalA,orbitalB);
                    if(isSymmetricOverlapAOs)
                        overlapAOs(nu,mu) = diatomicOverlapAOs(orbitalA,orbitalB);
                    end
                end
            end
        end
        
        function value = GetAuxiliaryA(~, k, rho)
            tmp1 = 0.0;
            value = exp(-1.0*rho)*factorial(k);
            for mu = 1:k+1
                tmp2=rho^mu; % tmp2 = pow(rho,mu)
                tmp1 = tmp1 + 1.0/(tmp2*factorial(k-mu+1));
            end
            value = value * tmp1;
        end
        
        function value = GetAuxiliaryB(~, k, rho)
            tmp1 = 0.0;
            tmp2 = 0.0;
            if(abs(rho)>0)
                pre1 = -1.0*exp(-1.0*rho);
                pre2 = -1.0*exp(rho);
                for mu = 1:k+1
                    tmp3 = rho^mu; % tmp3 = pow(rho,mu)
                    tmp4 = -1.0; % tmp4 = pow(-1.0,k-mu)
                    if(mod((k-mu),2)==0)
                        tmp4=1.0;
                    end
                    tmp1 = tmp1 + factorial(k)/factorial(k-mu+1)     /tmp3;
                    tmp2 = tmp2 + factorial(k)/factorial(k-mu+1)*tmp4/tmp3;
                end
                value = pre1*tmp1 + pre2*tmp2;
            else
                if(mod(k,2) == 0)
                    value = 2.0/(1.0+k);
                else
                    value = 0;
                end
            end
        end
        
        function value = GetAuxiliaryD(~, la, lb, m)
            if(m<0)
                throw(MException('Cndo2:GetAuxiliaryD', 'm < 0'));
            end
            tmp = factorial(m+1)/8.0;
            pre = tmp*tmp;
            termA = ( (2.0*la+1.0)*factorial(la-m) ) / ( 2.0*factorial(la+m) );
            termB = ( (2.0*lb+1.0)*factorial(lb-m) ) / ( 2.0*factorial(lb+m) );
            value = pre*sqrt(termA)*sqrt(termB);
        end
        
        %    double GetAuxiliaryA1stDerivative(int k, double rho) const;
        %    double GetAuxiliaryA2ndDerivative(int k, double rho) const;
        %    double GetAuxiliaryB1stDerivative(int k, double rho) const;
        %    double GetAuxiliaryB2ndDerivative(int k, double rho) const;
        %    void DoDamp(double rmsDensity,
        %                bool&  hasAppliedDamping,
        %                double** orbitalElectronPopulation,
        %                double const* const* oldOrbitalElectronPopulation,
        %                const MolDS_base::Molecule& molecule) const;
        %    void DoDIIS(double** orbitalElectronPopulation,
        %                double const* const* oldOrbitalElectronPopulation,
        %                double*** diisStoredDensityMatrix,
        %                double*** diisStoredErrorVect,
        %                double**  diisErrorProducts,
        %                double**  tmpDiisErrorProducts,
        %                double*   diisErrorCoefficients,
        %                double&   diisError,
        %                bool&     hasAppliedDIIS,
        %                int       diisNumErrorVect,
        %                const MolDS_base::Molecule& molecule,
        %                int step) const;
        
        function CheckEnableAtomType(obj, molecule)
            for i = 1:length(molecule.GetAtomVect())
                atomvect = molecule.GetAtomVect();
                atomType = atomvect{i}.GetAtomType();
                enable = false;
                for j = 1:length(obj.enableAtomTypes)
                    if(atomType == obj.enableAtomTypes{j})
                        enable = true;
                        break;
                    end
                end
                if(~enable)
                    throw(MException('Cndo2:CheckEnableAtomType', 'Atom type not enabled for this theory.'));
                end
            end
        end
        
        function CheckNumberValenceElectrons(~, molecule)
            if(mod(molecule.GetTotalNumberValenceElectrons(), 2) == 1)
                throw(MException('Cndo2:CheckNumberValenceElectrons', 'Odd number of electrons.'));
            end
        end
        
        %    void FreeDiatomicOverlapAOsAndRotatingMatrix(double*** diatomicOverlapAOs,
        %                                                 double*** rotatingMatrix) const;
        
        function elecSCFEnergy = CalcElecSCFEnergy(obj, molecule,...
                ~, ~, ~, coreRepulsionEnergy,...
                coreEpcCoulombEnergy, vdWCorrectionEnergy)
            % use density matrix for electronic energy
            totalNumberAOs = obj.molecule.GetTotalNumberAOs();
            
            dammyOrbitalElectronPopulation = zeros(totalNumberAOs);
            dammyAtomicElectronPopulation = zeros(1, length(molecule.GetAtomVect()));
            
            isGuess = false;
            fMatrix = obj.CalcFockMatrix(molecule, ...
                obj.overlapAOs, ...
                obj.gammaAB,...
                obj.orbitalElectronPopulation, ...
                obj.atomicElectronPopulation,...
                obj.twoElecsTwoAtomCores,...
                isGuess);
            hMatrix = obj.CalcFockMatrix(molecule, ...
                obj.overlapAOs, ...
                obj.gammaAB,...
                dammyOrbitalElectronPopulation, ...
                dammyAtomicElectronPopulation,...
                obj.twoElecsTwoAtomCores,...
                isGuess);
            electronicEnergy = sum(sum(obj.orbitalElectronPopulation .* (hMatrix+fMatrix)));
            electronicEnergy = electronicEnergy * 0.5;
            
            elecSCFEnergy = electronicEnergy...
                +coreRepulsionEnergy...
                +coreEpcCoulombEnergy...
                +vdWCorrectionEnergy;
        end
        
        
        
    end
    
end




