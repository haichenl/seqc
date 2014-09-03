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
        bondingAdjustParameterK = zeros(2,1); %see (3.79) in J. A. Pople book
        gammaAB;
        enableAtomTypesVdW = {};
        
        ReducedOverlapAOsParameters_Z;
        ReducedOverlapAOsParameters_Y;
   
    end
    
    methods (Access = public)
        
        function obj = Cndo2()
            Arguments.GetInstance().SetCurrentTheory(EnumTheory.CNDO2);
            
            %protected variables
            obj.molecule = [];
            obj.theory = EnumTheory.CNDO2;
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
            obj.CheckNumberValenceElectrons();
            obj.CheckEnableAtomType();
            obj.CheckEnableAtomTypeVdW();
            nbf = obj.molecule.totalNumberAOs;
            obj.fockMatrix = zeros(nbf);
            obj.energiesMO = zeros(nbf, 1);
            obj.orbitalElectronPopulation = zeros(nbf);
            obj.atomicElectronPopulation = zeros(length(obj.molecule.atomVect), 1);
            obj.overlapAOs = zeros(nbf);
            obj.cartesianMatrix = zeros(nbf, nbf, 3);
            electronicTransitionDipoleMomentsDim = 1;
%             if(Parameters::GetInstance()->RequiresCIS()){
%                 electronicTransitionDipoleMomentsDim += Parameters::GetInstance()->GetNumberExcitedStatesCIS();
%                 }
            obj.electronicTransitionDipoleMoments = zeros(electronicTransitionDipoleMomentsDim);
            obj.coreDipoleMoment = zeros(3, 1);
            
            % calculate electron integral
            if(uint8(obj.theory) == uint8(EnumTheory.CNDO2) || uint8(obj.theory) == uint8(EnumTheory.INDO))
                obj.gammaAB = obj.CalcGammaAB();
            end
            obj.overlapAOs = obj.CalcOverlapAOs();
            obj.cartesianMatrix = obj.CalcCartesianMatrixByGTOExpansion(uint8(EnumSTOnG.STO6G));
%             [obj.twoElecsTwoAtomCores, obj.twoElecsAtomEpcCores] = obj.CalcTwoElecsTwoCores(obj.molecule);
            
        end
        
        function DoSCF(obj)
            requiresGuess = true;
            
            % Cndo2::MallocSCFTemporaryMatrices
            diisNumErrorVect = Arguments.GetInstance().diisNumErrorVectSCF;
            if(0<diisNumErrorVect)
                diisStoredDensityMatrix = zeros(obj.molecule.totalNumberAOs* obj.molecule.totalNumberAOs, diisNumErrorVect);
                diisStoredErrorVect = zeros(obj.molecule.totalNumberAOs()* obj.molecule.totalNumberAOs, diisNumErrorVect);
                diisErrorProducts = zeros(diisNumErrorVect+1);
                diisErrorCoefficients = zeros(diisNumErrorVect+1, 1);
            end
            
            % SCF
            maxIterationsSCF = Arguments.GetInstance().maxIterationsSCF;
            for iterationStep = 0:maxIterationsSCF-1
                obj.atomicElectronPopulation = obj.CalcAtomicElectronPopulation();
                oldOrbitalElectronPopulation = obj.orbitalElectronPopulation;
                
                isGuess = (iterationStep==0 && requiresGuess);
                % continue here
                obj.fockMatrix = obj.CalcFockMatrix(isGuess);
                
                % diagonalization of the Fock matrix
                [orbital, eigenvalue] = eig(obj.fockMatrix);
                [obj.energiesMO, order] = sort(diag(eigenvalue));
                orbital = orbital(:, order);
                numberTotalValenceElectrons = obj.molecule.totalNumberValenceElectrons;
                obj.orbitalElectronPopulation = 2 .* orbital(:, 1:numberTotalValenceElectrons/2) ...
                    *orbital(:, 1:numberTotalValenceElectrons/2)';
                
                % check convergence
                hasConverged = obj.SatisfyConvergenceCriterion(oldOrbitalElectronPopulation, ...
                    obj.orbitalElectronPopulation);
                
                if(hasConverged)
                    obj.CalcSCFProperties();
                    break;
                else
                    if(~isGuess)
                        [obj.orbitalElectronPopulation, ...
                            diisStoredDensityMatrix,...
                            diisStoredErrorVect,...
                            diisErrorProducts,...
                            diisErrorCoefficients] = obj.DoDIIS(obj.orbitalElectronPopulation,...
                            oldOrbitalElectronPopulation,...
                            diisStoredDensityMatrix,...
                            diisStoredErrorVect,...
                            diisErrorProducts,...
                            diisErrorCoefficients,...
                            diisNumErrorVect,...
                            iterationStep);
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
        
%         function res = GetEnumTheory(obj)
%             res = obj.theory;
%         end
        
    end
    
    methods (Access = protected)
        
        function SetEnableAtomTypes(obj)
            obj.enableAtomTypes = {};
            obj.enableAtomTypes{end+1} = EnumAtom.H;
            obj.enableAtomTypes{end+1} = EnumAtom.Li;
            %obj.enableAtomTypes{end+1} = EnumAtom.Be;
            %obj.enableAtomTypes{end+1} = EnumAtom.B;
            obj.enableAtomTypes{end+1} = EnumAtom.C;
            obj.enableAtomTypes{end+1} = EnumAtom.N;
            obj.enableAtomTypes{end+1} = EnumAtom.O;
            obj.enableAtomTypes{end+1} = EnumAtom.F;
            %obj.enableAtomTypes{end+1} = EnumAtom.Na;
            %obj.enableAtomTypes{end+1} = EnumAtom.Mg;
            %obj.enableAtomTypes{end+1} = EnumAtom.Al;
            %obj.enableAtomTypes{end+1} = EnumAtom.Si;
            %obj.enableAtomTypes{end+1} = EnumAtom.P;
            obj.enableAtomTypes{end+1} = EnumAtom.S;
            obj.enableAtomTypes{end+1} = EnumAtom.Cl;
        end
        
        function CalcSCFProperties(obj)
            obj.atomicElectronPopulation = obj.CalcAtomicElectronPopulation();
            obj.CalcCoreRepulsionEnergy();
            if(Arguments.GetInstance().requiresVdWSCF)
                obj.CalcVdWCorrectionEnergy();
            end
            obj.elecSCFEnergy = obj.CalcElecSCFEnergy();
            obj.coreDipoleMoment = obj.CalcCoreDipoleMoment();
            obj.electronicTransitionDipoleMoments = obj.CalcElectronicDipoleMomentGroundState();
            %    groundState = 0;
            %    if(Parameters::GetInstance()->RequiresFrequencies() &&
            %       Parameters::GetInstance()->GetElectronicStateIndexFrequencies() == groundState){
            %       this->CalcNormalModes(this->normalModes, this->normalForceConstants, *this->molecule);
            %    }
        end
        
        %         CalcNormalModes(double** normalModes, double* normalForceConstants, const MolDS_base::Molecule& molecule) const;
        
        function transitionDipoleMoment = CalcElectronicTransitionDipoleMoment(obj, to, from)
            groundState = 1;
            if(from == groundState && to == groundState)
                dipoleCenter = obj.molecule.GetXyzDipoleCenter();
                transitionDipoleMoment = zeros(3, 1);
                transitionDipoleMoment(1) = transitionDipoleMoment(1) - sum(sum(obj.orbitalElectronPopulation.*obj.cartesianMatrix(:,:,1)));
                transitionDipoleMoment(2) = transitionDipoleMoment(2) - sum(sum(obj.orbitalElectronPopulation.*obj.cartesianMatrix(:,:,2)));
                transitionDipoleMoment(3) = transitionDipoleMoment(3) - sum(sum(obj.orbitalElectronPopulation.*obj.cartesianMatrix(:,:,3)));
                
                % set orign of dipole
                temp = sum(sum(obj.orbitalElectronPopulation.*obj.overlapAOs));
                transitionDipoleMoment = transitionDipoleMoment + dipoleCenter.*temp;
            else
                throw(MException('Cndo2:CalcElectronicTransitionDipoleMoment', 'Excited transition moment not enabled for this theory.'));
            end
        end
        
        function value = GetBondingAdjustParameterK(obj, shellA, shellB)
             value=1.0;
             if(shellA >= 3 || shellB >= 3)
                 value = obj.bondingAdjustParameterK(2); % notice c++->matlab index conversion
             end
        end
        
        %         GetAtomCoreEpcCoulombEnergy (const MolDS_base_atoms::Atom& atom,
        %                                                const MolDS_base_atoms::Atom& epc) const;
        
        function res = GetDiatomCoreRepulsionEnergy(obj, atomA, atomB)
            distance = obj.molecule.GetDistanceAtoms(atomA, atomB);
            res = atomA.coreCharge*atomB.coreCharge/distance;
        end
        
        % not tested
        function valueVec = GetDiatomCoreRepulsion1stDerivative(obj, atomA, atomB) % valueVec 1 x 3
            distance = obj.molecule.GetDistanceAtoms(atomA, atomB);
            valueVec = atomA.coreCharge*atomB.coreCharge;
            valueVec = valueVec .* ((atomA.xyz - atomB.xyz)/distance);
            valueVec = valueVec .* (-1.0/(distance*distance));
        end
        
        %    virtual double GetDiatomCoreRepulsion2ndDerivative(const MolDS_base_atoms::Atom& atomA,
        %                                                       const MolDS_base_atoms::Atom& atomB,
        %                                                       MolDS_base::CartesianType axisA1,
        %                                                       MolDS_base::CartesianType axisA2) const;
        
        % not tested
        function res = GetDiatomVdWCorrectionEnergy(obj, atomA, atomB)
            distance = obj.molecule.GetDistanceAtoms(atomA, atomB);
            vdWDistance = atomA.vdWRadii + atomB.vdWRadii;
            tmpSum = atomA.vdWCoefficient+atomB.vdWCoefficient;
            if(tmpSum<=0)
                res = 0;
                return;
            end
            vdWCoefficients = 2.0*atomA.vdWCoefficient*atomB.vdWCoefficient/tmpSum;
            damping = obj.GetVdWDampingValue(vdWDistance, distance);
            scalingFactor = Arguments.GetInstance().vdWScalingFactorSCF;
            res =  -1.0*scalingFactor*vdWCoefficients*damping ...
                /(distance^6);
        end
        % not tested
        function valueVec = GetDiatomVdWCorrection1stDerivative(obj, atomA, atomB)
            distance = obj.molecule.GetDistanceAtoms(atomA, atomB);
            vdWDistance = atomA.vdWRadii + atomB.vdWRadii;
            tmpSum = atomA.vdWCoefficient+atomB.vdWCoefficient;
            if(tmpSum<=0e0)
                valueVec = 0;
                return;
            end
            vdWCoefficients = 2.0*atomA.vdWCoefficient*atomB.vdWCoefficient/tmpSum;
%             dampingFactor = Arguments.GetInstance().vdWDampingFactorSCF;
            damping = obj.GetVdWDampingValue(vdWDistance, distance);
            damping1stDerivative = obj.GetVdWDampingValue1stDerivative(vdWDistance, distance);
            valueVec = zeros(3, 1);
            tmp = distance ^ 6; %distance*distance*distance*distance*distance*distance;
            valueVec = valueVec + 6.0*damping/(tmp*distance) ...
                -damping1stDerivative/tmp;
            valueVec = valueVec * vdWCoefficients;
            valueVec = valueVec * Arguments.GetInstance().vdWScalingFactorSCF;
            valueVec = valueVec * (atomA.xyz - atomB.xyz)/distance;
        end
        % not tested
        function value = GetDiatomVdWCorrection2ndDerivative(obj, atomA, atomB, axisA1, axisA2)
            distance = obj.molecule.GetDistanceAtoms(atomA, atomB);
            dCartesian1 = atomA.xyz(axisA1) - atomB.xyz(axisA1);
            dCartesian2 = atomA.xyz(axisA2) - atomB.xyz(axisA2);
            vdWDistance = atomA.vdWRadii + atomB.vdWRadii;
            vdWScalingFacotor = Arguments.GetInstance().vdWScalingFactorSCF;
            tmpSum = atomA.vdWCoefficient+atomB.vdWCoefficient;
            if(tmpSum<=0e0)
                value = 0;
                return;
            end
            vdWCoefficients = 2.0*atomA.vdWCoefficient*atomB.vdWCoefficient/tmpSum;
%             dampingFactor = Arguments.GetInstance().vdWDampingFactorSCF;
            damping = obj.GetVdWDampingValue(vdWDistance, distance);
            damping1stDerivative = obj.GetVdWDampingValue1stDerivative(vdWDistance, distance);
            damping2ndDerivative = obj.GetVdWDampingValue2ndDerivative(vdWDistance, distance);
            
            dis6 = distance ^ 6; %distance*distance*distance*distance*distance*distance;
            tmp1 = -6.0*damping             /(dis6*distance) ...
                +    damping1stDerivative/dis6;
            tmp2 = 42.0*damping             /(dis6*distance*distance) ...
                -12.0*damping1stDerivative/(dis6*distance)...
                +     damping2ndDerivative/dis6;
            
            if(axisA1 ~= axisA2)
                pre1 = -dCartesian1*dCartesian2/(distance*distance*distance);
                pre2 =  dCartesian1*dCartesian2/(distance*distance);
            else
                pre1 = 1.0/distance - dCartesian1*dCartesian1/(distance*distance*distance);
                pre2 = (dCartesian1*dCartesian1)/(distance*distance);
            end
            
            value = pre1*tmp1 + pre2*tmp2;
            value = value * (-1.0*vdWScalingFacotor*vdWCoefficients);

        end
        
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
%                 I = 2*double(EnumShell.ShellType_end)-1;
%                 J = 2*double(EnumShell.ShellType_end)-1;
                I = 9;
                J = 9;
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
        
        % not tested
        function value = GetReducedOverlapAOs1stDerivativeAlpha(obj, na, la, m, nb, lb, alpha, beta)
            value = 0.0;
%             I = 2*ShellType_end-1;
%             J = 2*ShellType_end-1;
            I = 9;
            J = 9;
            
            for i = 0:I-1
                for j = 0:J-1
                    if(0e0<abs(obj.ReducedOverlapAOsParameters_Y(na+1,nb+1,la+1,lb+1,m+1,i+1,j+1)))
                        temp1 = obj.GetAuxiliaryA1stDerivative(i, 0.5*(alpha+beta))...
                            *obj.GetAuxiliaryB(j, 0.5*(alpha-beta));
                        temp2 = obj.GetAuxiliaryA(i, 0.5*(alpha+beta))...
                            *obj.GetAuxiliaryB1stDerivative(j, 0.5*(alpha-beta));
                        value = value + obj.ReducedOverlapAOsParameters_Y(na+1,nb+1,la+1,lb+1,m+1,i+1,j+1)*(temp1 + temp2);
                    end
                end
            end
            value = value * 0.5*obj.GetAuxiliaryD(la, lb, m);
        end
        % not tested
        function value = GetReducedOverlapAOs1stDerivativeBeta(obj, na, la, m, nb, lb, alpha, beta)
            value = 0.0;
%             I = 2*ShellType_end-1;
%             J = 2*ShellType_end-1;
            I = 9;
            J = 9;
            for i = 0:I-1
                for j = 0:J-1
                    if(0e0<abs(obj.ReducedOverlapAOsParameters_Y(na+1,nb+1,la+1,lb+1,m+1,i+1,j+1)))
                        temp1 = obj.GetAuxiliaryA1stDerivative(i, 0.5*(alpha+beta))...
                            *obj.GetAuxiliaryB(j, 0.5*(alpha-beta));
                        temp2 = obj.GetAuxiliaryA(i, 0.5*(alpha+beta))...
                            *obj.GetAuxiliaryB1stDerivative(j, 0.5*(alpha-beta));
                        value = value + obj.ReducedOverlapAOsParameters_Y(na+1,nb+1,la+1,lb+1,m+1,i+1,j+1)*(temp1 - temp2);
                    end
                end
            end
            value = value * 0.5*obj.GetAuxiliaryD(la, lb, m);
        end
        % not tested
        function value = GetReducedOverlapAOs2ndDerivativeAlpha(obj, na, la, m, nb, lb, alpha, beta)
            value = 0.0;
%             I = 2*ShellType_end-1;
%             J = 2*ShellType_end-1;
            I = 9;
            J = 9;
            for i = 0:I-1
                for j = 0:J-1
                    if(0e0<fabs(obj.ReducedOverlapAOsParameters_Y(na+1,nb+1,la+1,lb+1,m+1,i+1,j+1)))
                        temp1 = obj.GetAuxiliaryA2ndDerivative(i, 0.5*(alpha+beta))...
                            *obj.GetAuxiliaryB(j, 0.5*(alpha-beta));
                        temp2 = obj.GetAuxiliaryA(i, 0.5*(alpha+beta))...
                            *obj.GetAuxiliaryB2ndDerivative(j, 0.5*(alpha-beta));
                        temp3 = obj.GetAuxiliaryA1stDerivative(i, 0.5*(alpha+beta))...
                            *obj.GetAuxiliaryB1stDerivative(j, 0.5*(alpha-beta));
                        value = value + obj.ReducedOverlapAOsParameters_Y(na,nb,la+1,lb+1,m,i,j)*(temp1 + temp2 + 2.0*temp3);
                    end
                end
            end
            value = value * 0.25*obj.GetAuxiliaryD(la, lb, m);
        end
        % not tested
        function value = GetReducedOverlapAOs2ndDerivativeBeta(obj, na, la, m, nb, lb, alpha, beta)
            value = 0.0;
%             I = 2*ShellType_end-1;
%             J = 2*ShellType_end-1;
            I = 9;
            J = 9;
            for i = 0:I-1
                for j = 0:J-1
                    if(0e0<abs(obj.ReducedOverlapAOsParameters_Y(na+1,nb+1,la+1,lb+1,m+1,i+1,j+1)))
                        temp1 = obj.GetAuxiliaryA2ndDerivative(i, 0.5*(alpha+beta))...
                            *obj.GetAuxiliaryB(j, 0.5*(alpha-beta));
                        temp2 = obj.GetAuxiliaryA(i, 0.5*(alpha+beta))...
                            *obj.GetAuxiliaryB2ndDerivative(j, 0.5*(alpha-beta));
                        temp3 = obj.GetAuxiliaryA1stDerivative(i, 0.5*(alpha+beta))...
                            *obj.GetAuxiliaryB1stDerivative(j, 0.5*(alpha-beta));
                        value = value + obj.ReducedOverlapAOsParameters_Y(na+1,nb+1,la+1,lb+1,m+1,i+1,j+1)*(temp1 + temp2 - 2.0*temp3);
                    end
                end
            end
            value = value * 0.25*obj.GetAuxiliaryD(la, lb, m);
        end
        % not tested
        function value = GetReducedOverlapAOs2ndDerivativeAlphaBeta(obj, na, la, m, nb, lb, alpha, beta)
            value = 0.0;
%             I = 2*ShellType_end-1;
%             J = 2*ShellType_end-1;
            I = 9;
            J = 9;
            
            for i = 0:I-1
                for j = 0:J-1
                    if(0e0<abs(obj.ReducedOverlapAOsParameters.Y(na+1,nb+1,la+1,lb+1,m+1,i+1,j+1)))
                        temp1 = obj.GetAuxiliaryA2ndDerivative(i, 0.5*(alpha+beta))...
                            *obj.GetAuxiliaryB(j, 0.5*(alpha-beta));
                        temp2 = obj.GetAuxiliaryA(i, 0.5*(alpha+beta))...
                            *obj.GetAuxiliaryB2ndDerivative(j, 0.5*(alpha-beta));
                        value = value + obj.ReducedOverlapAOsParameters.Y(na+1,nb+1,la+1,lb+1,m+1,i+1,j+1)*(temp1 - temp2);
                    end
                end
            end
            value = value * 0.25*obj.GetAuxiliaryD(la, lb, m);
        end
        
        %    double GetOverlapAOsElement1stDerivativeByGTOExpansion(const MolDS_base_atoms::Atom& atomA,
        %                                                           int valenceIndexA,
        %                                                           const MolDS_base_atoms::Atom& atomB,
        %                                                           int valenceIndexB,
        %                                                           MolDS_base::STOnGType stonG,
        %                                                           MolDS_base::CartesianType axisA) const; // See [DY_1977].
        
        % see J. Mol. Struc. (Theochem), 419, 19 (1997) (ref. [BFB_1997])
        % we set gamma=0 always.
        function rotatingMatrix = CalcRotatingMatrix(~, atomA, atomB)
%             rotatingMatrix = zeros(double(EnumOrbital.OrbitalType_end));
            rotatingMatrix = zeros(9);
            eulerAngle = EulerAngle(atomB.xyz - atomA.xyz);
            alpha = eulerAngle.alpha;
            beta  = eulerAngle.beta;
            
            s = 1;
            py = 2;
            pz = 3;
            px = 4;
            dxy = 5;
            dyz = 6;
            dzz = 7;
            dzx = 8;
            dxxyy = 9;
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
%             dMatrix = zeros(double(EnumOrbital.OrbitalType_end));
            dMatrix = zeros(19);
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
        
        function gammaAB = CalcGammaAB(obj)
            totalAtomNumber = length(obj.molecule.atomVect);
            gammaAB = zeros(totalAtomNumber);
            atomvect = obj.molecule.atomVect;
            for A = 1:totalAtomNumber
                atomA = atomvect{A};
                na = double(atomA.valenceShellType);
%                 orbitalExponentA = atomA.GetOrbitalExponent(atomA.GetValenceShellType(), EnumOrbital.s, obj.theory);
                orbitalExponentA = obj.AtomGetOrbitalExponent(atomA, atomA.valenceShellType, 1);
                for B = A:totalAtomNumber
                    atomB = atomvect{B};
                    nb = double(atomB.valenceShellType);
%                     orbitalExponentB = atomB.GetOrbitalExponent(atomB.GetValenceShellType(), EnumOrbital.s, obj.theory);
                    orbitalExponentB = obj.AtomGetOrbitalExponent(atomB, atomB.valenceShellType, 1);
                    
                    R = obj.molecule.distanceAtoms(A, B);
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
%                             temp = temp * factorial(2*na+2*nb-l-1);
%                             temp = temp / factorial(2*nb-l);
                            temp = temp * prod(2*nb-l+1:2*na+2*nb-l-1);
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
        
        function value = GetFockDiagElement(obj, atomA, mu, isGuess)
            indexAtomA = atomA.index;
            firstAOIndexA = atomA.firstAOIndex;
            value = obj.AtomGetCoreIntegral(atomA, atomA.valence(mu-firstAOIndexA+1), ...
                obj.gammaAB(indexAtomA, indexAtomA), ...
                isGuess);
            if(~isGuess)
                atomvect = obj.molecule.atomVect;
                temp = obj.atomicElectronPopulation(indexAtomA) ...
                    -0.5*diag(obj.orbitalElectronPopulation(mu, mu));
                value = value + temp*obj.gammaAB(indexAtomA, indexAtomA);
                
                temp = 0.0;
                for BB = 1:length(atomvect)
                    if(BB ~= indexAtomA)
                        atomBB = atomvect{BB};
                        temp = temp + ( obj.atomicElectronPopulation(BB) - atomBB.coreCharge  )...
                            *obj.gammaAB(indexAtomA, BB);
                    end
                end
                value = value + temp;
            end
        end
        
        function value = GetFockOffDiagElement(obj, atomA, atomB, mu, nu, isGuess)
            K = obj.GetBondingAdjustParameterK(atomA.valenceShellType, atomB.valenceShellType);
            bondParameter = 0.5*K*(obj.AtomGetBondingParameter(atomA) + obj.AtomGetBondingParameter(atomB));
            value =  bondParameter*obj.overlapAOs(mu, nu);
            if(~isGuess)
                value = value - 0.5*obj.orbitalElectronPopulation(mu, nu)*obj.gammaAB(atomA.index, atomB.index);
            end
        end
        
        %    void TransposeFockMatrixMatrix(double** transposedFockMatrix) const;
        
        function diatomicOverlapAOs = CalcDiatomicOverlapAOsInDiatomicFrame(obj, atomA, atomB)
            na = double(atomA.valenceShellType);
            nb = double(atomB.valenceShellType);
            
%             diatomicOverlapAOs = zeros(double(EnumOrbital.OrbitalType_end));
            diatomicOverlapAOs = zeros(9);
            rAB = obj.molecule.GetDistanceAtoms(atomA, atomB); % Inter nuclear distance between aton A and B.
            
            for a = 1:length(atomA.valence)
                valenceOrbitalA = atomA.valence(a);
                realShpericalHarmonicsA = atomA.realSphericalHarmonicsIndices{a};
                orbitalExponentA = obj.AtomGetOrbitalExponent(atomA, atomA.valenceShellType, valenceOrbitalA);

                for b = 1:length(atomB.valence)
                    valenceOrbitalB = atomB.valence(b);
                    realShpericalHarmonicsB = atomB.realSphericalHarmonicsIndices{b};
                    orbitalExponentB = obj.AtomGetOrbitalExponent(atomB, atomB.valenceShellType, valenceOrbitalB);

                    if(realShpericalHarmonicsA.m == realShpericalHarmonicsB.m)
                        m = abs(realShpericalHarmonicsA.m);
                        alpha = orbitalExponentA * rAB;
                        beta =  orbitalExponentB * rAB;

                        reducedOverlapAOs = obj.GetReducedOverlapAOs(na, realShpericalHarmonicsA.l, m,...
                            nb, realShpericalHarmonicsB.l, alpha, beta);


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
        
        % not tested
        function diatomicOverlapAOsDeri = CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(obj, atomA, atomB)
            na = double(atomA.valenceShellType);
            nb = double(atomB.valenceShellType);
            
            diatomicOverlapAOsDeri = zeros(9);
            
            R = obj.molecule.GetDistanceAtoms(atomA, atomB);
            
            for a = 1:length(atomA.valence)
                valenceOrbitalA = atomA.valence(a);
                realShpericalHarmonicsA = RealSphericalHarmonicsIndex(valenceOrbitalA);
                orbitalExponentA = obj.AtomGetOrbitalExponent(atomA, ...
                    atomA.valenceShellType, ...
                    valenceOrbitalA);
                
                for b = 1:length(atomB.valence)
                    valenceOrbitalB = atomB.valence(b);
                    realShpericalHarmonicsB = RealSphericalHarmonicsIndex(valenceOrbitalB);
                    orbitalExponentB = obj.AtomGetOrbitalExponent(atomB. ...
                        atomB.valenceShellType, ...
                        valenceOrbitalB);
                    
                    if(realShpericalHarmonicsA.m == realShpericalHarmonicsB.m)
                        m = abs(realShpericalHarmonicsA.m);
                        alpha = orbitalExponentA * R;
                        beta =  orbitalExponentB * R;
                        
                        reducedOverlapAOs = obj.GetReducedOverlapAOs(na, realShpericalHarmonicsA.l, m,...
                            nb, realShpericalHarmonicsB.l, alpha, beta);
                        reducedOverlapAOs1stDerivAlpha = obj.GetReducedOverlapAOs1stDerivativeAlpha(...
                            na, ...
                            realShpericalHarmonicsA.l, ...
                            m,...
                            nb, ...
                            realShpericalHarmonicsB.l, ...
                            alpha, ...
                            beta);
                        reducedOverlapAOs1stDerivBeta  = obj.GetReducedOverlapAOs1stDerivativeBeta(...
                            na, ...
                            realShpericalHarmonicsA.l, ...
                            m,...
                            nb, ...
                            realShpericalHarmonicsB.l, ...
                            alpha, ...
                            beta);
                        
                        temp1 = (na+nb+1)*power(R,na+nb)*reducedOverlapAOs;
                        temp2 = power(R,na+nb+1)*(orbitalExponentA*reducedOverlapAOs1stDerivAlpha...
                            +orbitalExponentB*reducedOverlapAOs1stDerivBeta);
                        
                        pre =  power(2.0*orbitalExponentA, na+0.5);
                        pre = pre * power(2.0*orbitalExponentB, nb+0.5);
                        factorials = factorial(2*na)*factorial(2*nb);
                        pre = pre / sqrt(factorials);
                        pre = pre / power(2.0, na+nb+1.0);
                        
                        diatomicOverlapAOsDeri(valenceOrbitalA,valenceOrbitalB) = pre*(temp1+temp2);
                    end
                    
                end
            end
        end
        % not tested
        function diatomicOverlapAOs2ndDeri = CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(obj, atomA, atomB)
            na = double(atomA.valenceShellType);
            nb = double(atomB.valenceShellType);
            
            diatomicOverlapAOs2ndDeri = zeros(9);
            R = obj.molecule.GetDistanceAtoms(atomA, atomB);
            
            for a = 1:length(atomA.valence)
                valenceOrbitalA = atomA.valence(a);
                realShpericalHarmonicsA = RealSphericalHarmonicsIndex(valenceOrbitalA);
                orbitalExponentA = obj.AtomGetOrbitalExponent(atomA, atomA.GetValenceShellType(),...
                    valenceOrbitalA);
                
                for b = 1:length(atomB.valence)
                    valenceOrbitalB = atomB.valence(b);
                    realShpericalHarmonicsB = RealSphericalHarmonicsIndex(valenceOrbitalB);
                    orbitalExponentB = obj.AtomGetOrbitalExponent(atomB, atomB.GetValenceShellType(),...
                        valenceOrbitalB);
                    
                    if(realShpericalHarmonicsA.m == realShpericalHarmonicsB.m)
                        m = abs(realShpericalHarmonicsA.m);
                        alpha = orbitalExponentA * R;
                        beta =  orbitalExponentB * R;
                        
                        reducedOverlapAOs = obj.GetReducedOverlapAOs(na,...
                            realShpericalHarmonicsA.l,...
                            m,...
                            nb,...
                            realShpericalHarmonicsB.l,...
                            alpha,...
                            beta);
                        reducedOverlapAOs1stDerivAlpha...
                            = obj.GetReducedOverlapAOs1stDerivativeAlpha(na,...
                            realShpericalHarmonicsA.l,...
                            m,...
                            nb,...
                            realShpericalHarmonicsB.l,...
                            alpha,...
                            beta);
                        reducedOverlapAOs1stDerivBeta...
                            = obj.GetReducedOverlapAOs1stDerivativeBeta(na,...
                            realShpericalHarmonicsA.l,...
                            m,...
                            nb,...
                            realShpericalHarmonicsB.l,...
                            alpha,...
                            beta);
                        reducedOverlapAOs2ndDerivAlpha...
                            = obj.GetReducedOverlapAOs2ndDerivativeAlpha(na,...
                            realShpericalHarmonicsA.l,...
                            m,...
                            nb,...
                            realShpericalHarmonicsB.l,...
                            alpha,...
                            beta);
                        reducedOverlapAOs2ndDerivBeta...
                            = obj.GetReducedOverlapAOs2ndDerivativeBeta(na,...
                            realShpericalHarmonicsA.l,...
                            m,...
                            nb,...
                            realShpericalHarmonicsB.l,...
                            alpha,...
                            beta);
                        reducedOverlapAOs2ndDerivAlphaBeta...
                            = obj.GetReducedOverlapAOs2ndDerivativeAlphaBeta(na,...
                            realShpericalHarmonicsA.l,...
                            m,...
                            nb,...
                            realShpericalHarmonicsB.l,...
                            alpha,...
                            beta);
                        
                        temp1 = (na+nb+1)...
                            *(na+nb)...
                            *power(R,na+nb-1)*reducedOverlapAOs;
                        temp2 = 2.0*(na+nb+1)*power(R,na+nb)...
                            *(orbitalExponentA*reducedOverlapAOs1stDerivAlpha...
                            +orbitalExponentB*reducedOverlapAOs1stDerivBeta);
                        temp3 = power(R,na+nb+1)...
                            *(power(orbitalExponentA,2.0)*reducedOverlapAOs2ndDerivAlpha...
                            +power(orbitalExponentB,2.0)*reducedOverlapAOs2ndDerivBeta...
                            +2.0*orbitalExponentA*orbitalExponentB*reducedOverlapAOs2ndDerivAlphaBeta);
                        
                        pre =  power(2.0*orbitalExponentA, na+0.5);
                        pre = pre * power(2.0*orbitalExponentB, nb+0.5);
                        factorials = factorial(2*na)*factorial(2*nb);
                        pre = pre / sqrt(factorials);
                        pre = pre / power(2.0, na+nb+1.0);
                        
                        diatomicOverlapAOs2ndDeri(valenceOrbitalA,valenceOrbitalB) = pre*(temp1+temp2+temp3);
                    end
                    
                end
            end
        end
        % not tested
        function diatomicOverlapAOs1stDerivs = CalcDiatomicOverlapAOs1stDerivatives(obj, atomA, atomB)
            if(isnumeric(atomA) || isnumeric(atomB))
                atomA = obj.molecule.atomVect{atomA};
                atomB = obj.molecule.atomVect{atomB};
            end
            
            diatomicOverlapAOs1stDerivs = zeros(9,9,3);
            cartesian = atomA.xyz - atomB.xyz;
            
            R = obj.molecule.GetDistanceAtoms(atomA, atomB);
            
            tmpDiaOverlapAOsInDiaFrame = obj.CalcDiatomicOverlapAOsInDiatomicFrame(atomA, atomB);
            tmpDiaOverlapAOs1stDerivInDiaFrame = obj.CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(atomA, atomB);
            tmpRotMat = obj.CalcRotatingMatrix(atomA, atomB);
            tmpRotMat1stDerivs  = obj.CalcRotatingMatrix1stDerivatives(atomA, atomB);
            
            % rotate (fast algorithm, see also slow algorithm shown later)
            for c = 1:3
                tmpRotMat1stDeriv = tmpRotMat1stDerivs(:,:,c);
                
                alpha = cartesian(c)/R;
                
                tmpRotatedDiatomicOverlap = alpha * tmpRotMat * tmpDiaOverlapAOs1stDerivInDiaFrame * tmpRotMat;
                
                tmpRotatedDiatomicOverlap = tmpRotMat1stDeriv * tmpDiaOverlapAOsInDiaFrame * tmpRotMat + tmpRotatedDiatomicOverlap;
                
                tmpRotatedDiatomicOverlap = tmpRotMat * tmpDiaOverlapAOsInDiaFrame * tmpRotMat1stDeriv + tmpRotatedDiatomicOverlap;
                
                diatomicOverlapAOs1stDerivs(:,:,c) = tmpRotatedDiatomicOverlap;
            end
        end
        % not tested
        function diatomicOverlapAOs2ndDerivs = CalcDiatomicOverlapAOs2ndDerivatives(obj, atomA, atomB)
            if(isnumeric(atomA) || isnumeric(atomB))
                atomA = obj.molecule.atomVect{atomA};
                atomB = obj.molecule.atomVect{atomB};
            end
            
            diatomicOverlapAOs2ndDerivs = zeros(9,9,3,3);
            cartesian = atomA.xyz - atomB.xyz;
            R = obj.molecule.GetDistanceAtoms(atomA, atomB);
            
            tmpDiaOverlapAOsInDiaFrame = obj.CalcDiatomicOverlapAOsInDiatomicFrame(atomA, atomB);
            tmpDiaOverlapAOs1stDerivInDiaFrame = obj.CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(atomA, atomB);
            tmpDiaOverlapAOs2ndDerivInDiaFrame = obj.CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(atomA, atomB);
            tmpRotMat = obj.CalcRotatingMatrix(atomA, atomB);
            tmpRotMat1stDerivs = obj.CalcRotatingMatrix1stDerivatives(atomA, atomB);
            tmpRotMat2ndDerivs = obj.CalcRotatingMatrix2ndDerivatives(atomA, atomB);
            
            % calculate each element of first derivatives
            tmpDiaOverlapAOs1stDerivs = zeros(9,9,3);
            for i = 1:9
                for j = 1:9
                    for dimA1 = 1:3
                        tmpDiaOverlapAOs1stDerivs(i,j,dimA1) = (cartesian(dimA1)/R)*tmpDiaOverlapAOs1stDerivInDiaFrame(i,j);
                    end
                end
            end
            
            % calculate each element of second derivatives
            tmpDiaOverlapAOs2ndDerivs = zeros(9,9,3,3);
            for i = 1:9
                for j = 1:9
                    for dimA1 = 1:3
                        for dimA2 = 1:3
                            tmpDiaOverlapAOs2ndDerivs(i,j,dimA1,dimA2) ...
                                = obj.Get2ndDerivativeElementFromDistanceDerivatives(tmpDiaOverlapAOs1stDerivInDiaFrame(i,j),...
                                tmpDiaOverlapAOs2ndDerivInDiaFrame(i,j),...
                                static_cast<CartesianType>(dimA1),...
                                static_cast<CartesianType>(dimA2),...
                                cartesian,...
                                R);...
                        end
                    end
                end
            end
            
            % rotate
            for i = 1:9
                for j = 1:9
                    for dimA1 = 1:3
                        for dimA2 = 1:3
                            diatomicOverlapAOs2ndDerivs(i,j,dimA1,dimA2) = 0.0;
                            
                            temp1 = 0.0; temp2=0.0; temp3 = 0.0;
                            temp4 = 0.0; temp5=0.0; temp6 = 0.0;
                            temp7 = 0.0; temp8=0.0; temp9 = 0.0;
                            for k = 1:9
                                for l = 1:9
                                    
                                    temp1 = temp1 + tmpRotMat2ndDerivs        (i,k,dimA1,dimA2)...
                                        *tmpRotMat                 (j,l)...
                                        *tmpDiaOverlapAOsInDiaFrame(k,l);
                                    temp2 = temp2 + tmpRotMat                 (i,k)...
                                        *tmpRotMat2ndDerivs        (j,l,dimA1,dimA2)...
                                        *tmpDiaOverlapAOsInDiaFrame(k,l);
                                    temp3 = temp3 + tmpRotMat                 (i,k)...
                                        *tmpRotMat                 (j,l)...
                                        *tmpDiaOverlapAOs2ndDerivs (k,l,dimA1,dimA2);
                                    temp4 = temp4 + tmpRotMat1stDerivs        (i,k,dimA1) ...
                                        *tmpRotMat1stDerivs        (j,l,dimA2)...
                                        *tmpDiaOverlapAOsInDiaFrame(k,l);
                                    temp5 = temp5 + tmpRotMat1stDerivs        (i,k,dimA1) ...
                                        *tmpRotMat                 (j,l)...
                                        *tmpDiaOverlapAOs1stDerivs (k,l,dimA2);
                                    temp6 = temp6 + tmpRotMat1stDerivs        (i,k,dimA2) ...
                                        *tmpRotMat1stDerivs        (j,l,dimA1)...
                                        *tmpDiaOverlapAOsInDiaFrame(k,l);
                                    temp7 = temp7 + tmpRotMat                 (i,k) ...
                                        *tmpRotMat1stDerivs        (j,l,dimA1)...
                                        *tmpDiaOverlapAOs1stDerivs (k,l,dimA2);
                                    temp8 = temp8 + tmpRotMat1stDerivs        (i,k,dimA2) ...
                                        *tmpRotMat                 (j,l)...
                                        *tmpDiaOverlapAOs1stDerivs (k,l,dimA1);
                                    temp9 = temp9 + tmpRotMat                 (i,k) ...
                                        *tmpRotMat1stDerivs        (j,l,dimA2)...
                                        *tmpDiaOverlapAOs1stDerivs (k,l,dimA1);
                                end
                            end
                            
                            diatomicOverlapAOs2ndDerivs(i,j,dimA1,dimA2) = temp1+temp2+temp3...
                                +temp4+temp5+temp6...
                                +temp7+temp8+temp9;
                        end
                    end
                end
            end
        end
        % not tested
        function value = Get2ndDerivativeElementFromDistanceDerivatives(~, firstDistanceDeri,...
                secondDistanceDeri,...
                axisA1,...
                axisA2,...
                cartesian,...
                rAB)              
            if(axisA1 ~= axisA2)
                value = -1.0*firstDistanceDeri/(rAB*rAB*rAB);
                value = value + secondDistanceDeri/(rAB*rAB);
                value = value * cartesian(axisA1)*cartesian(axisA2);
            else
                value = (rAB*rAB - cartesian(axisA1)*cartesian(axisA1))*firstDistanceDeri/(rAB*rAB*rAB);
                value = value + cartesian(axisA1)*cartesian(axisA1)*secondDistanceDeri/(rAB*rAB);
            end
        end
        
        %    virtual double GetMolecularIntegralElement(int moI, int moJ, int moK, int moL,
        %                                               const MolDS_base::Molecule& molecule,
        %                                               double const* const* fockMatrix,
        %                                               double const* const* gammaAB) const;
        
%         function [twoElecsTwoAtomCores, twoElecsAtomEpcCores] = CalcTwoElecsTwoCores(~, ~)
%             % do nothing for CNDO, INDO, and ZINDO/S.
%             % two electron two core integrals are not needed for CNDO, INDO, and ZINDO/S.
%             twoElecsTwoAtomCores = [];
%             twoElecsAtomEpcCores = [];
%         end
        
        %    virtual void CalcForce(const std::vector<int>& elecStates);
        
        % not tested
        function rotMat1stDerivatives = CalcRotatingMatrix1stDerivatives(~, atomA, atomB)
            rotMat1stDerivatives = zeros(9, 9, 3);
            x = atomB.xyz(1) - atomA.xyz(1);
            y = atomB.xyz(2) - atomA.xyz(2);
            z = atomB.xyz(3) - atomA.xyz(3);
            r = sqrt( x*x + y*y );
            R = sqrt( x*x + y*y + z*z );
            
            if(r==0e0)
                return;
            end
            
            % for s-function
            rotMat1stDerivatives(s,s,1) = 0.0;
            rotMat1stDerivatives(s,s,2) = 0.0;
            rotMat1stDerivatives(s,s,3) = 0.0;
            
            % for p-function
            rotMat1stDerivatives(py,py,1) = -1.0/r + x*x/(r*r*r);
            rotMat1stDerivatives(py,pz,1) = x*y/(R*R*R);
            rotMat1stDerivatives(py,px,1) = (1.0/(r*r*r*R) + 1.0/(R*R*R*r))*x*y*z;
            
            rotMat1stDerivatives(pz,py,1) = 0.0;
            rotMat1stDerivatives(pz,pz,1) = x*z/(R*R*R);
            rotMat1stDerivatives(pz,px,1) = x/(r*R) - x*r/(R*R*R);
            
            rotMat1stDerivatives(px,py,1) = -1.0*x*y/(r*r*r);
            rotMat1stDerivatives(px,pz,1) = -1.0/R + x*x/(R*R*R);
            rotMat1stDerivatives(px,px,1) = -1.0*z/(r*R) + ...
                (1.0/(r*r*r*R) + 1.0/(R*R*R*r))*x*x*z;
            
            rotMat1stDerivatives(py,py,2) = x*y/(r*r*r);
            rotMat1stDerivatives(py,pz,2) = -1.0/R + y*y/(R*R*R);
            rotMat1stDerivatives(py,px,2) = -1.0*z/(r*R) +...
                (1.0/(r*r*r*R) + 1.0/(R*R*R*r))*y*y*z;
            
            rotMat1stDerivatives(pz,py,2) = 0.0;
            rotMat1stDerivatives(pz,pz,2) = y*z/(R*R*R);
            rotMat1stDerivatives(pz,px,2) = y/(r*R) - y*r/(R*R*R);
            
            rotMat1stDerivatives(px,py,2) = 1.0/r - y*y/(r*r*r);
            rotMat1stDerivatives(px,pz,2) = x*y/(R*R*R);
            rotMat1stDerivatives(px,px,2) = (1.0/(r*r*r*R) + 1.0/(R*R*R*r))*x*y*z;
            
            rotMat1stDerivatives(py,py,3) = 0.0;
            rotMat1stDerivatives(py,pz,3) = y*z/(R*R*R);
            rotMat1stDerivatives(py,px,3) = -1.0*y/(r*R) + y*z*z/(r*R*R*R);
            
            rotMat1stDerivatives(pz,py,3) = 0.0;
            rotMat1stDerivatives(pz,pz,3) = -1.0/R + z*z/(R*R*R);
            rotMat1stDerivatives(pz,px,3) = -1.0*z*r/(R*R*R);
            
            rotMat1stDerivatives(px,py,3) = 0.0;
            rotMat1stDerivatives(px,pz,3) = x*z/(R*R*R);
            rotMat1stDerivatives(px,px,3) = -1.0*x/(r*R) + x*z*z/(r*R*R*R);
            
            % for d-function
            % ToDo: First derivative of rotating matrix for d-orbital...
        end
        % not tested
        function rotMat2ndDerivatives = CalcRotatingMatrix2ndDerivatives(~, atomA, atomB)
            rotMat2ndDerivatives = zeros(9,9,3,3);
            x = atomB.xyz(1) - atomA.xyz(1);
            y = atomB.xyz(2) - atomA.xyz(2);
            z = atomB.xyz(3) - atomA.xyz(3);
            r = sqrt( x*x + y*y );
            R = sqrt( x*x + y*y + z*z );
            
            if(r==0e0)
                return;
            end
            
            temp1 = 1.0/(r*r*r*R) + 1.0/(r*R*R*R);
            temp2 = 2.0/(r*r*r*R*R*R) + 3.0/(r*r*r*r*r*R) + 3.0/(r*R*R*R*R*R);
            temp3 = 1.0/(r*r*r*R*R*R) + 3.0/(r*R*R*R*R*R);
            
            % for s-function
            rotMat2ndDerivatives(s,s,1,1) = 0.0;
            rotMat2ndDerivatives(s,s,1,2) = 0.0;
            rotMat2ndDerivatives(s,s,1,3) = 0.0;
            rotMat2ndDerivatives(s,s,2,1) = 0.0;
            rotMat2ndDerivatives(s,s,2,2) = 0.0;
            rotMat2ndDerivatives(s,s,2,3) = 0.0;
            rotMat2ndDerivatives(s,s,3,1) = 0.0;
            rotMat2ndDerivatives(s,s,3,2) = 0.0;
            rotMat2ndDerivatives(s,s,3,3) = 0.0;
            
            % for p-function, xx-derivatives
            rotMat2ndDerivatives(py,py,1,1) = -3.0*x/(r*r*r) + 3.0*x*x*x/(r*r*r*r*r);
            rotMat2ndDerivatives(py,pz,1,1) = -1.0*y/(R*R*R) + 3.0*x*x*y/(R*R*R*R*R);
            rotMat2ndDerivatives(py,px,1,1) = -1.0*temp1*y*z+temp2*x*x*y*z;
            
            rotMat2ndDerivatives(pz,py,1,1) = 0.0;
            rotMat2ndDerivatives(pz,pz,1,1) = -1.0*z/(R*R*R) + 3.0*x*x*z/(R*R*R*R*R);
            rotMat2ndDerivatives(pz,px,1,1) = -1.0/(r*R) + temp1*x*x...
                +r/(R*R*R) - 3.0*x*x*r/(R*R*R*R*R) + x*x/(r*R*R*R);
            
            rotMat2ndDerivatives(px,py,1,1) = y/(r*r*r) - 3.0*x*x*y/(r*r*r*r*r);
            rotMat2ndDerivatives(px,pz,1,1) = -3.0*x/(R*R*R) + 3.0*x*x*x/(R*R*R*R*R);
            rotMat2ndDerivatives(px,px,1,1) = -3.0*temp1*x*z+temp2*x*x*x*z;
            
            % for p-function, xy-derivatives
            rotMat2ndDerivatives(py,py,1,2) = -1.0*y/(r*r*r) + 3.0*x*x*y/(r*r*r*r*r);
            rotMat2ndDerivatives(py,pz,1,2) = -1.0*x/(R*R*R) + 3.0*x*y*y/(R*R*R*R*R);
            rotMat2ndDerivatives(py,px,1,2) = -1.0*temp1*x*z+temp2*x*y*y*z;
            
            rotMat2ndDerivatives(pz,py,1,2) = 0.0;
            rotMat2ndDerivatives(pz,pz,1,2) = 3.0*x*y*z/(R*R*R*R*R);
            rotMat2ndDerivatives(pz,px,1,2) = temp1*x*y + x*y/(r*R*R*R) - 3.0*x*y*r/(R*R*R*R*R);
            
            rotMat2ndDerivatives(px,py,1,2) = x/(r*r*r) - 3.0*x*y*y/(r*r*r*r*r);
            rotMat2ndDerivatives(px,pz,1,2) = rotMat2ndDerivatives(py,pz,1,1);
            rotMat2ndDerivatives(px,px,1,2) = rotMat2ndDerivatives(py,px,1,1);
            
            % for p-function, yx-derivatives
            for i = 2:4
                for j = 2:4
                    rotMat2ndDerivatives(i,j,2,1) = rotMat2ndDerivatives(i,j,1,2);
                end
            end
            % for p-function, xz-derivatives
            rotMat2ndDerivatives(py,py,1,3) = 0.0;
            rotMat2ndDerivatives(py,pz,1,3) = rotMat2ndDerivatives(pz,pz,1,2);
            rotMat2ndDerivatives(py,px,1,3) = -1.0*temp1*x*y +temp3*x*y*z*z;
            
            rotMat2ndDerivatives(pz,py,1,3) = 0.0;
            rotMat2ndDerivatives(pz,pz,1,3) = -1.0*x/(R*R*R) + 3.0*x*z*z/(R*R*R*R*R);
            rotMat2ndDerivatives(pz,px,1,3) = x*z/(r*R*R*R) - 3.0*x*z*r/(R*R*R*R*R);
            
            rotMat2ndDerivatives(px,py,1,3) = 0.0;
            rotMat2ndDerivatives(px,pz,1,3) = rotMat2ndDerivatives(pz,pz,1,1);
            rotMat2ndDerivatives(px,px,1,3) = 1.0/(r*R) - z*z/(r*R*R*R)...
                -1.0*temp1*x*x+temp3*x*x*z*z;
            
            
            % for p-function, zx-derivatives
            for i = 2:4
                for j = 2:4
                    rotMat2ndDerivatives(i,j,3,1) = rotMat2ndDerivatives(i,j,1,3);
                end
            end
            
            % for p-function, yy-derivatives
            rotMat2ndDerivatives(py,py,2,2) = -1.0*x/(r*r*r) + 3.0*x*y*y/(r*r*r*r*r);
            rotMat2ndDerivatives(py,pz,2,2) = -3.0*y/(R*R*R) + 3.0*y*y*y/(R*R*R*R*R);
            rotMat2ndDerivatives(py,px,2,2) = -3.0*temp1*y*z+temp2*y*y*y*z;
            
            rotMat2ndDerivatives(pz,py,2,2) = 0.0;
            rotMat2ndDerivatives(pz,pz,2,2) = -1.0*z/(R*R*R) + 3.0*y*y*z/(R*R*R*R*R);
            rotMat2ndDerivatives(pz,px,2,2) = -1.0/(r*R) + temp1*y*y...
            +r/(R*R*R) - 3.0*y*y*r/(R*R*R*R*R) + y*y/(r*R*R*R);
            
            rotMat2ndDerivatives(px,py,2,2) = 3.0*y/(r*r*r) - 3.0*y*y*y/(r*r*r*r*r);
            rotMat2ndDerivatives(px,pz,2,2) = rotMat2ndDerivatives(py,pz,1,2);
            rotMat2ndDerivatives(px,px,2,2) = -1.0*temp1*x*z+temp2*x*y*y*z;
            
            % for p-function, yz-derivatives
            rotMat2ndDerivatives(py,py,2,3) = 0.0;
            rotMat2ndDerivatives(py,pz,2,3) = rotMat2ndDerivatives(pz,pz,2,2);
            rotMat2ndDerivatives(py,px,2,3) = 1.0/(r*R) - z*z/(r*R*R*R)...
            -1.0*temp1*y*y+temp3*y*y*z*z;
            
            rotMat2ndDerivatives(pz,py,2,3) = 0.0;
            rotMat2ndDerivatives(pz,pz,2,3) = -1.0*y/(R*R*R) + 3.0*y*z*z/(R*R*R*R*R);
            rotMat2ndDerivatives(pz,px,2,3) = y*z/(r*R*R*R) - 3.0*y*z*r/(R*R*R*R*R);
            
            rotMat2ndDerivatives(px,py,2,3) = 0.0;
            rotMat2ndDerivatives(px,pz,2,3) = rotMat2ndDerivatives(pz,pz,1,2);
            rotMat2ndDerivatives(px,px,2,3) = -1.0*temp1*x*y+temp3*x*y*z*z;
            
            
            % for p-function, zy-derivatives
            for i = 2:4
                for j = 2:4
                    rotMat2ndDerivatives(i,j,3,2) = rotMat2ndDerivatives(i,j,2,3);
                end
            end
            
            % for p-function, zz-derivatives
            rotMat2ndDerivatives(py,py,3,3) = 0.0;
            rotMat2ndDerivatives(py,pz,3,3) = rotMat2ndDerivatives(pz,pz,2,3);
            rotMat2ndDerivatives(py,px,3,3) = -3.0*y*z/(r*R*R*R) + 3.0*y*z*z*z/(r*R*R*R*R*R);
            
            rotMat2ndDerivatives(pz,py,3,3) = 0.0;
            rotMat2ndDerivatives(pz,pz,3,3) = -3.0*z/(R*R*R) + 3.0*z*z*z/(R*R*R*R*R);
            rotMat2ndDerivatives(pz,px,3,3) = -3.0*z*z*r/(R*R*R*R*R) + r/(R*R*R);
            
            rotMat2ndDerivatives(px,py,3,3) = 0.0;
            rotMat2ndDerivatives(px,pz,3,3) = rotMat2ndDerivatives(pz,pz,1,3);
            rotMat2ndDerivatives(px,px,3,3) = -3.0*x*z/(r*R*R*R) + 3.0*x*z*z*z/(r*R*R*R*R*R);
            
            % for d-function
            % ToDo: Second derivative of rotating matrix for d-orbital...
        end
        
    end
    
    methods (Access = private)
        
        function CalcCoreRepulsionEnergy(obj)
            energy = 0.0;
            atomvect = obj.molecule.atomVect;
            for i = 1:length(atomvect)
                for j = i+1:length(atomvect)
                    energy = energy + obj.GetDiatomCoreRepulsionEnergy(atomvect{i}, atomvect{j});
                end
            end
            obj.coreRepulsionEnergy = energy;
        end
        
        % not tested
        function SetEnableAtomTypesVdW(obj)
            obj.enableAtomTypesVdW = {};
            obj.enableAtomTypesVdW{end+1} = EnumAtom.H;
            obj.enableAtomTypesVdW{end+1} = EnumAtom.C;
            obj.enableAtomTypesVdW{end+1} = EnumAtom.N;
            obj.enableAtomTypesVdW{end+1} = EnumAtom.O;
            obj.enableAtomTypesVdW{end+1} = EnumAtom.F;
            obj.enableAtomTypesVdW{end+1} = EnumAtom.S;
            obj.enableAtomTypesVdW{end+1} = EnumAtom.Cl;
        end
        % not tested
        function CheckEnableAtomTypeVdW(obj)
            if(Arguments.GetInstance().requiresVdWSCF)
                if(obj.theory == EnumTheory.AM1D || obj.theory == EnumTheory.PM3D)
                    return;
                end
            else
                return;
            end
            for i = 1:length(obj.molecule.atomVect)
                atomvect = obj.molecule.atomVect;
                atomType = atomvect{i}.atomType;
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
        % not tested
        function CalcVdWCorrectionEnergy(obj)
            value = 0.0;
            atomvect = obj.molecule.atomVect;
            for i = 1:length(atomvect)
                atomI = atomvect{i};
                for j = i+1:length(atomvect)
                    atomJ = atomvect{j};
                    value = value + obj.GetDiatomVdWCorrectionEnergy(atomI, atomJ);
                end
            end
            obj.vdWCorrectionEnergy = value;
        end
        % not tested
        function res = GetVdWDampingValue(~, vdWDistance, distance)
            dampingFactor = Arguments.GetInstance().vdWDampingFactorSCF;
            res = 1.0/(1.0+exp(-1.0*dampingFactor*(distance/vdWDistance - 1.0)));
        end
        % not tested
        function res = GetVdWDampingValue1stDerivative(~, vdWDistance, distance)
            dampingFactor = Arguments.GetInstance().vdWDampingFactorSCF;
            res = (dampingFactor/vdWDistance)...
                *exp(-1.0*dampingFactor*(distance/vdWDistance - 1.0))...
                /(1.0+exp(-1.0*dampingFactor*(distance/vdWDistance - 1.0)))...
                /(1.0+exp(-1.0*dampingFactor*(distance/vdWDistance - 1.0)));
        end
        % not tested
        function res = GetVdWDampingValue2ndDerivative(~, vdWDistance, distance)
            dampingFactor = Arguments.GetInstance().vdWDampingFactorSCF;
            exponent = -1.0*dampingFactor*(distance/vdWDistance - 1.0);
            pre = dampingFactor/vdWDistance;
            dominator = 1.0+exp(exponent);
            res = 2.0*pre*pre*exp(2.0*exponent)/(dominator*dominator*dominator) ...
                -    pre*pre*exp(    exponent)/(dominator*dominator);
        end
        
        function electronicTransitionDipoleMoments = CalcElectronicDipoleMomentGroundState(obj)
            groundState = 1;
            electronicTransitionDipoleMoments(groundState, groundState, :) = obj.CalcElectronicTransitionDipoleMoment(...
                groundState,...
                groundState);
        end
        
        function satisfy = SatisfyConvergenceCriterion(~, oldOrbitalElectronPopulation, ...
                orbitalElectronPopulation)
            satisfy = false;
            change = sum(sum((oldOrbitalElectronPopulation - orbitalElectronPopulation).^2));
            
            change = change / numel(oldOrbitalElectronPopulation);
            rmsDensity = sqrt(change);
            
            if(rmsDensity < Arguments.GetInstance().thresholdSCF)
                satisfy = true;
            end
        end
        
        function atomicElectronPopulation = CalcAtomicElectronPopulation(obj)
            atomvect = obj.molecule.atomVect;
            totalNumberAtoms = length(atomvect);
            atomicElectronPopulation = zeros(totalNumberAtoms, 1);
            for A = 1:totalNumberAtoms
                firstAOIndex = atomvect{A}.firstAOIndex;
                numberAOs    = length(atomvect{A}.valence);
                i = firstAOIndex:firstAOIndex+numberAOs-1;
                atomicElectronPopulation(A) = atomicElectronPopulation(A) + sum(diag(obj.orbitalElectronPopulation(i, i)));
            end
        end
        
        function coreDipoleMoment = CalcCoreDipoleMoment(obj)
             dipoleCenter = obj.molecule.GetXyzDipoleCenter();
             atomvect = obj.molecule.atomVect;
             coreDipoleMoment = zeros(3, 1);
             for A = 1:length(atomvect)
                 coreDipoleMoment = coreDipoleMoment + atomvect{A}.coreCharge...
                     .*(atomvect{A}.xyz - dipoleCenter);
             end
        end
        
        function cartesianMatrix = CalcCartesianMatrixByGTOExpansion(obj, stonG)
            atomvect = obj.molecule.atomVect;
            totalAtomNumber = length(atomvect);
            cartesianMatrix = zeros(obj.molecule.totalNumberAOs, obj.molecule.totalNumberAOs, 3);
            for A = 1:totalAtomNumber
                atomA  = atomvect{A};
                firstAOIndexA  = atomA.firstAOIndex;
                numValenceAOsA = length(atomA.valence);
                for a = 1:numValenceAOsA
                    mu = firstAOIndexA + a - 1;
                    for B = 1:totalAtomNumber
                        atomB  = atomvect{B};
                        firstAOIndexB  = atomB.firstAOIndex;
                        numValenceAOsB = length(atomB.valence);
                        for b = 1:numValenceAOsB
                            nu = firstAOIndexB + b - 1;
                            cartesianMatrix(mu, nu, :) = obj.CalcCartesianMatrixElementsByGTOExpansion(atomA, a, atomB, b, stonG);
                        end
                    end 
                end
            end
        
        end
        
        function xyz = CalcCartesianMatrixElementsByGTOExpansion(obj, atomA, valenceIndexA, atomB, valenceIndexB, stonG)
            xyz = zeros(3, 1);
            shellTypeA = atomA.valenceShellType;
            shellTypeB = atomB.valenceShellType;
            valenceOrbitalA = atomA.valence(valenceIndexA);
            valenceOrbitalB = atomB.valence(valenceIndexB);
            orbitalExponentA = obj.AtomGetOrbitalExponent(atomA, atomA.valenceShellType, valenceOrbitalA);
            orbitalExponentB = obj.AtomGetOrbitalExponent(atomB, atomB.valenceShellType, valenceOrbitalB);
            dxyz = atomA.xyz - atomB.xyz;
            rAB = norm(dxyz);
            for i = 1:stonG
                for j = 1:stonG
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
                    tempX = obj.GetGaussianCartesianMatrix(atomA.atomType, valenceOrbitalA, gaussianExponentA, atomA.xyz,...
                        atomB.atomType, valenceOrbitalB, gaussianExponentB, atomB.xyz, ...
                        rAB, overlapSASB, 1);
                    tempY = obj.GetGaussianCartesianMatrix(atomA.atomType, valenceOrbitalA, gaussianExponentA, atomA.xyz,...
                        atomB.atomType, valenceOrbitalB, gaussianExponentB, atomB.xyz, ...
                        rAB, overlapSASB, 2);
                    tempZ = obj.GetGaussianCartesianMatrix(atomA.atomType, valenceOrbitalA, gaussianExponentA, atomA.xyz,...
                        atomB.atomType, valenceOrbitalB, gaussianExponentB, atomB.xyz, ...
                        rAB, overlapSASB, 3);
                    xyz = xyz + temp.*[tempX; tempY; tempZ];
                end
            end
        end
        
        value = GetGaussianCartesianMatrix(obj, ...
                atomTypeA, valenceOrbitalA, gaussianExponentA, xyzA,...
                atomTypeB, valenceOrbitalB, gaussianExponentB, xyzB,...
                rAB, overlapSASB, axis)
        
        function overlapAOs = CalcOverlapAOs(obj)
            totalAONumber = obj.molecule.totalNumberAOs;
            totalAtomNumber = length(obj.molecule.atomVect);
            overlapAOs = eye(totalAONumber);
            atomvect = obj.molecule.atomVect;
            for A = 1:totalAtomNumber
                atomA = atomvect{A};
                symmetrize = false;
                for B = A+1:totalAtomNumber
                    atomB = atomvect{B};
                    diatomicOverlapAOs = obj.CalcDiatomicOverlapAOsInDiatomicFrame(atomA, atomB);
                    rotatingMatrix = obj.CalcRotatingMatrix(atomA, atomB);
%                     diatomicOverlapAOs = obj.RotateDiatmicOverlapAOsToSpaceFrame(diatomicOverlapAOs, rotatingMatrix);
                    diatomicOverlapAOs = rotatingMatrix * diatomicOverlapAOs * rotatingMatrix; % rotate
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
        
        % not tested
        value = GetGaussianOverlapAOs1stDerivative(obj, valenceOrbitalA, ...
            gaussianExponentA, ...
            valenceOrbitalB, ...
            gaussianExponentB,...
            dxyz,  rAB, ...
            axisA)
        
        function fockMatrix = CalcFockMatrix(obj, isGuess)
            totalNumberAOs   = obj.molecule.totalNumberAOs;
            atomvect = obj.molecule.atomVect;
            totalNumberAtoms = length(atomvect);
            fockMatrix = zeros(totalNumberAOs);
            
            fockoffdiag = zeros(totalNumberAOs);
            for A = 1:totalNumberAtoms
                atomA = atomvect{A};
                firstAOIndexA = atomA.firstAOIndex;
                lastAOIndexA  = atomA.GetLastAOIndex();
                
                % diag
                mu = firstAOIndexA:lastAOIndexA;
                fockMatrix(mu, mu) = diag(obj.GetFockDiagElement(atomA, mu, isGuess));
                
                % offdiag
                for B = A:totalNumberAtoms % upper part
                    atomB = atomvect{B};
                    firstAOIndexB = atomB.firstAOIndex;
                    lastAOIndexB  = atomB.GetLastAOIndex();
                    nu = firstAOIndexB:lastAOIndexB;
                    fockoffdiag(mu, nu) = obj.GetFockOffDiagElement(atomA, atomB, mu, nu, isGuess);
                end
            end % end of loop A
            fockoffdiag = fockoffdiag - tril(fockoffdiag);
            fockMatrix = fockMatrix + fockoffdiag + fockoffdiag';
        end
        
        function h1Matrix = CalcH1Matrix(obj)
            bkupOrbitalElectronPopulation = obj.orbitalElectronPopulation;
            bkupAtomicElectronPopulation = obj.atomicElectronPopulation;
            obj.orbitalElectronPopulation = zeros(size(obj.orbitalElectronPopulation));
            obj.atomicElectronPopulation = zeros(size(obj.atomicElectronPopulation));
            h1Matrix = obj.CalcFockMatrix(false);
            obj.orbitalElectronPopulation = bkupOrbitalElectronPopulation;
            obj.atomicElectronPopulation = bkupAtomicElectronPopulation;
        end
        
        function res = RotateDiatmicOverlapAOsToSpaceFrame(~, diatomicOverlapAOs, rotatingMatrix)
            res = rotatingMatrix * diatomicOverlapAOs * rotatingMatrix;
        end
        
        function overlapAOs = SetOverlapAOsElement(~, overlapAOs, diatomicOverlapAOs, atomA, atomB, isSymmetricOverlapAOs)
            if(nargin < 5)
                isSymmetricOverlapAOs = true;
            end
            firstAOIndexAtomA = atomA.firstAOIndex;
            firstAOIndexAtomB = atomB.firstAOIndex;
            
            for i = 1:length(atomA.valence)
                orbitalA = atomA.valence(i);
                for j = 1:length(atomB.valence)
                    orbitalB = atomB.valence(j);
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
%             value = exp(-1.0*rho)*factorial(k);
            value = exp(-1.0*rho);
            for mu = 1:k+1
                tmp2=rho^mu; % tmp2 = pow(rho,mu)
%                 tmp1 = tmp1 + 1.0/(tmp2*factorial(k-mu+1));
                tmp1 = tmp1 + 1.0/(tmp2/prod(k-mu+2:k));
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
%                     tmp1 = tmp1 + factorial(k)/factorial(k-mu+1)     /tmp3;
%                     tmp2 = tmp2 + factorial(k)/factorial(k-mu+1)*tmp4/tmp3;
                    tmp1 = tmp1 + prod(k-mu+2:k)     /tmp3;
                    tmp2 = tmp2 + prod(k-mu+2:k)*tmp4/tmp3;
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
        
        % not tested
        function res = GetAuxiliaryA1stDerivative(obj, k, rho)
            res = -1.0*obj.GetAuxiliaryA(k+1, rho);
        end
        % not tested
        function res = GetAuxiliaryA2ndDerivative(obj, k, rho)
            res = obj.GetAuxiliaryA(k+2, rho);
        end
        % not tested
        function res = GetAuxiliaryB1stDerivative(obj, k, rho)
            res = -1.0*obj.GetAuxiliaryB(k+1, rho);
        end
        % not tested
        function res = GetAuxiliaryB2ndDerivative(obj, k, rho)
            res = obj.GetAuxiliaryB(k+2, rho);
        end
        
        %    void DoDamp(double rmsDensity,
        %                bool&  hasAppliedDamping,
        %                double** orbitalElectronPopulation,
        %                double const* const* oldOrbitalElectronPopulation,
        %                const MolDS_base::Molecule& molecule) const;
        
        function [orbitalElectronPopulation, ...
                diisStoredDensityMatrix, ...
                diisStoredErrorVect, ...
                diisErrorProducts, ...
                diisErrorCoefficients]...
                = DoDIIS(obj, orbitalElectronPopulation,...
                oldOrbitalElectronPopulation,...
                diisStoredDensityMatrix,...
                diisStoredErrorVect,...
                diisErrorProducts,...
                diisErrorCoefficients,...
                diisNumErrorVect,...
                iterStep)
            % diis start
            totalNumberAOs = obj.molecule.totalNumberAOs;
            diisStartError = Arguments.GetInstance().diisStartErrorSCF;
            diisEndError = Arguments.GetInstance().diisEndErrorSCF;
            
            if( 0 < diisNumErrorVect)
                diisStoredDensityMatrix(:, 1:diisNumErrorVect-1) = diisStoredDensityMatrix(:, 2:diisNumErrorVect);
                diisStoredDensityMatrix(:, diisNumErrorVect) = 0*diisStoredDensityMatrix(:, diisNumErrorVect);
                diisStoredErrorVect(:, 1:diisNumErrorVect-1) = diisStoredErrorVect(:, 2:diisNumErrorVect);
                diisStoredErrorVect(:, diisNumErrorVect) = 0*diisStoredErrorVect(:, diisNumErrorVect);
            end
            
            diisStoredDensityMatrix(:, diisNumErrorVect) = reshape(orbitalElectronPopulation, totalNumberAOs*totalNumberAOs, 1);
            diisStoredErrorVect(:, diisNumErrorVect) = reshape(orbitalElectronPopulation - oldOrbitalElectronPopulation, totalNumberAOs*totalNumberAOs, 1);
            
%             for mi = 1:diisNumErrorVect-1
%                 for mj = 1:diisNumErrorVect-1
%                     diisErrorProducts(mi,mj) = diisErrorProducts(mi+1,mj+1);
%                 end
%             end
            
            mi = (1:diisNumErrorVect-1)';
            diisErrorProducts(mi,mi) = diisErrorProducts(mi+1,mi+1);
            
            diisErrorProducts(1:end-1, diisNumErrorVect) = diisStoredErrorVect' * diisStoredErrorVect(:, diisNumErrorVect);
            
            for mi = 1:diisNumErrorVect
                diisErrorProducts(diisNumErrorVect, mi) = diisErrorProducts(mi, diisNumErrorVect);
                diisErrorProducts(diisNumErrorVect+1, mi) = -1.0;
                diisErrorProducts(mi, diisNumErrorVect+1) = -1.0;
                diisErrorCoefficients(mi) = 0.0;
            end
            
            diisErrorProducts(diisNumErrorVect+1,diisNumErrorVect+1) = 0.0;
            diisErrorCoefficients(diisNumErrorVect+1) = -1.0;
            diisError = max(abs(diisStoredErrorVect(diisNumErrorVect,:)));
            if(diisNumErrorVect <= iterStep && diisEndError<diisError && diisError<diisStartError)
                tmpDiisErrorProducts = diisErrorProducts;
%                 if(abs(det(tmpDiisErrorProducts)) < 1e-12)
%                     hasAppliedDIIS = false;
%                     continue;
%                 end
                diisErrorCoefficients = tmpDiisErrorProducts \ diisErrorCoefficients;
                orbitalElectronPopulation = reshape(diisStoredDensityMatrix * diisErrorCoefficients(1:end-1), totalNumberAOs, totalNumberAOs);
            end
            % diis end
        end
        
        function CheckEnableAtomType(obj)
            for i = 1:length(obj.molecule.atomVect)
                atomvect = obj.molecule.atomVect;
                atomType = atomvect{i}.atomType;
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
        
        function CheckNumberValenceElectrons(obj)
            if(mod(obj.molecule.totalNumberValenceElectrons, 2) == 1)
                throw(MException('Cndo2:CheckNumberValenceElectrons', 'Odd number of electrons.'));
            end
        end
        
        function elecSCFEnergy = CalcElecSCFEnergy(obj)
            % use density matrix for electronic energy
            isGuess = false;
            fMatrix = obj.CalcFockMatrix(isGuess);
            hMatrix = obj.CalcH1Matrix();
            electronicEnergy = sum(sum(obj.orbitalElectronPopulation .* (hMatrix+fMatrix)));
            electronicEnergy = electronicEnergy * 0.5;
            
            elecSCFEnergy = electronicEnergy...
                +obj.coreRepulsionEnergy...
                +obj.coreEpcCoulombEnergy...
                +obj.vdWCorrectionEnergy;
        end
        
        function res = AtomGetOrbitalExponent(~, atom, shellType, orbitalType)
            if(shellType == 1 && orbitalType == 1) % kShell, s
                res = atom.paramPool.effectiveNuclearChargeK/atom.GetEffectivePrincipalQuantumNumber(shellType);
            elseif(shellType == 2 && (orbitalType == 1  || ...
                    orbitalType == 4 || ...
                    orbitalType == 2 || ...
                    orbitalType == 3)) % lShell, s py pz px
                res = atom.paramPool.effectiveNuclearChargeL/atom.GetEffectivePrincipalQuantumNumber(shellType);
            elseif(shellType == 3 && (orbitalType == 1  || ...
                    orbitalType == 4 || ...
                    orbitalType == 2 || ...
                    orbitalType == 3 )) % mShell, s py pz px
                res = atom.paramPool.effectiveNuclearChargeMsp/atom.GetEffectivePrincipalQuantumNumber(shellType);
            elseif(shellType == 3 && (orbitalType == 5  || ...
                    orbitalType == 6 ||...
                    orbitalType == 7 ||...
                    orbitalType == 8 ||...
                    orbitalType == 9)) % mShell, dxy dyz dzz dzx dxxyy
                res = atom.paramPool.effectiveNuclearChargeMd/atom.GetEffectivePrincipalQuantumNumber(shellType);
            elseif(shellType == 4 && (1  || ...
                    orbitalType == 4 || ...
                    orbitalType == 2 || ...
                    orbitalType == 3 )) % nShell, dxy dyz dzz dzx dxxyy
                res = atom.paramPool.effectiveNuclearChargeNsp/atom.GetEffectivePrincipalQuantumNumber(shellType);
            else
                throw(MException('Cndo2:AtomGetOrbitalExponent', 'Shell/Orbital type wrong.'));
            end
        end
        
%         function res = AtomGetCoreIntegral(~, atom, orbital, gamma, isGuess) % OrbitalType orbital
%             if(orbital == 1) % s
%                 res = -1.0*atom.imuAmuS;
%             elseif(orbital == 4 || orbital == 2 || orbital == 3) % py pz px
%                 res = -1.0*atom.imuAmuP;
%             elseif(orbital == 5 || ...
%                     orbital == 6 || ...
%                     orbital == 7 || ...
%                     orbital == 8 || ...
%                     orbital == 9 ) % dxy dyz dzz dzx dxxyy
%                 res = -1.0*atom.imuAmuD;
%             else
%                 throw(MException('CNDO2:AtomGetCndo2CoreIntegral', 'Orbital type wrong.'));
%             end
%             if(~isGuess)
%                 res = res - (atom.coreCharge - 0.5)*gamma;
%             end
%         end
        
        % vectorized version of the above function
        function res = AtomGetCoreIntegral(~, atom, orbital, gamma, isGuess) % OrbitalType orbital
            res = zeros(length(orbital), 1);
            for i = 1:length(orbital)
                if(orbital(i) == 1) % s
                    res(i) = -1.0*atom.paramPool.imuAmuS;
                elseif(orbital(i) == 4 || orbital(i) == 2 || orbital(i) == 3) % py pz px
                    res(i) = -1.0*atom.paramPool.imuAmuP;
                elseif(orbital(i) == 5 || ...
                        orbital(i) == 6 || ...
                        orbital(i) == 7 || ...
                        orbital(i) == 8 || ...
                        orbital(i) == 9 ) % dxy dyz dzz dzx dxxyy
                    res(i) = -1.0*atom.paramPool.imuAmuD;
                else
                    throw(MException('CNDO2:AtomGetCndo2CoreIntegral', 'Orbital type wrong.'));
                end
            end
            if(~isGuess)
                res = res - (atom.coreCharge - 0.5)*gamma;
            end
        end
        
        function res = AtomGetBondingParameter(~, atom)
            res = atom.paramPool.bondingParameter;
        end
        
    end
    
end




