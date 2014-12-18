classdef LowLevel < FITSEQC.HighLevel
    
    properties (Access = private)
        
        seqc;
        
    end
    
    methods
        
        function obj = LowLevel(input)
            import SEQC.AtomFactory SEQC.TheoryFactory ...
                SEQC.MoleculeSEQC SEQC.EnumAtom;
            
            molecule = input.molecule;
            seqcTheory = input.seqcTheory;
            fieldStrengths = input.fieldStrengths;
            paramManager = input.paramManager;
            
            cartesianInBohr = molecule.CartesianGeometryInBohrRadius();
            moleculeSEQC = MoleculeSEQC();
            for iAtom = 1:molecule.NumAtoms()
                atom = AtomFactory.Create( ...
                    EnumAtom(molecule.cartesian(iAtom,1)), ...
                    iAtom, cartesianInBohr(iAtom, :));
                moleculeSEQC.AddAtom(atom);
            end
            moleculeSEQC.CalcBasics();
            
            obj.seqc = TheoryFactory.Create(seqcTheory);
            obj.seqc.SetMolecule(moleculeSEQC);
            obj.seqc.SetEnvironmentFieldStrengths(fieldStrengths)
            
            for atom = obj.seqc.molecule.atomVect
                paramManager.AddParamPool(atom{1}.paramPool);
            end
        end
        
        function DoSCF(obj)
            obj.seqc.DoSCF();
            obj.energyTotal = obj.seqc.elecSCFEnergy;
            numElectrons = obj.seqc.molecule.totalNumberValenceElectrons;
            obj.energyHOMO = obj.seqc.energiesMO(numElectrons./2);
            obj.energyLUMO = obj.seqc.energiesMO(numElectrons./2 + 1);
            obj.energyEnvironment = obj.seqc.EnvironmentFieldEnergy();
        end
        
    end
    
end