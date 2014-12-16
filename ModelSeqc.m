classdef ModelSeqc < handle % wrapper of seqc
    
    properties
        
        seqc;
        
        % the below 2's are for computing homo/lumo
        nelec;
        Eorb;
        
        % fake properties
        nenv = 0;
        nenvExt = 0;
        frag; % need frag.config.zmat
        aType;
        mixers;
        index;
        densitySave;
        
    end
    
    methods
        
        function obj = ModelSeqc(zmat_in, class, paramManager)
            obj.frag.config.zmat = zmat_in;
            mol = SEQC.Molecule();
            cart = zmat_in.GauCart();
            for atomi = 1:length(zmat_in.atoms)
                xyz = cart(atomi, :) / 0.5291772;
                mol.AddAtom(SEQC.AtomFactory.Create(SEQC.EnumAtom(zmat_in.atoms{atomi}.z), atomi, xyz));
            end
            mol.CalcBasics();
            mol.AddParamPoolsIntoManager(paramManager);
            
            obj.seqc = class();
            obj.seqc.SetMolecule(mol);
            
            obj.nelec = mol.totalNumberValenceElectrons;
        end
        
        function solveHF(obj, ~, ~, ~, ~)
            obj.seqc.DoSCF();
            obj.Eorb = obj.seqc.energiesMO;
        end
        
        function res = Etot(obj, ~)
            res = obj.seqc.elecSCFEnergy;
        end
        
        function res = E1(obj, ~)
            res = sum(sum(obj.seqc.h1Matrix.*obj.seqc.orbitalElectronPopulation));
        end
        
        function res = density(obj, ~)
            res = obj.seqc.orbitalElectronPopulation;
        end
        
        % fake methods
        function clearModifiers(~)
        end
        function cacheOperators(~, ~)
        end
        function res = getCache(~)
            res = [];
        end
        function setCache(~, ~)
        end
        function clearCache(~, ~)
        end
        
    end
    
end