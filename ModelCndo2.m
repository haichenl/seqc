classdef ModelCndo2 < handle % wrapper of cndo2
    
    properties
        
        cndo2;
        
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
        
        function obj = ModelCndo2(zmat_in)
            obj.frag.config.zmat = zmat_in;
            mol = Molecule();
            cart = zmat_in.GauCart();
            for atomi = 1:length(zmat_in.atoms)
                xyz = cart(atomi, :) / 0.5291772;
                mol.AddAtom(AtomFactory.Create(EnumAtom(zmat_in.atoms{atomi}.z), atomi, xyz));
            end
            mol.CalcBasics();
            mol.AddParamPoolsIntoManager(ParamPoolManagerCndo2.GetInstance());
            
            obj.cndo2 = Cndo2();
            obj.cndo2.SetMolecule(mol);
            
            obj.nelec = mol.totalNumberValenceElectrons;
        end
        
        function solveHF(obj, ~, ~, ~, ~)
            obj.cndo2.DoSCF();
            obj.Eorb = obj.cndo2.energiesMO;
        end
        
        function res = Etot(obj, ~)
            res = obj.cndo2.elecSCFEnergy;
        end
        
        function res = E1(obj, ~)
            res = sum(sum(obj.cndo2.h1Matrix.*obj.cndo2.orbitalElectronPopulation));
        end
        
        function res = density(obj, ~)
            res = obj.cndo2.orbitalElectronPopulation;
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