classdef Molecule < handle
    
    properties (SetAccess = private)
        
        atomVect = {};
        epcVect = {};      % Vector of Environmental Point Charges
        xyzCOM; % x, y, z coordinates of the center of atomic mass;
        xyzCOC; % x, y, z coordinates of the center of core's mass;
        distanceAtoms;    % distance between each atom;
        distanceEpcs;     % distance between each environmental point charge;
        distanceAtomsEpcs;% distance between each atom and environmental point charge;
        totalNumberShells;
        totalNumberAOs;
        totalNumberValenceElectrons;
        totalCoreMass;
        
    end
    
    methods
        
        function obj = Molecule()
        end
        
        function AddAtom(obj, atom)
            obj.atomVect{end+1} = atom;
        end
        
        function AddEpc(obj, epc)
            obj.epcVect{end+1} = epc;
        end
        
        function AddParamPoolsIntoManager(obj, paramPoolManager)
            for i = 1:length(obj.atomVect)
                paramPoolManager.AddParamPool(obj.atomVect{i}.paramPool);
            end
        end
        
        function res = GetXyzDipoleCenter(obj)
            res = obj.xyzCOC;
        end
        
        function CalcBasics(obj)
            obj.CalcTotalNumberAOs();
            obj.CalcTotalNumberShells();
            obj.CalcTotalNumberValenceElectrons();
            obj.CalcTotalCoreMass();
            obj.CalcBasicsConfiguration();
        end
        
        function CalcBasicsConfiguration(obj)
            obj.CalcXyzCOM();
            obj.CalcXyzCOC();
            obj.CalcDistanceAtoms();
            obj.CalcDistanceEpcs();
            obj.CalcDistanceAtomsEpcs();
        end
        
        function [inertiaMoments, inertiaTensor] = CalcPrincipalAxes(obj)
            obj.CalcXyzCOM();
            inertiaTensorOrigin = obj.xyzCOM;
            inertiaTensor = obj.CalcInertiaTensor(inertiaTensorOrigin);
            inertiaMoments = eig(inertiaTensor);
        end
        
        function res = GetDistanceAtoms(obj, atomA, atomB)
            res = obj.distanceAtoms(atomA.index, atomB.index);
        end
        
%         function res = GetDistanceAtoms(obj, arg1, arg2)
%             if(isnumeric(arg1) && isnumeric(arg2))
%                 indexAtomA = arg1;
%                 indexAtomB = arg2;
%             elseif(isa(arg1, 'Atom') && isa(arg2, 'Atom'))
%                 indexAtomA = arg1.index;
%                 indexAtomB = arg2.index;
%             else
%                 throw(MException('Molecule:GetDistanceAtoms', 'Input argument type wrong.'));
%             end
%             res = obj.distanceAtoms(indexAtomA, indexAtomB);
%         end
        
        function res = GetDistanceEpcs(obj, arg1, arg2)
            if(isnumeric(arg1) && isnumeric(arg2))
                indexEpcA = arg1;
                indexEpcB = arg2;
            elseif(isa(arg1, 'Atom') && isa(arg2, 'Atom'))
                indexEpcA = arg1.index;
                indexEpcB = arg2.index;
            else
                throw(MException('Molecule:GetDistanceEpcs', 'Input argument type wrong.'));
            end
            res = obj.distanceEpcs(indexEpcA, indexEpcB);
        end
        
        function res = GetDistanceAtomEpc(obj, arg1, arg2)
            if(isnumeric(arg1) && isnumeric(arg2))
                indexAtom = arg1;
                indexEpc = arg2;
            elseif(isa(arg1, 'Atom') && isa(arg2, 'Atom'))
                indexAtom = arg1.index;
                indexEpc = arg2.index;
            else
                throw(MException('Molecule:GetDistanceAtomEpc', 'Input argument type wrong.'));
            end
            res = obj.distanceAtomsEpcs(indexAtom, indexEpc);
        end
        
    end
    
    methods (Access = private)
        
        function CalcTotalNumberAOs(obj)
            obj.totalNumberAOs = 0; 
            for i = 1:length(obj.atomVect)
                obj.atomVect{i}.SetFirstAOIndex(obj.totalNumberAOs + 1);
                obj.totalNumberAOs = obj.totalNumberAOs + obj.atomVect{i}.GetValenceSize();
            end
        end
        
        function CalcTotalNumberShells(obj)
            obj.totalNumberShells = 0; 
            for i = 1:length(obj.atomVect)
                obj.atomVect{i}.SetFirstShellIndex(obj.totalNumberShells + 1);
                obj.totalNumberShells = obj.totalNumberShells + obj.atomVect{i}.nShell;
            end
        end
        
        function CalcTotalNumberValenceElectrons(obj)
            obj.totalNumberValenceElectrons = 0;
            for i = 1:length(obj.atomVect)
                obj.totalNumberValenceElectrons = obj.totalNumberValenceElectrons + obj.atomVect{i}.numberValenceElectrons;
            end
        end
        
        function CalcTotalCoreMass(obj)
            obj.totalCoreMass = 0; 
            for i = 1:length(obj.atomVect)
                obj.totalCoreMass = obj.totalCoreMass + obj.atomVect{i}.GetCoreMass();
            end
        end
        
        function CalcXyzCOM(obj)
            totalAtomicMass = 0.0;
            obj.xyzCOM = zeros(3, 1);
            for i = 1:length(obj.atomVect)
                atomicXyz = obj.atomVect{i}.xyz;
                atomicMass = obj.atomVect{i}.atomicMass;
                totalAtomicMass = totalAtomicMass + atomicMass;
                obj.xyzCOM = obj.xyzCOM + atomicXyz .* atomicMass;
            end
            obj.xyzCOM = obj.xyzCOM ./ totalAtomicMass;
        end
        
        function CalcXyzCOC(obj)
            totalCoreMass_ = 0.0;
            obj.xyzCOC = zeros(3, 1);
            for i = 1:length(obj.atomVect)
                atomicXyz = obj.atomVect{i}.xyz;
                coreMass = obj.atomVect{i}.GetCoreMass();
                totalCoreMass_ = totalCoreMass_ + coreMass;
                obj.xyzCOC = obj.xyzCOC + atomicXyz .* coreMass;
            end
            obj.xyzCOC = obj.xyzCOC ./ totalCoreMass_;
        end
        
        function CalcDistanceAtoms(obj)
            obj.distanceAtoms = zeros(length(obj.atomVect));
            for a = 1:length(obj.atomVect)
                atomA = obj.atomVect{a};
                for b = a:length(obj.atomVect)
                    atomB = obj.atomVect{b};
                    obj.distanceAtoms(a, b) = norm(atomA.xyz - atomB.xyz);
                end
            end
            obj.distanceAtoms = obj.distanceAtoms + obj.distanceAtoms' - diag(diag(obj.distanceAtoms));
        end
        
        function CalcDistanceEpcs(obj)
            obj.distanceEpcs = zeros(length(obj.epcVect));
            for a = 1:length(obj.epcVect)
                epcA = obj.epcVect{a};
                for b = a:length(obj.epcVect)
                    epcB = obj.epcVect{b};
                    obj.distanceEpcs(a, b) = norm(epcA.xyz - epcB.xyz);
                end
            end
            obj.distanceEpcs = obj.distanceEpcs + obj.distanceEpcs' - diag(diag(obj.distanceEpcs));
        end
        
        function CalcDistanceAtomsEpcs(obj)
            obj.distanceAtomsEpcs = zeros(length(obj.atomVect), length(obj.epcVect));
            for a = 1:length(obj.atomVect)
                atom = obj.atomVect{a};
                for b = 1:length(obj.epcVect)
                    epc = obj.epcVect{b};
                    obj.distanceAtomsEpcs(a, b) = norm(atom.xyz - epc.xyz);
                end
            end
        end
        
        function inertiaTensor = CalcInertiaTensor(obj, inertiaTensorOrigin)
            inertiaTensor = zeros(3, 3);
            for a = 1:length(obj.atomVect)
                atomicMass = obj.atomVect{a}.atomicMass;
                xyz = obj.atomVect{a}.xyz;
                x = xyz(1) - inertiaTensorOrigin(1);
                y = xyz(2) - inertiaTensorOrigin(2);
                z = xyz(3) - inertiaTensorOrigin(3);
                
                inertiaTensor(1, 1) = inertiaTensor(1, 1) + atomicMass*(y*y + z*z);
                inertiaTensor(1, 2) = inertiaTensor(1, 2) - atomicMass*x*y;
                inertiaTensor(1, 3) = inertiaTensor(1, 3) - atomicMass*x*z;
                
                inertiaTensor(2, 1) = inertiaTensor(2, 1) - atomicMass*y*x;
                inertiaTensor(2, 2) = inertiaTensor(2, 2) + atomicMass*(x*x + z*z);
                inertiaTensor(2, 3) = inertiaTensor(2, 3) - atomicMass*y*z;
                
                inertiaTensor(3, 1) = inertiaTensor(3, 1) - atomicMass*z*x;
                inertiaTensor(3, 2) = inertiaTensor(3, 2) - atomicMass*z*y;
                inertiaTensor(3, 3) = inertiaTensor(3, 3) + atomicMass*(x*x + y*y);

            end
        end
        
    end
    
end



