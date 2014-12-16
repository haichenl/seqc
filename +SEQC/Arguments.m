classdef (Sealed) Arguments < handle
    
    properties (Constant)
        
        eV2AU         = 0.03674903;
        j2AU          = 1.0e18/4.35974394;
        kcalMolin2AU  = 0.00159360175;
        angstrom2AU   = 1.0/0.5291772;
        nm2AU         = 10.0*1.0/0.5291772;
        kayser2AU     = 4.556336e-6;
        nmin2AU       = 4.556336e01;
        fs2AU         = 1.0/(2.418884326505e-2);
        gMolin2AU     = 1.0e5/(6.0221415*9.1095);
        degree2Radian = pi / 180.0;
        boltzmann     = 3.166791e-6;
        avogadro      = 6.0221415e23;
        debye2AU      = 0.393430191;
        
        vdWScalingFactorSCFPM3DAM1D = 1.40;
        vdWDampingFactorSCFPM3DAM1D = 23.0;
        
    end
    
    properties (SetAccess = private)
        
        currentTheory;
        
        % scf
        thresholdSCF;
        maxIterationsSCF;
        dampingThreshSCF;
        dampingWeightSCF;
        diisNumErrorVectSCF;
        diisStartErrorSCF;
        diisEndErrorSCF;
        requiresVdWSCF;
        sumChargesIndexPairsSCF = {};
        vdWScalingFactorSCF;
        vdWDampingFactorSCF;
        
    end
    
    methods (Access = public)
        function res = GetEV2AU(obj)
            res = obj.eV2AU;
        end
        function res = GetJ2AU(obj)
            res = obj.j2AU;
        end
        function res = GetKcalMolin2AU(obj)
            res = obj.kcalMolin2AU;
        end
        function res = GetAngstrom2AU(obj)
            res = obj.angstrom2AU;
        end
        function res = GetNm2AU(obj)
            res = obj.nm2AU;
        end
        function res = GetKayser2AU(obj)
            res = obj.kayser2AU;
        end
        function res = GetNmin2AU(obj)
            res = obj.nmin2AU;
        end
        function res = GetFs2AU(obj)
            res = obj.fs2AU;
        end
        function res = GetGMolin2AU(obj)
            res = obj.gMolin2AU;
        end
        function res = GetDegree2Radian(obj)
            res = obj.degree2Radian;
        end
        function res = GetBoltzmann(obj)
            res = obj.boltzmann;
        end
        function res = GetAvogadro(obj)
            res = obj.avogadro;
        end
        function res = GetDebye2AU(obj)
            res = obj.debye2AU;
        end
        function res = GetThresholdSCF(obj)
            res = obj.thresholdSCF;
        end
        function SetThresholdSCF(obj, threshold)
            obj.thresholdSCF = threshold;
        end
        function res = GetMaxIterationsSCF(obj)
            res = obj.maxIterationsSCF;
        end
        function SetMaxIterationsSCF(obj, maxIter)
            obj.maxIterationsSCF = maxIter;
        end
        function res = GetDampingThreshSCF(obj)
            res = obj.dampingThreshSCF;
        end
        function SetDampingThreshSCF(obj, dThresh)
            obj.dampingThreshSCF = dThresh;
        end
        function res = GetDampingWeightSCF(obj)
            res = obj.dampingWeightSCF;
        end
        function SetDampingWeightSCF(obj, dWeight)
            obj.dampingWeightSCF = dWeight;
        end
        function res = GetDiisNumErrorVectSCF(obj)
            res = obj.diisNumErrorVectSCF;
        end
        function SetDiisNumErrorVectSCF(obj, numEVect)
            obj.diisNumErrorVectSCF = numEVect;
        end
        function res = GetDiisStartErrorSCF(obj)
            res = obj.diisStartErrorSCF;
        end
        function SetDiisStartErrorSCF(obj, sError)
            obj.diisStartErrorSCF = sError;
        end
        function res = GetDiisEndErrorSCF(obj)
            res = obj.diisEndErrorSCF;
        end
        function SetDiisEndErrorSCF(obj, eError)
            obj.diisEndErrorSCF = eError;
        end
        function res = RequiresSumChargesSCF(obj)
            res = (~isempty(obj.sumChargesIndexPairsSCF)) && (0<length(obj.sumChargesIndexPairsSCF));
        end
        function res = GetSumChargesIndexPairsSCF(obj)
            res = obj.sumChargesIndexPairsSCF;
        end
        function AddSumChargesIndexPairsSCF(obj, firstAtomIndex, lastAtomIndex)
            atomIndexPair.firstAtomIndex = firstAtomIndex;
            atomIndexPair.lastAtomIndex = lastAtomIndex;
            obj.sumChargesIndexPairsSCF{end+1} = atomIndexPair;
        end
        function res = RequiresVdWSCF(obj)
            res = obj.requiresVdWSCF;
        end
        function SetRequiresVdWSCF(obj, requires)
            obj.requiresVdWSCF = requires;
        end
        function res = GetVdWScalingFactorSCF(obj)
            res = obj.vdWScalingFactorSCF;
        end
        function SetVdWScalingFactorSCF(obj, vdWScal)
            if(nargin < 1)
                obj.vdWScalingFactorSCF = obj.vdWScalingFactorSCFPM3DAM1D;
            else
                obj.vdWScalingFactorSCF = vdWScal;
            end
        end
        function SetVdWDampingFactorSCF(obj, vdWDamp)
            if(nargin < 1)
                obj.vdWDampingFactorSCF = obj.vdWDampingFactorSCFPM3DAM1D;
            else
                obj.vdWDampingFactorSCF = vdWDamp;
            end
        end
        
        function SetCurrentTheory(obj, currTheory)
            obj.currentTheory = currTheory;
        end
        function res = GetCurrentTheory(obj)
            res = obj.currentTheory;
        end
        
    end
    
    methods (Access = private)
        
        function obj = Arguments()
            obj.SetDefaultValues();
        end
        
        function SetDefaultValues(obj)
            obj.currentTheory = SEQC.EnumTheory.CNDO2;
            % SCF
            obj.thresholdSCF        = 1.0e-8;
            obj.maxIterationsSCF    = 100;
            obj.dampingThreshSCF    = 1.0;
            obj.dampingWeightSCF    = 0.8;
            obj.diisNumErrorVectSCF = 5;
            obj.diisStartErrorSCF   = 1.0e-2;
            obj.diisEndErrorSCF     = 1.0e-8;
            obj.requiresVdWSCF      = false;
            obj.vdWScalingFactorSCF = 1.40;
            obj.vdWDampingFactorSCF = 23.0;
        end
        
    end
    
    methods (Static)
        
        function singleObj = GetInstance()
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = SEQC.Arguments();
            end
            singleObj = localObj;
        end
        
    end
    
end