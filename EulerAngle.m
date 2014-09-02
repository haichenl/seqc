classdef EulerAngle < handle
    properties (SetAccess = private)
        alpha;
        beta;
        gamma;
    end
    
    methods
        
        function obj = EulerAngle(xyz)
            r = norm(xyz);
            obj.beta = acos(xyz(3) / r);
            if(xyz(1) == 0 && xyz(2) == 0)
                obj.alpha = 0;
            else
                r = norm(xyz(1:2));
                obj.alpha = atan2(xyz(2)/r, xyz(1)/r);
            end
            obj.gamma = 0;
        end
        
%         function res = GetAlpha(obj)
%             res = obj.alpha;
%         end
%         
%         function res = GetBeta(obj)
%             res = obj.beta;
%         end
%         
%         function res = GetGamma(obj)
%             res = obj.gamma;
%         end
        
        function SetAlpha(obj, alpha_in)
            obj.alpha = alpha_in;
        end
        
        function SetBeta(obj, beta_in)
            obj.beta = beta_in;
        end
        
        function SetGamma(obj, gamma_in)
            obj.gamma = gamma_in;
        end
        
    end
end