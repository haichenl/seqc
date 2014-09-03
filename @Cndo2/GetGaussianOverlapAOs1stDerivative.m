function value = GetGaussianOverlapAOs1stDerivative(obj, valenceOrbitalA, ...
    gaussianExponentA, ...
    valenceOrbitalB, ...
    gaussianExponentB,...
    dxyz,  rAB, ...
    axisA)
value = 0.0;
gauPlusAB = gaussianExponentA+gaussianExponentB;
gauMultAB = gaussianExponentA*gaussianExponentB;
dx = dxyz(1);
dy = dxyz(2);
dz = dxyz(3);
if(valenceOrbitalA == 1 && valenceOrbitalB == 1)
    temp = -2.0*gauMultAB/gauPlusAB;
    value = temp;
    if(axisA == 1)
        value = value * dx;
    elseif(axisA == 2)
        value = value * dy;
    elseif(axisA == 3)
        value = value * dz;
    end
elseif(valenceOrbitalA == 1 && valenceOrbitalB == 4)
    temp1 = 4.0*(gaussianExponentA*gaussianExponentA)*gaussianExponentB*sqrt(gaussianExponentB)...
        /(gauPlusAB*gauPlusAB);
    if(axisA == 1)
        temp2 = 2.0*gaussianExponentA*sqrt(gaussianExponentB)...
            /gauPlusAB;
        value = temp2-temp1*dx*dx;
    elseif(axisA == 2)
        value = -1.0*temp1*dx*dy;
    elseif(axisA == 3)
        value = -1.0*temp1*dx*dz;
    end
elseif(valenceOrbitalA == 1 && valenceOrbitalB == 2)
    temp1 = 4.0*(gaussianExponentA*gaussianExponentA)*gaussianExponentB*sqrt(gaussianExponentB)...
        /(gauPlusAB*gauPlusAB);
    if(axisA == 1)
        value = -1.0*temp1*dx*dy;
    elseif(axisA == 2)
        temp2 = 2.0*gaussianExponentA*sqrt(gaussianExponentB)...
            /gauPlusAB;
        value = temp2-temp1*dy*dy;
    elseif(axisA == 3)
        value = -1.0*temp1*dy*dz;
    end
elseif(valenceOrbitalA == 1 && valenceOrbitalB == 3)
    temp1 = 4.0*(gaussianExponentA*gaussianExponentA)*gaussianExponentB*sqrt(gaussianExponentB)...
        /(gauPlusAB*gauPlusAB);
    if(axisA == 1)
        value = -1.0*temp1*dx*dz;
    elseif(axisA == 2)
        value = -1.0*temp1*dy*dz;
    elseif(axisA == 3)
        temp2 = 2.0*gaussianExponentA*sqrt(gaussianExponentB)/gauPlusAB;
        value = temp2-temp1*dz*dz;
    end
elseif(valenceOrbitalA == 4 && valenceOrbitalB == 1)
    temp1 = 4.0*gaussianExponentA*sqrt(gaussianExponentA)...
        *gaussianExponentB*gaussianExponentB...
        /(gauPlusAB*gauPlusAB);
    if(axisA == 1)
        temp2 = 2.0*sqrt(gaussianExponentA)*gaussianExponentB/gauPlusAB;
        value = -1.0*temp2+temp1*dx*dx;
    elseif(axisA == 2)
        value = temp1*dx*dy;
    elseif(axisA == 3)
        value = temp1*dx*dz;
    end
elseif(valenceOrbitalA == 2 && valenceOrbitalB == 1)
    temp1 = 4.0*gaussianExponentA*sqrt(gaussianExponentA)...
        *gaussianExponentB*gaussianExponentB...
        /(gauPlusAB*gauPlusAB);
    if(axisA == 1)
        value = temp1*dx*dy;
    elseif(axisA == 2)
        temp2 = 2.0*sqrt(gaussianExponentA)*gaussianExponentB/gauPlusAB;
        value = -1.0*temp2+temp1*dy*dy;
    elseif(axisA == 3)
        value = temp1*dy*dz;
    end
elseif(valenceOrbitalA == 3 && valenceOrbitalB == 1)
    temp1 = 4.0*gaussianExponentA*sqrt(gaussianExponentA)...
        *gaussianExponentB*gaussianExponentB...
        /(gauPlusAB*gauPlusAB);
    if(axisA == 1)
        value = temp1*dx*dz;
    elseif(axisA == 2)
        value = temp1*dy*dz;
    elseif(axisA == 3)
        temp2 = 2.0*sqrt(gaussianExponentA)*gaussianExponentB/gauPlusAB;
        value = -1.0*temp2+temp1*dz*dz;
    end
elseif(valenceOrbitalA == 4 && valenceOrbitalB == 2)
    temp1 = 8.0*(gauMultAB*gauMultAB*sqrt(gauMultAB))...
        /(gauPlusAB*gauPlusAB*gauPlusAB);
    temp2 = 4.0*gauMultAB*sqrt(gauMultAB)...
        /(gauPlusAB*gauPlusAB);
    if(axisA == 1)
        value = -1.0*temp2*dy+temp1*dx*dx*dy;
    elseif(axisA == 2)
        value = -1.0*temp2*dx+temp1*dx*dy*dy;
    elseif(axisA == 3)
        value = temp1*dx*dy*dz;
    end
elseif(valenceOrbitalA == 2 && valenceOrbitalB == 4)
    temp1 = 8.0*(gauMultAB*gauMultAB*sqrt(gauMultAB))...
        /(gauPlusAB*gauPlusAB*gauPlusAB);
    temp2 = 4.0*gauMultAB*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB);
    if(axisA == 1)
        value = -1.0*temp2*dy+temp1*dy*dx*dx;
    elseif(axisA == 2)
        value = -1.0*temp2*dx+temp1*dy*dy*dx;
    elseif(axisA == 3)
        value = temp1*dx*dy*dz;
    end
elseif(valenceOrbitalA == 4 && valenceOrbitalB == 3)
    temp1 = 8.0*(gauMultAB*gauMultAB*sqrt(gauMultAB))...
        /(gauPlusAB*gauPlusAB*gauPlusAB);
    temp2 = 4.0*gauMultAB*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB);
    if(axisA == 1)
        value = -1.0*temp2*dz+temp1*dx*dx*dz;
    elseif(axisA == 2)
        value = temp1*dx*dy*dz;
    elseif(axisA == 3)
        value = -1.0*temp2*dx+temp1*dx*dz*dz;
    end
elseif(valenceOrbitalA == 3 && valenceOrbitalB == 4)
    temp1 = 8.0*(gauMultAB*gauMultAB*sqrt(gauMultAB))...
        /(gauPlusAB*gauPlusAB*gauPlusAB);
    temp2 = 4.0*gauMultAB*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB);
    if(axisA == 1)
        value = -1.0*temp2*dz+temp1*dz*dx*dx;
    elseif(axisA == 2)
        value = temp1*dx*dy*dz;
    elseif(axisA == 3)
        value = -1.0*temp2*dx+temp1*dz*dz*dx;
    end
elseif(valenceOrbitalA == 2 && valenceOrbitalB == 3)
    temp1 = 8.0*(gauMultAB*gauMultAB*sqrt(gauMultAB))...
        /(gauPlusAB*gauPlusAB*gauPlusAB);
    temp2 = 4.0*gauMultAB*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB);
    if(axisA == 1)
        value = temp1*dx*dy*dz;
    elseif(axisA == 2)
        value = -1.0*temp2*dz+temp1*dy*dy*dz;
    elseif(axisA == 3)
        value = -1.0*temp2*dy+temp1*dy*dz*dz;
    end
elseif(valenceOrbitalA == 3 && valenceOrbitalB == 2)
    temp1 = 8.0*(gauMultAB*gauMultAB*sqrt(gauMultAB))...
        /(gauPlusAB*gauPlusAB*gauPlusAB);
    temp2 = 4.0*gauMultAB*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB);
    if(axisA == 1)
        value = temp1*dx*dy*dz;
    elseif(axisA == 2)
        value = -1.0*temp2*dz+temp1*dz*dy*dy;
    elseif(axisA == 3)
        value = -1.0*temp2*dy+temp1*dz*dz*dy;
    end
elseif(valenceOrbitalA == 4 && valenceOrbitalB == 4)
    temp1 = 8.0*gauMultAB*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB);
    temp2 = gauMultAB/gauPlusAB;
    if(axisA == 1)
        value = -1.0*temp1*dx*(1.5-temp2*dx*dx);
    elseif(axisA == 2)
        value = -1.0*temp1*dy*(0.5-temp2*dx*dx);
    elseif(axisA == 3)
        value = -1.0*temp1*dz*(0.5-temp2*dx*dx);
    end
elseif(valenceOrbitalA == 2 && valenceOrbitalB == 2)
    temp1 = 8.0*gauMultAB*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB);
    temp2 = gauMultAB/gauPlusAB;
    if(axisA == 1)
        value = -1.0*temp1*dx*(0.5-temp2*dy*dy);
    elseif(axisA == 2)
        value = -1.0*temp1*dy*(1.5-temp2*dy*dy);
    elseif(axisA == 3)
        value = -1.0*temp1*dz*(0.5-temp2*dy*dy);
    end
elseif(valenceOrbitalA == 3 && valenceOrbitalB == 3)
    temp1 = 8.0*gauMultAB*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB);
    temp2 = gauMultAB/gauPlusAB;
    if(axisA == 1)
        value = -1.0*temp1*dx*(0.5-temp2*dz*dz);
    elseif(axisA == 2)
        value = -1.0*temp1*dy*(0.5-temp2*dz*dz);
    elseif(axisA == 3)
        value = -1.0*temp1*dz*(1.5-temp2*dz*dz);
    end
else
    throw(MException('Cndo2:GetGaussianOverlapAOs1stDerivative', 'Orbital type wrong.'));
end


overlapSASB = obj.GetGaussianOverlapAOsSASB(gaussianExponentA,...
    gaussianExponentB, rAB);
value = value * overlapSASB;
end
