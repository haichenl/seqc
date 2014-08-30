function value = GetGaussianOverlapAOs(~, ~, valenceOrbitalA, gaussianExponentA,...
    ~, valenceOrbitalB,gaussianExponentB,...
    dxyz, ~, overlapSASB)
dx = dxyz(1);
dy = dxyz(2);
dz = dxyz(3);
gauPlusAB = gaussianExponentA+gaussianExponentB;
gauMultAB = gaussianExponentA*gaussianExponentB;
if(valenceOrbitalA == OrbitalType.s && valenceOrbitalB == OrbitalType.s)
    value = 1.0;
elseif(valenceOrbitalA == OrbitalType.s && valenceOrbitalB == OrbitalType.px)
    value = 2.0*gaussianExponentA*sqrt(gaussianExponentB)*dx;
    value = value / gauPlusAB;
elseif(valenceOrbitalA == OrbitalType.s && valenceOrbitalB == OrbitalType.py)
    value = 2.0*gaussianExponentA*sqrt(gaussianExponentB)*dy;
    value = value / gauPlusAB;
elseif(valenceOrbitalA == OrbitalType.s && valenceOrbitalB == OrbitalType.pz)
    value = 2.0*gaussianExponentA*sqrt(gaussianExponentB)*dz;
    value = value / gauPlusAB;
elseif(valenceOrbitalA == OrbitalType.px && valenceOrbitalB == OrbitalType.s)
    value = -2.0*sqrt(gaussianExponentA)*gaussianExponentB*dx;
    value = value / gauPlusAB;
elseif(valenceOrbitalA == OrbitalType.py && valenceOrbitalB == OrbitalType.s)
    value = -2.0*sqrt(gaussianExponentA)*gaussianExponentB*dy;
    value = value / gauPlusAB;
elseif(valenceOrbitalA == OrbitalType.pz && valenceOrbitalB == OrbitalType.s)
    value = -2.0*sqrt(gaussianExponentA)*gaussianExponentB*dz;
    value = value / gauPlusAB;
elseif(valenceOrbitalA == OrbitalType.px && valenceOrbitalB == OrbitalType.px)
    temp = -1.0*(dx*dx)*gauMultAB;
    temp = temp / gauPlusAB;
    temp = temp + 0.5;
    value = 4.0*sqrt(gauMultAB);
    value = value / gauPlusAB;
    value = value * temp;
elseif(valenceOrbitalA == OrbitalType.px && valenceOrbitalB == OrbitalType.py)
    value = -4.0*gauMultAB*sqrt(gauMultAB);
    value = value / (gauPlusAB*gauPlusAB);
    value = value * dx*dy;
elseif(valenceOrbitalA == OrbitalType.px && valenceOrbitalB == OrbitalType.pz)
    value = -4.0*gauMultAB*sqrt(gauMultAB);
    value = value / (gauPlusAB*gauPlusAB);
    value = value * dx*dz;
elseif(valenceOrbitalA == OrbitalType.py && valenceOrbitalB == OrbitalType.px)
    value = -4.0*gauMultAB*sqrt(gauMultAB);
    value = value / (gauPlusAB*gauPlusAB);
    value = value * dy*dx;
elseif(valenceOrbitalA == OrbitalType.py && valenceOrbitalB == OrbitalType.py)
    temp = -1.0*(dy*dy)*gauMultAB;
    temp = temp / gauPlusAB;
    temp = temp + 0.5;
    value = 4.0*sqrt(gauMultAB);
    value = value / gauPlusAB;
    value = value * temp;
elseif(valenceOrbitalA == OrbitalType.py && valenceOrbitalB == OrbitalType.pz)
    value = -4.0*gauMultAB*sqrt(gauMultAB);
    value = value / (gauPlusAB*gauPlusAB);
    value = value * dy*dz;
elseif(valenceOrbitalA == OrbitalType.pz && valenceOrbitalB == OrbitalType.px)
    value = -4.0*gauMultAB*sqrt(gauMultAB);
    value = value / (gauPlusAB*gauPlusAB);
    value = value * dz*dx;
elseif(valenceOrbitalA == OrbitalType.pz && valenceOrbitalB == OrbitalType.py)
    value = -4.0*gauMultAB*sqrt(gauMultAB);
    value = value / (gauPlusAB*gauPlusAB);
    value = value * dz*dy;
elseif(valenceOrbitalA == OrbitalType.pz && valenceOrbitalB == OrbitalType.pz)
    temp = -1.0*(dz*dz)*gauMultAB;
    temp = temp / gauPlusAB;
    temp = temp + 0.5;
    value = 4.0*sqrt(gauMultAB);
    value = value / gauPlusAB;
    value = value * temp;
elseif(valenceOrbitalA == OrbitalType.dxy && valenceOrbitalB == OrbitalType.s)
    value = 4.0*gaussianExponentA;
    value = value / (gauPlusAB*gauPlusAB);
    value = value * gaussianExponentB*dx;
    value = value * gaussianExponentB*dy;
elseif(valenceOrbitalA == OrbitalType.dyz && valenceOrbitalB == OrbitalType.s)
    value = 4.0*gaussianExponentA;
    value = value / (gauPlusAB*gauPlusAB);
    value = value * gaussianExponentB*dy;
    value = value * gaussianExponentB*dz;
elseif(valenceOrbitalA == OrbitalType.dzx && valenceOrbitalB == OrbitalType.s)
    value = 4.0*gaussianExponentA;
    value = value / (gauPlusAB*gauPlusAB);
    value = value * gaussianExponentB*dz;
    value = value * gaussianExponentB*dx;
elseif(valenceOrbitalA == OrbitalType.s && valenceOrbitalB == OrbitalType.dxy)
    value = 4.0*gaussianExponentB;
    value = value / (gauPlusAB*gauPlusAB);
    value = value * gaussianExponentA*dx;
    value = value * gaussianExponentA*dy;
elseif(valenceOrbitalA == OrbitalType.s && valenceOrbitalB == OrbitalType.dyz)
    value = 4.0*gaussianExponentB;
    value = value / (gauPlusAB*gauPlusAB);
    value = value * gaussianExponentA*dy;
    value = value * gaussianExponentA*dz;
elseif(valenceOrbitalA == OrbitalType.s && valenceOrbitalB == OrbitalType.dzx)
    value = 4.0*gaussianExponentB;
    value = value / (gauPlusAB*gauPlusAB);
    value = value * gaussianExponentA*dz;
    value = value * gaussianExponentA*dx;
elseif(valenceOrbitalA == OrbitalType.dxy && valenceOrbitalB == OrbitalType.px)
    temp1 = -0.5*gaussianExponentB*dy;
    temp2 = gaussianExponentB*dx*gaussianExponentA*dx*gaussianExponentB*dy;
    temp2 = temp2 / gauPlusAB;
    value = temp1 + temp2;
    value = value * 8.0*gaussianExponentA*sqrt(gaussianExponentB);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dxy && valenceOrbitalB == OrbitalType.py)
    temp1 = -0.5*gaussianExponentB*dx;
    temp2 = gaussianExponentB*dy*gaussianExponentA*dy*gaussianExponentB*dx;
    temp2 = temp2 / gauPlusAB;
    value = temp1 + temp2;
    value = value * 8.0*gaussianExponentA*sqrt(gaussianExponentB);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dyz && valenceOrbitalB == OrbitalType.py)
    temp1 = -0.5*gaussianExponentB*dz;
    temp2 = gaussianExponentB*dy*gaussianExponentA*dy*gaussianExponentB*dz;
    temp2 = temp2 / gauPlusAB;
    value = temp1 + temp2;
    value = value * 8.0*gaussianExponentA*sqrt(gaussianExponentB);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dyz && valenceOrbitalB == OrbitalType.pz)
    temp1 = -0.5*gaussianExponentB*dy;
    temp2 = gaussianExponentB*dz*gaussianExponentA*dz*gaussianExponentB*dy;
    temp2 = temp2 / gauPlusAB;
    value = temp1 + temp2;
    value = value * 8.0*gaussianExponentA*sqrt(gaussianExponentB);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dzx && valenceOrbitalB == OrbitalType.pz)
    temp1 = -0.5*gaussianExponentB*dx;
    temp2 = gaussianExponentB*dz*gaussianExponentA*dz*gaussianExponentB*dx;
    temp2 = temp2 / gauPlusAB;
    value = temp1 + temp2;
    value = value * 8.0*gaussianExponentA*sqrt(gaussianExponentB);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dzx && valenceOrbitalB == OrbitalType.px)
    temp1 = -0.5*gaussianExponentB*dz;
    temp2 = gaussianExponentB*dx*gaussianExponentA*dx*gaussianExponentB*dz;
    temp2 = temp2 / gauPlusAB;
    value = temp1 + temp2;
    value = value * 8.0*gaussianExponentA*sqrt(gaussianExponentB);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.px && valenceOrbitalB == OrbitalType.dxy)
    temp1 = 0.5*gaussianExponentA*dy;
    temp2 = -1.0*gaussianExponentA*dx*gaussianExponentB*dx*gaussianExponentA*dy;
    temp2 = temp2 / gauPlusAB;
    value = temp1 + temp2;
    value = value * 8.0*gaussianExponentB*sqrt(gaussianExponentA);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.py && valenceOrbitalB == OrbitalType.dxy)
    temp1 = 0.5*gaussianExponentA*dx;
    temp2 = -1.0*gaussianExponentA*dy*gaussianExponentB*dy*gaussianExponentA*dx;
    temp2 = temp2 / gauPlusAB;
    value = temp1 + temp2;
    value = value * 8.0*gaussianExponentB*sqrt(gaussianExponentA);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.py && valenceOrbitalB == OrbitalType.dyz)
    temp1 = 0.5*gaussianExponentA*dz;
    temp2 = -1.0*gaussianExponentA*dy*gaussianExponentB*dy*gaussianExponentA*dz;
    temp2 = temp2 / gauPlusAB;
    value = temp1 + temp2;
    value = value * 8.0*gaussianExponentB*sqrt(gaussianExponentA);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.pz && valenceOrbitalB == OrbitalType.dyz)
    temp1 = 0.5*gaussianExponentA*dy;
    temp2 = -1.0*gaussianExponentA*dz*gaussianExponentB*dz*gaussianExponentA*dy;
    temp2 = temp2 / gauPlusAB;
    value = temp1 + temp2;
    value = value * 8.0*gaussianExponentB*sqrt(gaussianExponentA);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.pz && valenceOrbitalB == OrbitalType.dzx)
    temp1 = 0.5*gaussianExponentA*dx;
    temp2 = -1.0*gaussianExponentA*dz*gaussianExponentB*dz*gaussianExponentA*dx;
    temp2 = temp2 / gauPlusAB;
    value = temp1 + temp2;
    value = value * 8.0*gaussianExponentB*sqrt(gaussianExponentA);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.px && valenceOrbitalB == OrbitalType.dzx)
    temp1 = 0.5*gaussianExponentA*dz;
    temp2 = -1.0*gaussianExponentA*dx*gaussianExponentB*dx*gaussianExponentA*dz;
    temp2 = temp2 / gauPlusAB;
    value = temp1 + temp2;
    value = value * 8.0*gaussianExponentB*sqrt(gaussianExponentA);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dxy && valenceOrbitalB == OrbitalType.pz)
    value = 8.0*gaussianExponentA*sqrt(gaussianExponentB);
    value = value / (gauPlusAB*gauPlusAB*gauPlusAB);
    value = value * gaussianExponentB*dx*gaussianExponentB*dy*gaussianExponentA*dz;
elseif(valenceOrbitalA == OrbitalType.dyz && valenceOrbitalB == OrbitalType.px)
    value = 8.0*gaussianExponentA*sqrt(gaussianExponentB);
    value = value / (gauPlusAB*gauPlusAB*gauPlusAB);
    value = value * gaussianExponentB*dy*gaussianExponentB*dz*gaussianExponentA*dx;
elseif(valenceOrbitalA == OrbitalType.dzx && valenceOrbitalB == OrbitalType.py)
    value = 8.0*gaussianExponentA*sqrt(gaussianExponentB);
    value = value / (gauPlusAB*gauPlusAB*gauPlusAB);
    value = value * gaussianExponentB*dz*gaussianExponentB*dx*gaussianExponentA*dy;
elseif(valenceOrbitalA == OrbitalType.pz && valenceOrbitalB == OrbitalType.dxy)
    value = -8.0*gaussianExponentB*sqrt(gaussianExponentA);
    value = value / (gauPlusAB*gauPlusAB*gauPlusAB);
    value = value * gaussianExponentA*dx*gaussianExponentA*dy*gaussianExponentB*dz;
elseif(valenceOrbitalA == OrbitalType.px && valenceOrbitalB == OrbitalType.dyz)
    value = -8.0*gaussianExponentB*sqrt(gaussianExponentA);
    value = value / (gauPlusAB*gauPlusAB*gauPlusAB);
    value = value * gaussianExponentA*dy*gaussianExponentA*dz*gaussianExponentB*dx;
elseif(valenceOrbitalA == OrbitalType.py && valenceOrbitalB == OrbitalType.dzx)
    value = -8.0*gaussianExponentB*sqrt(gaussianExponentA);
    value = value / (gauPlusAB*gauPlusAB*gauPlusAB);
    value = value * gaussianExponentA*dz*gaussianExponentA*dx*gaussianExponentB*dy;
elseif(valenceOrbitalA == OrbitalType.dxxyy && valenceOrbitalB == OrbitalType.s)
    value = 2.0*gaussianExponentA;
    value = value / (gauPlusAB*gauPlusAB);
    value = value * ((gaussianExponentB*gaussianExponentB*dx*dx) - (gaussianExponentB*gaussianExponentB*dy*dy));
    
elseif(valenceOrbitalA == OrbitalType.s && valenceOrbitalB == OrbitalType.dxxyy)
    value = 2.0*gaussianExponentB;
    value = value / (gauPlusAB*gauPlusAB);
    value = value * ((gaussianExponentA*gaussianExponentA*dx*dx) - (gaussianExponentA*gaussianExponentA*dy*dy));
    
elseif(valenceOrbitalA == OrbitalType.dxxyy && valenceOrbitalB == OrbitalType.px)
    value = gaussianExponentB*dx;
    value = value - (gaussianExponentB*gaussianExponentB*dx*dx)*gaussianExponentA*dx/gauPlusAB;
    value = value + (gaussianExponentB*gaussianExponentB*dy*dy)*gaussianExponentA*dx/gauPlusAB;
    value = value * (-4.0*gaussianExponentA*sqrt(gaussianExponentB));
    value = value / (gauPlusAB*gauPlusAB);
    
elseif(valenceOrbitalA == OrbitalType.px && valenceOrbitalB == OrbitalType.dxxyy)
    value = gaussianExponentA*dx;
    value = value - (gaussianExponentA*gaussianExponentA*dx*dx)*gaussianExponentB*dx/gauPlusAB;
    value = value + (gaussianExponentA*gaussianExponentA*dy*dy)*gaussianExponentB*dx/gauPlusAB;
    value = value * 4.0*gaussianExponentB*sqrt(gaussianExponentA);
    value = value / (gauPlusAB*gauPlusAB);
    
elseif(valenceOrbitalA == OrbitalType.dxxyy && valenceOrbitalB == OrbitalType.py)
    value = gaussianExponentB*dy;
    value = value + (gaussianExponentB*gaussianExponentB*dx*dx)*gaussianExponentA*dy/gauPlusAB;
    value = value - (gaussianExponentB*gaussianExponentB*dy*dy)*gaussianExponentA*dy/gauPlusAB;
    value = value * 4.0*gaussianExponentA*sqrt(gaussianExponentB);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.py && valenceOrbitalB == OrbitalType.dxxyy)
    value = gaussianExponentA*dy;
    value = value + (gaussianExponentA*gaussianExponentA*dx*dx)*gaussianExponentB*dy/gauPlusAB;
    value = value - (gaussianExponentA*gaussianExponentA*dy*dy)*gaussianExponentB*dy/gauPlusAB;
    value = value * (-4.0*gaussianExponentB*sqrt(gaussianExponentA));
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dxxyy && valenceOrbitalB == OrbitalType.pz)
    value = (gaussianExponentB*gaussianExponentB*dx*dx) - (gaussianExponentB*gaussianExponentB*dy*dy);
    value = value * gaussianExponentA*dz;
    value = value * 4.0*gaussianExponentA*sqrt(gaussianExponentB);
    value = value / (gauPlusAB*gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.pz && valenceOrbitalB == OrbitalType.dxxyy)
    value = (gaussianExponentA*gaussianExponentA*dx*dx) - (gaussianExponentA*gaussianExponentA*dy*dy);
    value = value * gaussianExponentB*dz;
    value = value * (-4.0*gaussianExponentB*sqrt(gaussianExponentA));
    value = value / (gauPlusAB*gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dzz && valenceOrbitalB == OrbitalType.s)
    temp = 2.0*(gaussianExponentB*gaussianExponentB*dz*dz) ...
        -    (gaussianExponentB*gaussianExponentB*dx*dx) ...
        -    (gaussianExponentB*gaussianExponentB*dy*dy);
    value = 2.0*gaussianExponentA/sqrt(3.0);
    value = value / (gauPlusAB*gauPlusAB);
    value = value * temp;
elseif(valenceOrbitalA == OrbitalType.s && valenceOrbitalB == OrbitalType.dzz)
    temp = 2.0*(gaussianExponentA*gaussianExponentA*dz*dz) ...
        -    (gaussianExponentA*gaussianExponentA*dx*dx) ...
        -    (gaussianExponentA*gaussianExponentA*dy*dy);
    value = 2.0*gaussianExponentB/sqrt(3.0);
    value = value / (gauPlusAB*gauPlusAB);
    value = value * temp;
elseif(valenceOrbitalA == OrbitalType.dzz && valenceOrbitalB == OrbitalType.px)
    temp = gaussianExponentB*dx;
    temp = temp + 2.0*(gaussianExponentB*gaussianExponentB*dz*dz)*gaussianExponentA*dx/gauPlusAB;
    temp = temp -     (gaussianExponentB*gaussianExponentB*dx*dx)*gaussianExponentA*dx/gauPlusAB;
    temp = temp -     (gaussianExponentB*gaussianExponentB*dy*dy)*gaussianExponentA*dx/gauPlusAB;
    value = temp;
    value = value * 4.0*gaussianExponentA*sqrt(gaussianExponentB)/sqrt(3.0);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.px && valenceOrbitalB == OrbitalType.dzz)
    temp = gaussianExponentA*dx;
    temp = temp + 2.0*(gaussianExponentA*gaussianExponentA*dz*dz)*gaussianExponentB*dx/gauPlusAB;
    temp = temp -     (gaussianExponentA*gaussianExponentA*dx*dx)*gaussianExponentB*dx/gauPlusAB;
    temp = temp -     (gaussianExponentA*gaussianExponentA*dy*dy)*gaussianExponentB*dx/gauPlusAB;
    value = temp;
    value = value * (-4.0*gaussianExponentB*sqrt(gaussianExponentA)/sqrt(3.0));
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dzz && valenceOrbitalB == OrbitalType.py)
    temp = gaussianExponentB*dy;
    temp = temp + 2.0*(gaussianExponentB*gaussianExponentB*dz*dz)*gaussianExponentA*dy/gauPlusAB;
    temp = temp -     (gaussianExponentB*gaussianExponentB*dx*dx)*gaussianExponentA*dy/gauPlusAB;
    temp = temp -     (gaussianExponentB*gaussianExponentB*dy*dy)*gaussianExponentA*dy/gauPlusAB;
    value = temp;
    value = value * 4.0*gaussianExponentA*sqrt(gaussianExponentB)/sqrt(3.0);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.py && valenceOrbitalB == OrbitalType.dzz)
    temp = gaussianExponentA*dy;
    temp = temp + 2.0*(gaussianExponentA*gaussianExponentA*dz*dz)*gaussianExponentB*dy/gauPlusAB;
    temp = temp -     (gaussianExponentA*gaussianExponentA*dx*dx)*gaussianExponentB*dy/gauPlusAB;
    temp = temp -     (gaussianExponentA*gaussianExponentA*dy*dy)*gaussianExponentB*dy/gauPlusAB;
    value = temp;
    value = value * (-4.0*gaussianExponentB*sqrt(gaussianExponentA)/sqrt(3.0));
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dzz && valenceOrbitalB == OrbitalType.pz)
    temp = -2.0*gaussianExponentB*dz;
    temp = temp + 2.0*(gaussianExponentB*gaussianExponentB*dz*dz)*gaussianExponentA*dz/gauPlusAB;
    temp = temp -     (gaussianExponentB*gaussianExponentB*dx*dx)*gaussianExponentA*dz/gauPlusAB;
    temp = temp -     (gaussianExponentB*gaussianExponentB*dy*dy)*gaussianExponentA*dz/gauPlusAB;
    value = temp;
    value = value * 4.0*gaussianExponentA*sqrt(gaussianExponentB)/sqrt(3.0);
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.pz && valenceOrbitalB == OrbitalType.dzz)
    temp = -2.0*gaussianExponentA*dz;
    temp = temp + 2.0*(gaussianExponentA*gaussianExponentA*dz*dz)*gaussianExponentB*dz/gauPlusAB;
    temp = temp -     (gaussianExponentA*gaussianExponentA*dx*dx)*gaussianExponentB*dz/gauPlusAB;
    temp = temp -     (gaussianExponentA*gaussianExponentA*dy*dy)*gaussianExponentB*dz/gauPlusAB;
    value = temp;
    value = value * (-4.0*gaussianExponentB*sqrt(gaussianExponentA)/sqrt(3.0));
    value = value / (gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dxy && valenceOrbitalB == OrbitalType.dxy)
    temp = 0.25;
    temp = temp - 0.5*gaussianExponentB*dx*gaussianExponentA*dx/gauPlusAB;
    temp = temp - 0.5*gaussianExponentB*dy*gaussianExponentA*dy/gauPlusAB;
    temp = temp + gaussianExponentB*dx*gaussianExponentA*dx...
        *gaussianExponentB*dy*gaussianExponentA*dy...
        /(gauPlusAB*gauPlusAB);
    value = 16.0*temp*gauMultAB...
        /(gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dyz && valenceOrbitalB == OrbitalType.dyz)
    temp = 0.25;
    temp = temp - 0.5*gaussianExponentB*dy*gaussianExponentA*dy/gauPlusAB;
    temp = temp - 0.5*gaussianExponentB*dz*gaussianExponentA*dz/gauPlusAB;
    temp = temp + gaussianExponentB*dy*gaussianExponentA*dy...
        *gaussianExponentB*dz*gaussianExponentA*dz...
        /(gauPlusAB*gauPlusAB);
    value = 16.0*temp*gauMultAB...
        /(gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dzx && valenceOrbitalB == OrbitalType.dzx)
    temp = 0.25;
    temp = temp - 0.5*gaussianExponentB*dz*gaussianExponentA*dz/gauPlusAB;
    temp = temp - 0.5*gaussianExponentB*dx*gaussianExponentA*dx/gauPlusAB;
    temp = temp + gaussianExponentB*dz*gaussianExponentA*dz...
        *gaussianExponentB*dx*gaussianExponentA*dx...
        /(gauPlusAB*gauPlusAB);
    value = 16.0*temp*gauMultAB...
        /(gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dxxyy && valenceOrbitalB == OrbitalType.dxxyy)
    temp1 = 1.0;
    temp1 = temp1 -  2.0*gaussianExponentB*dx*gaussianExponentA*dx/gauPlusAB;
    temp1 = temp1 -  2.0*gaussianExponentB*dy*gaussianExponentA*dy/gauPlusAB;
    temp2 = gauMultAB*((dx*dx)-(dy*dy))...
        /gauPlusAB;
    temp1 = temp1 +  temp2*temp2;
    value = 4.0*temp1*gauMultAB...
        /(gauPlusAB*gauPlusAB);
elseif(valenceOrbitalA == OrbitalType.dzz && valenceOrbitalB == OrbitalType.dzz)
    temp = 3.0;
    temp = temp - gauMultAB...
        /gauPlusAB...
        *(8.0*(dz*dz)+2.0*(dx*dx)+2.0*(dy*dy));
    temp = temp + (gaussianExponentA*gauMultAB*gaussianExponentB)...
        *(4.0*(dz*dz*dz*dz)...
        +(dx*dx*dx*dx)...
        +(dy*dy*dy*dy)...
        -4.0*(dx*dx*dz*dz)...
        -4.0*(dy*dy*dz*dz)...
        +2.0*(dx*dx*dy*dy))...
        /(gauPlusAB*gauPlusAB);
    value = 4.0*temp*gauMultAB...
        /(3.0*gauPlusAB*gauPlusAB);
elseif((valenceOrbitalA == OrbitalType.dxy && valenceOrbitalB == OrbitalType.dyz) ||...
        (valenceOrbitalA == OrbitalType.dyz && valenceOrbitalB == OrbitalType.dxy))
    temp = 0.5;
    temp = temp - gauMultAB*(dy*dy)/gauPlusAB;
    value = -16.0*(gaussianExponentA*gauMultAB*gaussianExponentB)...
        *dx*dz*temp/(gauPlusAB*gauPlusAB*gauPlusAB);
elseif((valenceOrbitalA == OrbitalType.dyz && valenceOrbitalB == OrbitalType.dzx) ||...
        (valenceOrbitalA == OrbitalType.dzx && valenceOrbitalB == OrbitalType.dyz))
    temp = 0.5;
    temp = temp - gauMultAB*(dz*dz)/gauPlusAB;
    value = -16.0*(gaussianExponentA*gauMultAB*gaussianExponentB)...
        *dy*dx*temp/(gauPlusAB*gauPlusAB*gauPlusAB);
elseif((valenceOrbitalA == OrbitalType.dzx && valenceOrbitalB == OrbitalType.dxy) ||...
        (valenceOrbitalA == OrbitalType.dxy && valenceOrbitalB == OrbitalType.dzx))
    temp = 0.5;
    temp = temp - gauMultAB*(dx*dx)/gauPlusAB;
    value = -16.0*(gaussianExponentA*gauMultAB*gaussianExponentB)...
        *dz*dy*temp/(gauPlusAB*gauPlusAB*gauPlusAB);
elseif((valenceOrbitalA == OrbitalType.dxxyy && valenceOrbitalB == OrbitalType.dxy) ||...
        (valenceOrbitalA == OrbitalType.dxy && valenceOrbitalB == OrbitalType.dxxyy))
    temp = 2.0*gauMultAB;
    value = (temp*temp*temp)*(dy*(dx*dx*dx)-dx*(dy*dy*dy))...
        /(gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB);
elseif((valenceOrbitalA == OrbitalType.dxxyy && valenceOrbitalB == OrbitalType.dyz) ||...
        (valenceOrbitalA == OrbitalType.dyz && valenceOrbitalB == OrbitalType.dxxyy))
    temp = 2.0*gauMultAB;
    value = (temp*temp*temp)*(dy*dz*gauPlusAB...
        /gauMultAB...
        +((dx*dx)*dy*dz - (dy*dy*dy)*dz))...
        /(gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB);
elseif((valenceOrbitalA == OrbitalType.dxxyy && valenceOrbitalB == OrbitalType.dzx) ||...
        (valenceOrbitalA == OrbitalType.dzx && valenceOrbitalB == OrbitalType.dxxyy))
    temp = 2.0*gauMultAB;
    value = -1.0*(temp*temp*temp)*(dx*dz*gauPlusAB...
        /gauMultAB...
        +((dy*dy)*dx*dz - (dx*dx*dx)*dz))...
        /(gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB);
elseif((valenceOrbitalA == OrbitalType.dzz && valenceOrbitalB == OrbitalType.dxy) ||...
        (valenceOrbitalA == OrbitalType.dxy && valenceOrbitalB == OrbitalType.dzz))
    temp = 2.0*dx*dy*(dz*dz) - (dx*dx*dx)*dy - dx*(dy*dy*dy);
    temp = temp * gauMultAB/gauPlusAB;
    temp = temp + 2.0*dx*dy;
    value = 8.0*(gaussianExponentA*gauMultAB*gaussianExponentB)*temp;
    value = value / (sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB));
elseif((valenceOrbitalA == OrbitalType.dzz && valenceOrbitalB == OrbitalType.dyz) ||...
        (valenceOrbitalA == OrbitalType.dyz && valenceOrbitalB == OrbitalType.dzz))
    temp1 = -1.0*dy*dz;
    temp2 = 2.0*dy*(dz*dz*dz) - (dy*dy*dy)*dz - (dx*dx)*dy*dz;
    temp2 = temp2 * gauMultAB/gauPlusAB;
    temp1 = temp1 +  temp2;
    value = 8.0*(gaussianExponentA*gauMultAB*gaussianExponentB)*temp1;
    value = value / (sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB));
elseif((valenceOrbitalA == OrbitalType.dzz && valenceOrbitalB == OrbitalType.dzx) ||...
        (valenceOrbitalA == OrbitalType.dzx && valenceOrbitalB == OrbitalType.dzz))
    temp1 = -1.0*dx*dz;
    temp2 = 2.0*dx*(dz*dz*dz) - (dx*dx*dx)*dz - (dy*dy)*dx*dz;
    temp2 = temp2 * gauMultAB/gauPlusAB;
    temp1 = temp1 +  temp2;
    value = 8.0*(gaussianExponentA*gauMultAB*gaussianExponentB)*temp1;
    value = value / (sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB));
elseif((valenceOrbitalA == OrbitalType.dxxyy && valenceOrbitalB == OrbitalType.dzz) ||...
        (valenceOrbitalA == OrbitalType.dzz && valenceOrbitalB == OrbitalType.dxxyy))
    temp = 2.0*(dz*dz)-(dx*dx)-(dy*dy);
    temp = temp * gauMultAB/gauPlusAB;
    temp = temp + 2.0;
    value = 4.0*(gaussianExponentA*gauMultAB*gaussianExponentB);
    value = value / (sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB));
    value = value * ((dx*dx)-(dy*dy))*temp;
else
    throw(MException('Cndo2:GetGaussianOverlapAOs', 'Orbital type wrong.'));
end

value = value * overlapSASB;
end