#include "mex.h"
#include <math.h>



double Cpp_MndoGetSemiEmpiricalMultipoleInteraction(double* DAvec, double* rhoAvec,
        double* DBvec, double* rhoBvec,
        int multipoleA,
        int multipoleB,
        double rAB) {
    int dIndexA=0;
    switch(multipoleA){
        case 1:
            dIndexA = 0;
            break;
        case 8:
        case 9:
        case 10:
            dIndexA = 1;
            break;
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
            dIndexA = 2;
            break;
        default:
            mexErrMsgTxt("Cpp_MndoGetSemiEmpiricalMultipoleInteractionl: Multipole type wrong.");
    }
    
    int dIndexB=0;
    switch(multipoleB){
        case 1:
            dIndexB = 0;
            break;
        case 8:
        case 9:
        case 10:
            dIndexB = 1;
            break;
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
            dIndexB = 2;
            break;
        default:
            mexErrMsgTxt("Cpp_MndoGetSemiEmpiricalMultipoleInteractionl: Multipole type wrong.");
    }
    
    double DA = DAvec[dIndexA];
    double DB = DBvec[dIndexB];
    double rhoA = rhoAvec[dIndexA];
    double rhoB = rhoBvec[dIndexB];
    
    double value = 0.0;
    double a = rhoA + rhoB;
    
    // Eq. (52) in [DT_1977]
    if(multipoleA == 1 && multipoleB == 1){
        value = 1.0/sqrt(rAB*rAB + a*a);
    }
    // Eq. (53) in [DT_1977]
    else if(multipoleA == 1 && multipoleB == 10){
        double temp1 = ((rAB+DB)*(rAB+DB)) + (a*a);
        double temp2 = ((rAB-DB)*(rAB-DB)) + (a*a);
        value = 1.0/sqrt(temp1)/2.0 - 1.0/sqrt(temp2)/2.0;
    }
    else if(multipoleA == 10 && multipoleB == 1){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DBvec, rhoBvec, DAvec, rhoAvec, multipoleB, multipoleA, rAB);
        value *= -1.0;
    }
    // Eq. (54) in [DT_1977]
    else if(multipoleA == 1 && multipoleB == 2){
        double temp1 = (rAB*rAB) + (4.0*DB*DB) + (a*a);
        double temp2 = (rAB*rAB) + (a*a);
        value = 1.0/sqrt(temp1)/2.0 - 1.0/sqrt(temp2)/2.0;
    }
    else if(multipoleA == 2 && multipoleB == 1){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DBvec, rhoBvec, DAvec, rhoAvec, multipoleB, multipoleA, rAB);
    }
    else if(multipoleA == 1 && multipoleB == 3){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, multipoleA, 2, rAB);
    }
    else if(multipoleA == 3 && multipoleB == 1){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DBvec, rhoBvec, DAvec, rhoAvec, multipoleB, multipoleA, rAB);
    }
    // Eq. (55) in [DT_1977]
    else if(multipoleA == 1 && multipoleB == 4){
        double temp1 = ((rAB+2.0*DB)*(rAB+2.0*DB)) + (a*a);
        double temp2 = (rAB*rAB) + (a*a);
        double temp3 = ((rAB-2.0*DB)*(rAB-2.0*DB)) + (a*a);
        value = 1.0/sqrt(temp1)/4.0 - 1.0/sqrt(temp2)/2.0 + 1.0/sqrt(temp3)/4.0;
    }
    else if(multipoleA == 4 && multipoleB == 1){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DBvec, rhoBvec, DAvec, rhoAvec, multipoleB, multipoleA, rAB);
    }
    // Eq. (56) in [DT_1977]
    else if(multipoleA == 8 && multipoleB == 8){
        double temp1 = (rAB*rAB) + ((DA-DB)*(DA-DB)) + (a*a);
        double temp2 = (rAB*rAB) + ((DA+DB)*(DA+DB)) + (a*a);
        value = 1.0/sqrt(temp1)/2.0 - 1.0/sqrt(temp2)/2.0;
    }
    else if(multipoleA == 9 && multipoleB == 9){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 8, 8, rAB);
    }
    // Eq. (57) in [DT_1977]
    else if(multipoleA == 10 && multipoleB == 10){
        double temp1 = ((rAB+DA-DB)*(rAB+DA-DB)) + (a*a);
        double temp2 = ((rAB+DA+DB)*(rAB+DA+DB)) + (a*a);
        double temp3 = ((rAB-DA-DB)*(rAB-DA-DB)) + (a*a);
        double temp4 = ((rAB-DA+DB)*(rAB-DA+DB)) + (a*a);
        value = 1.0/sqrt(temp1)/4.0 - 1.0/sqrt(temp2)/4.0
                -1.0/sqrt(temp3)/4.0 + 1.0/sqrt(temp4)/4.0;
    }
    // Eq. (58) in [DT_1977]
    else if(multipoleA == 8 && multipoleB == 5){
        double temp1 = ((rAB-DB)*(rAB-DB)) + ((DA-DB)*(DA-DB)) + (a*a);
        double temp2 = ((rAB-DB)*(rAB-DB)) + ((DA+DB)*(DA+DB)) + (a*a);
        double temp3 = ((rAB+DB)*(rAB+DB)) + ((DA-DB)*(DA-DB)) + (a*a);
        double temp4 = ((rAB+DB)*(rAB+DB)) + ((DA+DB)*(DA+DB)) + (a*a);
        value =-1.0/sqrt(temp1)/4.0 + 1.0/sqrt(temp2)/4.0
                +1.0/sqrt(temp3)/4.0 - 1.0/sqrt(temp4)/4.0;
    }
    else if(multipoleA == 5 && multipoleB == 8){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DBvec, rhoBvec, DAvec, rhoAvec, multipoleB, multipoleA, rAB);
        value *= -1.0;
    }
    else if(multipoleA == 9 && multipoleB == 6){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 8, 5, rAB);
    }
    else if(multipoleA == 6 && multipoleB == 9){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DBvec, rhoBvec, DAvec, rhoAvec, multipoleB, multipoleA, rAB);
        value *= -1.0;
    }
    // Eq. (59) in [DT_1977]
    else if(multipoleA == 10 && multipoleB == 2){
        double temp1 = ((rAB+DA)*(rAB+DA)) + (4.0*DB*DB) + (a*a);
        double temp2 = ((rAB-DA)*(rAB-DA)) + (4.0*DB*DB) + (a*a);
        double temp3 = ((rAB+DA)*(rAB+DA)) + (a*a);
        double temp4 = ((rAB-DA)*(rAB-DA)) + (a*a);
        value =-1.0/sqrt(temp1)/4.0 + 1.0/sqrt(temp2)/4.0
                +1.0/sqrt(temp3)/4.0 - 1.0/sqrt(temp4)/4.0;
    }
    else if(multipoleA == 2 && multipoleB == 10){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DBvec, rhoBvec, DAvec, rhoAvec, multipoleB, multipoleA, rAB);
        value *= -1.0;
    }
    else if(multipoleA == 10 && multipoleB == 3){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 10, 2, rAB);
    }
    else if(multipoleA == 3 && multipoleB == 10){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DBvec, rhoBvec, DAvec, rhoAvec, multipoleB, multipoleA, rAB);
        value *= -1.0;
    }
    // Eq. (60) in [DT_1977]
    else if(multipoleA == 10 && multipoleB == 4){
        double temp1 = ((rAB+DA-2.0*DB)*(rAB+DA-2.0*DB)) + (a*a);
        double temp2 = ((rAB-DA-2.0*DB)*(rAB-DA-2.0*DB)) + (a*a);
        double temp3 = ((rAB+DA+2.0*DB)*(rAB+DA+2.0*DB)) + (a*a);
        double temp4 = ((rAB-DA+2.0*DB)*(rAB-DA+2.0*DB)) + (a*a);
        double temp5 = ((rAB+DA)*(rAB+DA)) + (a*a);
        double temp6 = ((rAB-DA)*(rAB-DA)) + (a*a);
        value =-1.0/sqrt(temp1)/8.0 + 1.0/sqrt(temp2)/8.0
                -1.0/sqrt(temp3)/8.0 + 1.0/sqrt(temp4)/8.0
                +1.0/sqrt(temp5)/4.0 - 1.0/sqrt(temp6)/4.0;
    }
    else if(multipoleA == 4 && multipoleB == 10){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DBvec, rhoBvec, DAvec, rhoAvec, multipoleB, multipoleA, rAB);
        value *= -1.0;
    }
    // Eq. (61) in [DT_1977]
    else if(multipoleA == 2 && multipoleB == 2){
        double temp1 = (rAB*rAB) + 4.0*((DA-DB)*(DA-DB)) + (a*a);
        double temp2 = (rAB*rAB) + 4.0*((DA+DB)*(DA+DB)) + (a*a);
        double temp3 = (rAB*rAB) + (4.0*DA*DA) + (a*a);
        double temp4 = (rAB*rAB) + (4.0*DB*DB) + (a*a);
        double temp5 = (rAB*rAB) + (a*a);
        value = 1.0/sqrt(temp1)/8.0 + 1.0/sqrt(temp2)/8.0
                -1.0/sqrt(temp3)/4.0 - 1.0/sqrt(temp4)/4.0
                +1.0/sqrt(temp5)/4.0;
    }
    else if(multipoleA == 3 && multipoleB == 3){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 2, 2, rAB);
    }
    // Eq. (62) in [DT_1977]
    else if(multipoleA == 2 && multipoleB == 3){
        double temp1 = (rAB*rAB) + (4.0*DA*DA) + (4.0*DB*DB)+ (a*a);
        double temp2 = (rAB*rAB) + (4.0*DA*DA) + (a*a);
        double temp3 = (rAB*rAB) + (4.0*DB*DB) + (a*a);
        double temp4 = (rAB*rAB) + (a*a);
        value = 1.0/sqrt(temp1)/4.0 - 1.0/sqrt(temp2)/4.0
                -1.0/sqrt(temp3)/4.0 + 1.0/sqrt(temp4)/4.0;
    }
    else if(multipoleA == 3 && multipoleB == 2){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DBvec, rhoBvec, DAvec, rhoAvec, multipoleB, multipoleA, rAB);
    }
    // Eq. (63) in [DT_1977]
    else if(multipoleA == 2 && multipoleB == 4){
        double temp1 = ((rAB-2.0*DB)*(rAB-2.0*DB)) + (4.0*DA*DA) + (a*a);
        double temp2 = ((rAB+2.0*DB)*(rAB+2.0*DB)) + (4.0*DA*DA) + (a*a);
        double temp3 = ((rAB-2.0*DB)*(rAB-2.0*DB)) + (a*a);
        double temp4 = ((rAB+2.0*DB)*(rAB+2.0*DB)) + (a*a);
        double temp5 = (rAB*rAB) + (4.0*DA*DA) + (a*a);
        double temp6 = (rAB*rAB) + (a*a);
        value = 1.0/sqrt(temp1)/8.0 + 1.0/sqrt(temp2)/8.0
                -1.0/sqrt(temp3)/8.0 - 1.0/sqrt(temp4)/8.0
                -1.0/sqrt(temp5)/4.0 + 1.0/sqrt(temp6)/4.0;
    }
    else if(multipoleA == 4 && multipoleB == 2){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DBvec, rhoBvec, DAvec, rhoAvec, multipoleB, multipoleA, rAB);
    }
    else if(multipoleA == 3 && multipoleB == 4){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 2, multipoleB, rAB);
    }
    else if(multipoleA == 4 && multipoleB == 3){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DBvec, rhoBvec, DAvec, rhoAvec, multipoleB, multipoleA, rAB);
    }
    // Eq. (64) in [DT_1977]
    else if(multipoleA == 4 && multipoleB == 4){
        double temp1 = ((rAB+2.0*DA-2.0*DB)*(rAB+2.0*DA-2.0*DB)) + (a*a);
        double temp2 = ((rAB+2.0*DA+2.0*DB)*(rAB+2.0*DA+2.0*DB)) + (a*a);
        double temp3 = ((rAB-2.0*DA-2.0*DB)*(rAB-2.0*DA-2.0*DB)) + (a*a);
        double temp4 = ((rAB-2.0*DA+2.0*DB)*(rAB-2.0*DA+2.0*DB)) + (a*a);
        double temp5 = ((rAB+2.0*DA)*(rAB+2.0*DA)) + (a*a);
        double temp6 = ((rAB-2.0*DA)*(rAB-2.0*DA)) + (a*a);
        double temp7 = ((rAB+2.0*DB)*(rAB+2.0*DB)) + (a*a);
        double temp8 = ((rAB-2.0*DB)*(rAB-2.0*DB)) + (a*a);
        double temp9 = (rAB*rAB) + (a*a);
        value = 1.0/sqrt(temp1)/16.0 + 1.0/sqrt(temp2)/16.0
                +1.0/sqrt(temp3)/16.0 + 1.0/sqrt(temp4)/16.0
                -1.0/sqrt(temp5)/8.0 - 1.0/sqrt(temp6)/8.0
                -1.0/sqrt(temp7)/8.0 - 1.0/sqrt(temp8)/8.0
                +1.0/sqrt(temp9)/4.0;
    }
    // Eq. (65) in [DT_1977]
    else if(multipoleA == 5 && multipoleB == 5){
        double temp1 = ((rAB+DA-DB)*(rAB+DA-DB)) + ((DA-DB)*(DA-DB)) + (a*a);
        double temp2 = ((rAB+DA-DB)*(rAB+DA-DB)) + ((DA+DB)*(DA+DB)) + (a*a);
        double temp3 = ((rAB+DA+DB)*(rAB+DA+DB)) + ((DA-DB)*(DA-DB)) + (a*a);
        double temp4 = ((rAB+DA+DB)*(rAB+DA+DB)) + ((DA+DB)*(DA+DB)) + (a*a);
        double temp5 = ((rAB-DA-DB)*(rAB-DA-DB)) + ((DA-DB)*(DA-DB)) + (a*a);
        double temp6 = ((rAB-DA-DB)*(rAB-DA-DB)) + ((DA+DB)*(DA+DB)) + (a*a);
        double temp7 = ((rAB-DA+DB)*(rAB-DA+DB)) + ((DA-DB)*(DA-DB)) + (a*a);
        double temp8 = ((rAB-DA+DB)*(rAB-DA+DB)) + ((DA+DB)*(DA+DB)) + (a*a);
        value = 1.0/sqrt(temp1)/8.0 - 1.0/sqrt(temp2)/8.0
                -1.0/sqrt(temp3)/8.0 + 1.0/sqrt(temp4)/8.0
                -1.0/sqrt(temp5)/8.0 + 1.0/sqrt(temp6)/8.0
                +1.0/sqrt(temp7)/8.0 - 1.0/sqrt(temp8)/8.0;
    }
    else if(multipoleA == 6 && multipoleB == 6){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 5, 5, rAB);
    }
    // Eq. (66) in [DT_1977]
    else if(multipoleA == 7 && multipoleB == 7){
        double temp1 = (rAB*rAB) + 2.0*((DA-DB)*(DA-DB)) + (a*a);
        double temp2 = (rAB*rAB) + 2.0*((DA+DB)*(DA+DB)) + (a*a);
        double temp3 = (rAB*rAB) + 2.0*(DA*DA) + 2.0*(DB*DB) + (a*a);
        value = 1.0/sqrt(temp1)/4.0 + 1.0/sqrt(temp2)/4.0
                -1.0/sqrt(temp3)/2.0;
    }
    else{
        mexErrMsgTxt("Cpp_MndoGetSemiEmpiricalMultipoleInteractionl: Multipole type wrong.");
    }
    return value;
}

double Cpp_MndoGetNddoRepulsionIntegral(double* DAvec, double* rhoAvec,
        int mu, int nu,
        double* DBvec, double* rhoBvec,
        int lambda, int sigma,
        double rAB) {
    double value = 0.0;
    
    if(mu == 1 && nu == 1 && lambda == 1 && sigma == 1){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
    }
    // (29) in [DT_1977]
    else if(mu == 1 && nu == 1 && lambda == 4 && sigma == 4){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 2, rAB);
        value = temp1 + temp2;
    }
    else if(mu == 1 && nu == 1 && lambda == 2 && sigma == 2){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 3, rAB);
        value = temp1 + temp2;
    }
    // (30) in [DT_1977]
    else if(mu == 1 && nu == 1 && lambda == 3 && sigma == 3){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 4, rAB);
        value = temp1 + temp2;
    }
    // (31) in [DT_1977]
    else if(mu == 4 && nu == 4 && lambda == 1 && sigma == 1){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 2, 1, rAB);
        value = temp1 + temp2;
    }
    else if(mu == 2 && nu == 2 && lambda == 1 && sigma == 1){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 3, 1, rAB);
        value = temp1 + temp2;
    }
    // (32) in [DT_1977]
    else if(mu == 3 && nu == 3 && lambda == 1 && sigma == 1){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 4, 1, rAB);
        value = temp1 + temp2;
    }
    // (33) in [DT_1977]
    else if(mu == 4 && nu == 4 && lambda == 4 && sigma == 4){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 2, rAB);
        double temp3 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 2, 1, rAB);
        double temp4 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 2, 2, rAB);
        value = temp1 + temp2 + temp3 + temp4;
    }
    else if(mu == 2 && nu == 2 && lambda == 2 && sigma == 2){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 3, rAB);
        double temp3 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 3, 1, rAB);
        double temp4 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 3, 3, rAB);
        value = temp1 + temp2 + temp3 + temp4;
    }
    // (34) in [DT_1977]
    else if(mu == 4 && nu == 4 && lambda == 2 && sigma == 2){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 3, rAB);
        double temp3 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 2, 1, rAB);
        double temp4 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 2, 3, rAB);
        value = temp1 + temp2 + temp3 + temp4;
    }
    else if(mu == 2 && nu == 2 && lambda == 4 && sigma == 4){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 2, rAB);
        double temp3 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 3, 1, rAB);
        double temp4 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 3, 2, rAB);
        value = temp1 + temp2 + temp3 + temp4;
    }
    // (35) in [DT_1977]
    else if(mu == 4 && nu == 4 && lambda == 3 && sigma == 3){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 4, rAB);
        double temp3 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 2, 1, rAB);
        double temp4 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 2, 4, rAB);
        value = temp1 + temp2 + temp3 + temp4;
    }
    else if(mu == 2 && nu == 2 && lambda == 3 && sigma == 3){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 4, rAB);
        double temp3 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 3, 1, rAB);
        double temp4 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 3, 4, rAB);
        value = temp1 + temp2 + temp3 + temp4;
    }
    // (36) in [DT_1977]
    else if(mu == 3 && nu == 3 && lambda == 4 && sigma == 4){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 2, rAB);
        double temp3 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 4, 1, rAB);
        double temp4 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 4, 2, rAB);
        value = temp1 + temp2 + temp3 + temp4;
    }
    else if(mu == 3 && nu == 3 && lambda == 2 && sigma == 2){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 3, rAB);
        double temp3 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 4, 1, rAB);
        double temp4 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 4, 3, rAB);
        value = temp1 + temp2 + temp3 + temp4;
    }
    // (37) in [DT_1977]
    else if(mu == 3 && nu == 3 && lambda == 3 && sigma == 3){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 4, rAB);
        double temp3 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 4, 1, rAB);
        double temp4 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 4, 4, rAB);
        value = temp1 + temp2 + temp3 + temp4;
    }
    // (38) in [DT_1977]
    else if(mu == 1 && nu == 3 && lambda == 1 && sigma == 1){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 10, 1, rAB);
        value = temp1;
    }
    else if(mu == 3 && nu == 1 && lambda == 1 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    // (39) in [DT_1977]
    else if(mu == 1 && nu == 3 && lambda == 4 && sigma == 4){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 10, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 10, 2, rAB);
        value = temp1 + temp2;
    }
    else if(mu == 3 && nu == 1 && lambda == 4 && sigma == 4){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    else if(mu == 1 && nu == 3 && lambda == 2 && sigma == 2){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 10, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 10, 3, rAB);
        value = temp1 + temp2;
    }
    else if(mu == 3 && nu == 1 && lambda == 2 && sigma == 2){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    // (40) in [DT_1977]
    else if(mu == 1 && nu == 3 && lambda == 3 && sigma == 3){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 10, 1, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 10, 4, rAB);
        value = temp1 + temp2;
    }
    else if(mu == 3 && nu == 1 && lambda == 3 && sigma == 3){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    // (41) in [DT_1977]
    else if(mu == 1 && nu == 1 && lambda == 1 && sigma == 3){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 10, rAB);
        value = temp1;
    }
    else if(mu == 1 && nu == 1 && lambda == 3 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    // (42) in [DT_1977]
    else if(mu == 4 && nu == 4 && lambda == 1 && sigma == 3){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 10, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 2, 10, rAB);
        value = temp1 + temp2;
    }
    else if(mu == 4 && nu == 4 && lambda == 3 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 2 && nu == 2 && lambda == 1 && sigma == 3){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 10, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 3, 10, rAB);
        value = temp1 + temp2;
    }
    else if(mu == 2 && nu == 2 && lambda == 3 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    // (43) in [DT_1977]
    else if(mu == 3 && nu == 3 && lambda == 1 && sigma == 3){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 1, 10, rAB);
        double temp2 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 4, 10, rAB);
        value = temp1 + temp2;
    }
    else if(mu == 3 && nu == 3 && lambda == 3 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    // (44) in [DT_1977]
    else if(mu == 1 && nu == 4 && lambda == 1 && sigma == 4){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 8, 8, rAB);
        value = temp1;
    }
    else if(mu == 4 && nu == 1 && lambda == 1 && sigma == 4){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    else if(mu == 1 && nu == 4 && lambda == 4 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 4 && nu == 1 && lambda == 4 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 1 && nu == 2 && lambda == 1 && sigma == 2){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 9, 9, rAB);
        value = temp1;
    }
    else if(mu == 2 && nu == 1 && lambda == 1 && sigma == 2){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    else if(mu == 1 && nu == 2 && lambda == 2 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 2 && nu == 1 && lambda == 2 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    // (45) in [DT_1977]
    else if(mu == 1 && nu == 3 && lambda == 1 && sigma == 3){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 10, 10, rAB);
        value = temp1;
    }
    else if(mu == 3 && nu == 1 && lambda == 1 && sigma == 3){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    else if(mu == 1 && nu == 3 && lambda == 3 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 3 && nu == 1 && lambda == 3 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    // (46) in [DT_1977]
    else if(mu == 1 && nu == 4 && lambda == 4 && sigma == 3){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 8, 5, rAB);
        value = temp1;
    }
    else if(mu == 4 && nu == 1 && lambda == 4 && sigma == 3){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    else if(mu == 1 && nu == 4 && lambda == 3 && sigma == 4){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 4 && nu == 1 && lambda == 3 && sigma == 4){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 1 && nu == 2 && lambda == 2 && sigma == 3){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 9, 6, rAB);
        value = temp1;
    }
    else if(mu == 2 && nu == 1 && lambda == 2 && sigma == 3){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    else if(mu == 1 && nu == 2 && lambda == 3 && sigma == 2){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 2 && nu == 1 && lambda == 3 && sigma == 2){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    // (47) in [DT_1977]
    else if(mu == 4 && nu == 3 && lambda == 1 && sigma == 4){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 5, 8, rAB);
        value = temp1;
    }
    else if(mu == 3 && nu == 4 && lambda == 1 && sigma == 4){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    else if(mu == 4 && nu == 3 && lambda == 4 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 3 && nu == 4 && lambda == 4 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 2 && nu == 3 && lambda == 1 && sigma == 2){
        double temp1 = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 6, 9, rAB);
        value = temp1;
    }
    else if(mu == 3 && nu == 2 && lambda == 1 && sigma == 2){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    else if(mu == 2 && nu == 3 && lambda == 2 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 3 && nu == 2 && lambda == 2 && sigma == 1){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    // (48) in [DT_1977]
    else if(mu == 4 && nu == 3 && lambda == 4 && sigma == 3){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 5, 5, rAB);
    }
    else if(mu == 3 && nu == 4 && lambda == 4 && sigma == 3){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    else if(mu == 4 && nu == 3 && lambda == 3 && sigma == 4){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 3 && nu == 4 && lambda == 3 && sigma == 4){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 2 && nu == 3 && lambda == 2 && sigma == 3){
        value = Cpp_MndoGetSemiEmpiricalMultipoleInteraction(DAvec, rhoAvec, DBvec, rhoBvec, 6, 6, rAB);
    }
    else if(mu == 3 && nu == 2 && lambda == 2 && sigma == 3){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    else if(mu == 2 && nu == 3 && lambda == 3 && sigma == 2){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 3 && nu == 2 && lambda == 3 && sigma == 2){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    // (49) in [DT_1977] and p19 in [MOPAC_1990]
    else if(mu == 4 && nu == 2 && lambda == 4 && sigma == 2){
        value = 0.5*(Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, mu, DBvec, rhoBvec, mu, mu, rAB)
        -Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, mu, DBvec, rhoBvec, nu, nu, rAB));
    }
    else if(mu == 2 && nu == 4 && lambda == 4 && sigma == 2){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, lambda, sigma, rAB);
    }
    else if(mu == 4 && nu == 2 && lambda == 2 && sigma == 4){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    else if(mu == 2 && nu == 4 && lambda == 2 && sigma == 4){
        value = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, nu, mu, DBvec, rhoBvec, sigma, lambda, rAB);
    }
    // d-orbitals
    else if(mu == 5 || mu == 6 || mu == 7 || mu == 8 || mu == 9 ||
            nu == 5 || nu == 6 || nu == 7 || nu == 8 || nu == 9 ||
            lambda == 5 || lambda == 6 || lambda == 7 || lambda  == 8 || lambda == 9 ||
            sigma == 5 || sigma == 6 || sigma == 7 || sigma  == 8 || sigma == 9){
        
        mexErrMsgTxt("Cpp_MndoGetNddoRepulsionIntegral: Orbital type wrong.");
    }
    else{
        value = 0.0;
    }
    return value;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nlhs != 1)
        mexErrMsgTxt("Cpp_MndoGetNddoRepulsionIntegral: 1 output expected.");
    if ( nrhs!=9 || 
            mxGetM(prhs[0])!=3 || mxGetN(prhs[0])!=1 || 
            mxGetM(prhs[1])!=3 || mxGetN(prhs[1])!=1 || 
            mxGetM(prhs[4])!=3 || mxGetN(prhs[4])!=1 || 
            mxGetM(prhs[5])!=3 || mxGetN(prhs[5])!=1 )
        mexErrMsgTxt("Cpp_MndoGetNddoRepulsionIntegral: 9 input expected.");
    double* DAvec = mxGetPr(prhs[0]);
    double* rhoAvec = mxGetPr(prhs[1]);
    int mu = (int)mxGetScalar(prhs[2]);
    int nu = (int)mxGetScalar(prhs[3]);
    double* DBvec = mxGetPr(prhs[4]);
    double* rhoBvec = mxGetPr(prhs[5]);
    int lambda = (int)mxGetScalar(prhs[6]);
    int sigma = (int)mxGetScalar(prhs[7]);
    double rAB = mxGetScalar(prhs[8]);
    
    plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL);
    *(mxGetPr(plhs[0])) = Cpp_MndoGetNddoRepulsionIntegral(DAvec, rhoAvec, mu, nu,
            DBvec, rhoBvec, lambda, sigma, rAB);
}

