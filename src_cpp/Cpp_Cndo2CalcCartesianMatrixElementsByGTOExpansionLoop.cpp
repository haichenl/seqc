#include "mex.h"
#include <math.h>


#define kShell 0
#define lShell 1
#define mShell 2
#define nShell 3
#define ShellType_end 4

#define sAzimuthal 0
#define pAzimuthal 1
#define dAzimuthal 2
#define AzimuthalType_end 3

#define s 0
#define py 1
#define pz 2
#define px 3
#define dxy 4
#define dyz 5
#define dzz 6
#define dzx 7
#define dxxyy 8
#define OrbitalType_end 9

#define STO1G 0
#define STO2G 1
#define STO3G 2
#define STO4G 3
#define STO5G 4
#define STO6G 5
#define STOnGType_end 6

class Uncopyable{
public:
protected:
   Uncopyable(){};
   ~Uncopyable(){};
private:
   Uncopyable(const Uncopyable&);
   Uncopyable& operator = (const Uncopyable&);
};

// GTOExpansionSTO is singleton
class GTOExpansionSTO: private Uncopyable{
public:
   static GTOExpansionSTO* GetInstance();
   static void DeleteInstance();
   double GetExponent(int stonG, 
                      int shellType, 
                      int orbitalType, 
                      int index) const;
   double GetCoefficient(int stonG, 
                         int shellType, 
                         int orbitalType, 
                         int index) const;

private:
   static GTOExpansionSTO* gTOExpansionSTO;
   GTOExpansionSTO();
   ~GTOExpansionSTO();
   
   void SetCoefficientsExponents();
   double exponents[STOnGType_end]
                   [ShellType_end]
                   [AzimuthalType_end]
                   [6];    
      //[N:expansion number][Shelltype][quasi orbital type:s, p, or d][expansion index]. 
      //This is alpha in (3) of [S_1970]. See Table I and II in [S_1970]
   double coefficients[STOnGType_end]
                      [ShellType_end]
                      [AzimuthalType_end]
                      [6]; 
      //[N:expansion number][Shelltype][quasi orbital type:s, p, or d][expansion index]. 
      //This is d in (3) of [S_1970]. See Table I and II in [S_1970]
};


GTOExpansionSTO* GTOExpansionSTO::gTOExpansionSTO = NULL;

GTOExpansionSTO::GTOExpansionSTO(){
   this->SetCoefficientsExponents();
}

GTOExpansionSTO::~GTOExpansionSTO(){
}

GTOExpansionSTO* GTOExpansionSTO::GetInstance(){
   if(gTOExpansionSTO == NULL){
      gTOExpansionSTO = new GTOExpansionSTO();
   }
   return gTOExpansionSTO;
}

void GTOExpansionSTO::DeleteInstance(){
   if(gTOExpansionSTO != NULL){
      delete gTOExpansionSTO; 
   }
   gTOExpansionSTO = NULL;
}

double GTOExpansionSTO::GetExponent(int stonG, int shellType, int orbitalType, int index) const{

   int azimuthalType;
   if(orbitalType == s){
      azimuthalType = sAzimuthal;
   }
   else if(orbitalType == px || orbitalType == py ||orbitalType == pz ){
      azimuthalType = pAzimuthal;
   }
   else if(orbitalType == dxy || orbitalType == dyz ||orbitalType == dzz || orbitalType == dzx ||orbitalType == dxxyy ){
      azimuthalType = dAzimuthal;
   }
   else{
      mexErrMsgTxt("GTOExpansionSTO::GetExponent: Orbital type wrong.");
   }

   return this->exponents[stonG][shellType][azimuthalType][index];
}

double GTOExpansionSTO::GetCoefficient(int stonG, int shellType, int orbitalType, int index) const{

   int azimuthalType;
   if(orbitalType == s){
      azimuthalType = sAzimuthal;
   }
   else if(orbitalType == px || orbitalType == py ||orbitalType == pz ){
      azimuthalType = pAzimuthal;
   }
   else if(orbitalType == dxy || orbitalType == dyz ||orbitalType == dzz || orbitalType == dzx ||orbitalType == dxxyy ){
      azimuthalType = dAzimuthal;
   }
   else{
      mexErrMsgTxt("GTOExpansionSTO::GetCoefficient: Orbital type wrong.");
   }

   return this->coefficients[stonG][shellType][azimuthalType][index];
}

//  see Table I and II in [S_1970]
void GTOExpansionSTO::SetCoefficientsExponents(){

   //STO-1G, 
   {   
      // 1s
      exponents[STO1G][kShell][sAzimuthal][0] = 2.709498091e-1;   coefficients[STO1G][kShell][sAzimuthal][0] = 1.0000;
      // 2s
      exponents[STO1G][lShell][sAzimuthal][0] = 1.012151084e-1;   coefficients[STO1G][lShell][sAzimuthal][0] = 1.0000;
      // 2p
      exponents[STO1G][lShell][pAzimuthal][0] = 1.759666885e-1;   coefficients[STO1G][lShell][pAzimuthal][0] = 1.0000;
      // 3s
      exponents[STO1G][mShell][sAzimuthal][0] = 5.296881757e-2;   coefficients[STO1G][mShell][sAzimuthal][0] = 1.0000;
      // 3p
      exponents[STO1G][mShell][pAzimuthal][0] = 9.113614253e-2;   coefficients[STO1G][mShell][pAzimuthal][0] = 1.0000;
      // 3d
      exponents[STO1G][mShell][dAzimuthal][0] = 1.302270363e-1;   coefficients[STO1G][mShell][dAzimuthal][0] = 1.0000;
      // 4s
      exponents[STO1G][nShell][sAzimuthal][0] = 3.264600274e-2;   coefficients[STO1G][nShell][sAzimuthal][0] = 1.0000;
      // 4p
      exponents[STO1G][nShell][pAzimuthal][0] = 5.578350235e-2;   coefficients[STO1G][nShell][pAzimuthal][0] = 1.0000;
      // 4d
      exponents[STO1G][nShell][dAzimuthal][0] = 7.941656339e-2;   coefficients[STO1G][nShell][dAzimuthal][0] = 1.0000;
   }

   //STO-2G, 
   {
      // 1s
      exponents[STO2G][kShell][sAzimuthal][0] = 8.518186635e-1;   coefficients[STO2G][kShell][sAzimuthal][0] = 4.301284983e-1; 
      exponents[STO2G][kShell][sAzimuthal][1] = 1.516232927e-1;   coefficients[STO2G][kShell][sAzimuthal][1] = 6.789135305e-1; 
      // 2s
      exponents[STO2G][lShell][sAzimuthal][0] = 1.292278611e-1;   coefficients[STO2G][lShell][sAzimuthal][0] = 7.470867124e-1; 
      exponents[STO2G][lShell][sAzimuthal][1] = 4.908584205e-2;   coefficients[STO2G][lShell][sAzimuthal][1] = 2.855980556e-1; 
      // 2p
      exponents[STO2G][lShell][pAzimuthal][0] = 4.323908358e-1;   coefficients[STO2G][lShell][pAzimuthal][0] = 4.522627513e-1; 
      exponents[STO2G][lShell][pAzimuthal][1] = 1.069439065e-1;   coefficients[STO2G][lShell][pAzimuthal][1] = 6.713122642e-1; 
      // 3s
      exponents[STO2G][mShell][sAzimuthal][0] = 6.694095822e-1;   coefficients[STO2G][mShell][sAzimuthal][0] =-1.529645716e-1; 
      exponents[STO2G][mShell][sAzimuthal][1] = 5.837135094e-2;   coefficients[STO2G][mShell][sAzimuthal][1] = 1.051370110; 
      // 3p
      exponents[STO2G][mShell][pAzimuthal][0] = 1.458620964e-1;   coefficients[STO2G][mShell][pAzimuthal][0] = 5.349653144e-1; 
      exponents[STO2G][mShell][pAzimuthal][1] = 5.664210742e-2;   coefficients[STO2G][mShell][pAzimuthal][1] = 5.299607212e-1; 
      // 3d
      exponents[STO2G][mShell][dAzimuthal][0] = 2.777427345e-1;   coefficients[STO2G][mShell][dAzimuthal][0] = 4.666137923e-1; 
      exponents[STO2G][mShell][dAzimuthal][1] = 8.336507714e-2;   coefficients[STO2G][mShell][dAzimuthal][1] = 6.644706516e-1; 
      // 4s
      exponents[STO2G][nShell][sAzimuthal][0] = 2.441785453e-1;   coefficients[STO2G][nShell][sAzimuthal][0] =-3.046656896e-1; 
      exponents[STO2G][nShell][sAzimuthal][1] = 4.051097664e-2;   coefficients[STO2G][nShell][sAzimuthal][1] = 1.146877294e00; 
      // 4p                                                               
      exponents[STO2G][nShell][pAzimuthal][0] = 6.190052680e-2;   coefficients[STO2G][nShell][pAzimuthal][0] = 8.743116767e-1; 
      exponents[STO2G][nShell][pAzimuthal][1] = 2.648418407e-2;   coefficients[STO2G][nShell][pAzimuthal][1] = 1.513640107e-1; 
      // 4d                                                               
      exponents[STO2G][nShell][dAzimuthal][0] = 1.330958892e-1;   coefficients[STO2G][nShell][dAzimuthal][0] = 4.932764167e-1; 
      exponents[STO2G][nShell][dAzimuthal][1] = 5.272119659e-2;   coefficients[STO2G][nShell][dAzimuthal][1] = 5.918727866e-1; 
   }

   //STO-3G, 
   {
      // 1s
      exponents[STO3G][kShell][sAzimuthal][0] = 2.227660584e00;   coefficients[STO3G][kShell][sAzimuthal][0] = 1.543289673e-1; 
      exponents[STO3G][kShell][sAzimuthal][1] = 4.057711562e-1;   coefficients[STO3G][kShell][sAzimuthal][1] = 5.353281423e-1; 
      exponents[STO3G][kShell][sAzimuthal][2] = 1.098175104e-1;   coefficients[STO3G][kShell][sAzimuthal][2] = 4.446345422e-1; 
      // 2s
      exponents[STO3G][lShell][sAzimuthal][0] = 2.581578398e00;   coefficients[STO3G][lShell][sAzimuthal][0] =-5.994474934e-2; 
      exponents[STO3G][lShell][sAzimuthal][1] = 1.567622104e-1;   coefficients[STO3G][lShell][sAzimuthal][1] = 5.960385398e-1; 
      exponents[STO3G][lShell][sAzimuthal][2] = 6.018332272e-2;   coefficients[STO3G][lShell][sAzimuthal][2] = 4.581786291e-1; 
      // 2p
      exponents[STO3G][lShell][pAzimuthal][0] = 9.192379002e-1;   coefficients[STO3G][lShell][pAzimuthal][0] = 1.623948553e-1; 
      exponents[STO3G][lShell][pAzimuthal][1] = 2.359194503e-1;   coefficients[STO3G][lShell][pAzimuthal][1] = 5.661708862e-1; 
      exponents[STO3G][lShell][pAzimuthal][2] = 8.009805746e-2;   coefficients[STO3G][lShell][pAzimuthal][2] = 4.223071752e-1; 
      // 3s
      exponents[STO3G][mShell][sAzimuthal][0] = 5.641487709e-1;   coefficients[STO3G][mShell][sAzimuthal][0] =-1.782577972e-1; 
      exponents[STO3G][mShell][sAzimuthal][1] = 6.924421391e-2;   coefficients[STO3G][mShell][sAzimuthal][1] = 8.612761663e-1; 
      exponents[STO3G][mShell][sAzimuthal][2] = 3.269529097e-2;   coefficients[STO3G][mShell][sAzimuthal][2] = 2.261841969e-1; 
      // 3p
      exponents[STO3G][mShell][pAzimuthal][0] = 2.692880368e00;   coefficients[STO3G][mShell][pAzimuthal][0] =-1.061945788e-2; 
      exponents[STO3G][mShell][pAzimuthal][1] = 1.489359592e-1;   coefficients[STO3G][mShell][pAzimuthal][1] = 5.218564264e-1; 
      exponents[STO3G][mShell][pAzimuthal][2] = 5.739585040e-2;   coefficients[STO3G][mShell][pAzimuthal][2] = 5.450015143e-1; 
      // 3d
      exponents[STO3G][mShell][dAzimuthal][0] = 5.229112225e-1;   coefficients[STO3G][mShell][dAzimuthal][0] = 1.686596060e-1; 
      exponents[STO3G][mShell][dAzimuthal][1] = 1.639595876e-1;   coefficients[STO3G][mShell][dAzimuthal][1] = 5.847984817e-1; 
      exponents[STO3G][mShell][dAzimuthal][2] = 6.386630021e-2;   coefficients[STO3G][mShell][dAzimuthal][2] = 4.056779523e-1; 
      // 4s
      exponents[STO3G][nShell][sAzimuthal][0] = 2.267938753e-1;   coefficients[STO3G][nShell][sAzimuthal][0] =-3.349048323e-1; 
      exponents[STO3G][nShell][sAzimuthal][1] = 4.448178019e-2;   coefficients[STO3G][nShell][sAzimuthal][1] = 1.056744667e00; 
      exponents[STO3G][nShell][sAzimuthal][2] = 2.195294664e-2;   coefficients[STO3G][nShell][sAzimuthal][2] = 1.256661680e-1; 
      // 4p                                                               
      exponents[STO3G][nShell][pAzimuthal][0] = 4.859692220e-1;   coefficients[STO3G][nShell][pAzimuthal][0] =-6.147823411e-2; 
      exponents[STO3G][nShell][pAzimuthal][1] = 7.430216918e-2;   coefficients[STO3G][nShell][pAzimuthal][1] = 6.604172234e-1; 
      exponents[STO3G][nShell][pAzimuthal][2] = 3.653340923e-2;   coefficients[STO3G][nShell][pAzimuthal][2] = 3.932639495e-1; 
      // 4d                                                               
      exponents[STO3G][nShell][dAzimuthal][0] = 1.777717219e-1;   coefficients[STO3G][nShell][dAzimuthal][0] = 2.308552718e-1; 
      exponents[STO3G][nShell][dAzimuthal][1] = 8.040647350e-2;   coefficients[STO3G][nShell][dAzimuthal][1] = 6.042409177e-1; 
      exponents[STO3G][nShell][dAzimuthal][2] = 3.949855551e-2;   coefficients[STO3G][nShell][dAzimuthal][2] = 2.595768926e-1; 
   }

   //STO-4G, 
   {
      // 1s
      exponents[STO4G][kShell][sAzimuthal][0] = 5.216844534e00;   coefficients[STO4G][kShell][sAzimuthal][0] = 5.675242080e-2; 
      exponents[STO4G][kShell][sAzimuthal][1] = 9.546182760e-1;   coefficients[STO4G][kShell][sAzimuthal][1] = 2.601413550e-1; 
      exponents[STO4G][kShell][sAzimuthal][2] = 2.652034102e-1;   coefficients[STO4G][kShell][sAzimuthal][2] = 5.328461143e-1; 
      exponents[STO4G][kShell][sAzimuthal][3] = 8.801862774e-2;   coefficients[STO4G][kShell][sAzimuthal][3] = 2.916254405e-1; 
      // 2s
      exponents[STO4G][lShell][sAzimuthal][0] = 1.161525551e01;   coefficients[STO4G][lShell][sAzimuthal][0] =-1.198411747e-2; 
      exponents[STO4G][lShell][sAzimuthal][1] = 2.000243111e00;   coefficients[STO4G][lShell][sAzimuthal][1] =-5.472052539e-2; 
      exponents[STO4G][lShell][sAzimuthal][2] = 1.607280687e-1;   coefficients[STO4G][lShell][sAzimuthal][2] = 5.805587176e-1; 
      exponents[STO4G][lShell][sAzimuthal][3] = 6.125744532e-2;   coefficients[STO4G][lShell][sAzimuthal][3] = 4.770079976e-1; 
      // 2p
      exponents[STO4G][lShell][pAzimuthal][0] = 1.798260992e00;   coefficients[STO4G][lShell][pAzimuthal][0] = 5.713170255e-2; 
      exponents[STO4G][lShell][pAzimuthal][1] = 4.662622228e-1;   coefficients[STO4G][lShell][pAzimuthal][1] = 2.857455515e-1; 
      exponents[STO4G][lShell][pAzimuthal][2] = 1.643718620e-1;   coefficients[STO4G][lShell][pAzimuthal][2] = 5.517873105e-1; 
      exponents[STO4G][lShell][pAzimuthal][3] = 6.543927065e-2;   coefficients[STO4G][lShell][pAzimuthal][3] = 2.632314924e-1; 
      // 3s
      exponents[STO4G][mShell][sAzimuthal][0] = 1.513265591e00;   coefficients[STO4G][mShell][sAzimuthal][0] =-3.295496352e-2; 
      exponents[STO4G][mShell][sAzimuthal][1] = 4.262497508e-1;   coefficients[STO4G][mShell][sAzimuthal][1] =-1.724516959e-1; 
      exponents[STO4G][mShell][sAzimuthal][2] = 7.643320863e-2;   coefficients[STO4G][mShell][sAzimuthal][2] = 7.518511194e-1; 
      exponents[STO4G][mShell][sAzimuthal][3] = 3.760545063e-2;   coefficients[STO4G][mShell][sAzimuthal][3] = 3.589627317e-1; 
      // 3p
      exponents[STO4G][mShell][pAzimuthal][0] = 1.853180239e00;   coefficients[STO4G][mShell][pAzimuthal][0] =-1.434249391e-2; 
      exponents[STO4G][mShell][pAzimuthal][1] = 1.915075719e-1;   coefficients[STO4G][mShell][pAzimuthal][1] = 2.755177589e-1; 
      exponents[STO4G][mShell][pAzimuthal][2] = 8.655487938e-2;   coefficients[STO4G][mShell][pAzimuthal][2] = 5.846750879e-1; 
      exponents[STO4G][mShell][pAzimuthal][3] = 4.184253862e-2;   coefficients[STO4G][mShell][pAzimuthal][3] = 2.144986514e-1; 
      // 3d
      exponents[STO4G][mShell][dAzimuthal][0] = 9.185846715e-1;   coefficients[STO4G][mShell][dAzimuthal][0] = 5.799057705e-2; 
      exponents[STO4G][mShell][dAzimuthal][1] = 2.920461109e-1;   coefficients[STO4G][mShell][dAzimuthal][1] = 3.045581349e-1; 
      exponents[STO4G][mShell][dAzimuthal][2] = 1.187568890e-1;   coefficients[STO4G][mShell][dAzimuthal][2] = 5.601358038e-1; 
      exponents[STO4G][mShell][dAzimuthal][3] = 5.286755896e-2;   coefficients[STO4G][mShell][dAzimuthal][3] = 2.432423313e-1; 
      // 4s
      exponents[STO4G][nShell][sAzimuthal][0] = 3.242212833e-1;   coefficients[STO4G][nShell][sAzimuthal][0] =-1.120682822e-1; 
      exponents[STO4G][nShell][sAzimuthal][1] = 1.663217177e-1;   coefficients[STO4G][nShell][sAzimuthal][1] =-2.845426863e-1; 
      exponents[STO4G][nShell][sAzimuthal][2] = 5.081097451e-2;   coefficients[STO4G][nShell][sAzimuthal][2] = 8.909873788e-1; 
      exponents[STO4G][nShell][sAzimuthal][3] = 2.829066600e-2;   coefficients[STO4G][nShell][sAzimuthal][3] = 3.517811205e-1; 
      // 4p                                                               
      exponents[STO4G][nShell][pAzimuthal][0] = 1.492607880e00;   coefficients[STO4G][nShell][pAzimuthal][0] =-6.035216774e-3; 
      exponents[STO4G][nShell][pAzimuthal][1] = 4.327619272e-1;   coefficients[STO4G][nShell][pAzimuthal][1] =-6.013310874e-2; 
      exponents[STO4G][nShell][pAzimuthal][2] = 7.553156064e-2;   coefficients[STO4G][nShell][pAzimuthal][2] = 6.451518200e-1; 
      exponents[STO4G][nShell][pAzimuthal][3] = 3.706272183e-2;   coefficients[STO4G][nShell][pAzimuthal][3] = 4.117923820e-1; 
      // 4d                                                               
      exponents[STO4G][nShell][dAzimuthal][0] = 1.995825422e00;   coefficients[STO4G][nShell][dAzimuthal][0] =-2.816702620e-3; 
      exponents[STO4G][nShell][dAzimuthal][1] = 1.823461280e-1;   coefficients[STO4G][nShell][dAzimuthal][1] = 2.177095871e-1; 
      exponents[STO4G][nShell][dAzimuthal][2] = 8.197240896e-2;   coefficients[STO4G][nShell][dAzimuthal][2] = 6.058047348e-1; 
      exponents[STO4G][nShell][dAzimuthal][3] = 4.000634951e-2;   coefficients[STO4G][nShell][dAzimuthal][3] = 2.717811257e-1; 
   }

   //STO-5G, 
   {
      // 1s
      exponents[STO5G][kShell][sAzimuthal][0] = 1.130563696e01;   coefficients[STO5G][kShell][sAzimuthal][0] = 2.214055312e-2; 
      exponents[STO5G][kShell][sAzimuthal][1] = 2.071728178e00;   coefficients[STO5G][kShell][sAzimuthal][1] = 1.135411520e-1; 
      exponents[STO5G][kShell][sAzimuthal][2] = 5.786484833e-1;   coefficients[STO5G][kShell][sAzimuthal][2] = 3.318161484e-1; 
      exponents[STO5G][kShell][sAzimuthal][3] = 1.975724573e-1;   coefficients[STO5G][kShell][sAzimuthal][3] = 4.825700713e-1; 
      exponents[STO5G][kShell][sAzimuthal][4] = 7.445271746e-2;   coefficients[STO5G][kShell][sAzimuthal][4] = 1.935721966e-1; 
      // 2s
      exponents[STO5G][lShell][sAzimuthal][0] = 8.984956862e00;   coefficients[STO5G][lShell][sAzimuthal][0] =-1.596349096e-2; 
      exponents[STO5G][lShell][sAzimuthal][1] = 1.673710636e00;   coefficients[STO5G][lShell][sAzimuthal][1] =-5.685884883e-2; 
      exponents[STO5G][lShell][sAzimuthal][2] = 1.944726668e-1;   coefficients[STO5G][lShell][sAzimuthal][2] = 3.698265599e-1; 
      exponents[STO5G][lShell][sAzimuthal][3] = 8.806345634e-2;   coefficients[STO5G][lShell][sAzimuthal][3] = 5.480512593e-1; 
      exponents[STO5G][lShell][sAzimuthal][4] = 4.249068522e-2;   coefficients[STO5G][lShell][sAzimuthal][4] = 1.472634893e-1; 
      // 2p
      exponents[STO5G][lShell][pAzimuthal][0] = 3.320386533e00;   coefficients[STO5G][lShell][pAzimuthal][0] = 2.079051117e-2; 
      exponents[STO5G][lShell][pAzimuthal][1] = 8.643257633e-1;   coefficients[STO5G][lShell][pAzimuthal][1] = 1.235472099e-1; 
      exponents[STO5G][lShell][pAzimuthal][2] = 3.079819284e-1;   coefficients[STO5G][lShell][pAzimuthal][2] = 3.667738986e-1; 
      exponents[STO5G][lShell][pAzimuthal][3] = 1.273309895e-1;   coefficients[STO5G][lShell][pAzimuthal][3] = 4.834930290e-1; 
      exponents[STO5G][lShell][pAzimuthal][4] = 5.606243164e-2;   coefficients[STO5G][lShell][pAzimuthal][4] = 1.653444074e-1; 
      // 3s
      exponents[STO5G][mShell][sAzimuthal][0] = 4.275877914e00;   coefficients[STO5G][mShell][sAzimuthal][0] =-3.920358850e-3; 
      exponents[STO5G][mShell][sAzimuthal][1] = 1.132409433e00;   coefficients[STO5G][mShell][sAzimuthal][1] =-4.168430506e-2; 
      exponents[STO5G][mShell][sAzimuthal][2] = 4.016256968e-1;   coefficients[STO5G][mShell][sAzimuthal][2] =-1.637440990e-1; 
      exponents[STO5G][mShell][sAzimuthal][3] = 7.732370620e-2;   coefficients[STO5G][mShell][sAzimuthal][3] = 7.419373723e-1; 
      exponents[STO5G][mShell][sAzimuthal][4] = 3.800708627e-2;   coefficients[STO5G][mShell][sAzimuthal][4] = 3.724364929e-1; 
      // 3p
      exponents[STO5G][mShell][pAzimuthal][0] = 6.466803859e00;   coefficients[STO5G][mShell][pAzimuthal][0] =-2.329023747e-3; 
      exponents[STO5G][mShell][pAzimuthal][1] = 1.555914802e00;   coefficients[STO5G][mShell][pAzimuthal][1] =-1.357395221e-2; 
      exponents[STO5G][mShell][pAzimuthal][2] = 1.955925255e-1;   coefficients[STO5G][mShell][pAzimuthal][2] = 2.632185383e-1; 
      exponents[STO5G][mShell][pAzimuthal][3] = 8.809647701e-2;   coefficients[STO5G][mShell][pAzimuthal][3] = 5.880427024e-1; 
      exponents[STO5G][mShell][pAzimuthal][4] = 4.234835707e-2;   coefficients[STO5G][mShell][pAzimuthal][4] = 2.242794445e-1; 
      // 3d
      exponents[STO5G][mShell][dAzimuthal][0] = 1.539033958e00;   coefficients[STO5G][mShell][dAzimuthal][0] = 2.020869128e-2; 
      exponents[STO5G][mShell][dAzimuthal][1] = 4.922090297e-1;   coefficients[STO5G][mShell][dAzimuthal][1] = 1.321157923e-1; 
      exponents[STO5G][mShell][dAzimuthal][2] = 2.029756928e-1;   coefficients[STO5G][mShell][dAzimuthal][2] = 3.911240346e-1; 
      exponents[STO5G][mShell][dAzimuthal][3] = 9.424112917e-2;   coefficients[STO5G][mShell][dAzimuthal][3] = 4.779609701e-1; 
      exponents[STO5G][mShell][dAzimuthal][4] = 4.569058269e-2;   coefficients[STO5G][mShell][dAzimuthal][4] = 1.463662294e-1; 
      // 4s
      exponents[STO5G][nShell][sAzimuthal][0] = 2.980263783e00;   coefficients[STO5G][nShell][sAzimuthal][0] = 1.513948997e-3; 
      exponents[STO5G][nShell][sAzimuthal][1] = 3.792228833e-1;   coefficients[STO5G][nShell][sAzimuthal][1] =-7.316801518e-2; 
      exponents[STO5G][nShell][sAzimuthal][2] = 1.789717224e-1;   coefficients[STO5G][nShell][sAzimuthal][2] =-3.143703799e-1; 
      exponents[STO5G][nShell][sAzimuthal][3] = 5.002110360e-2;   coefficients[STO5G][nShell][sAzimuthal][3] = 9.032615169e-1; 
      exponents[STO5G][nShell][sAzimuthal][4] = 2.789361681e-2;   coefficients[STO5G][nShell][sAzimuthal][4] = 3.294210848e-1; 
      // 4p                                                               
      exponents[STO5G][nShell][pAzimuthal][0] = 1.091977298e00;   coefficients[STO5G][nShell][pAzimuthal][0] =-1.143929558e-2; 
      exponents[STO5G][nShell][pAzimuthal][1] = 3.719985051e-1;   coefficients[STO5G][nShell][pAzimuthal][1] =-6.322651538e-2; 
      exponents[STO5G][nShell][pAzimuthal][2] = 8.590019352e-2;   coefficients[STO5G][nShell][pAzimuthal][2] = 4.398907721e-1; 
      exponents[STO5G][nShell][pAzimuthal][3] = 4.786503860e-2;   coefficients[STO5G][nShell][pAzimuthal][3] = 5.245859166e-1; 
      exponents[STO5G][nShell][pAzimuthal][4] = 2.730479990e-2;   coefficients[STO5G][nShell][pAzimuthal][4] = 1.017072253e-1; 
      // 4d                                                               
      exponents[STO5G][nShell][dAzimuthal][0] = 1.522122079e00;   coefficients[STO5G][nShell][dAzimuthal][0] =-3.673711876e-3; 
      exponents[STO5G][nShell][dAzimuthal][1] = 2.173041823e-1;   coefficients[STO5G][nShell][dAzimuthal][1] = 1.167122499e-1; 
      exponents[STO5G][nShell][dAzimuthal][2] = 1.084876577e-1;   coefficients[STO5G][nShell][dAzimuthal][2] = 4.216476416e-1; 
      exponents[STO5G][nShell][dAzimuthal][3] = 5.836797641e-2;   coefficients[STO5G][nShell][dAzimuthal][3] = 4.547673415e-1; 
      exponents[STO5G][nShell][dAzimuthal][4] = 3.206682246e-2;   coefficients[STO5G][nShell][dAzimuthal][4] = 1.037803318e-1; 
   }

   //STO-6G, 
   {
      // 1s
      exponents[STO6G][kShell][sAzimuthal][0] = 2.310303149e01;   coefficients[STO6G][kShell][sAzimuthal][0] = 9.163596280e-3; 
      exponents[STO6G][kShell][sAzimuthal][1] = 4.235915534e00;   coefficients[STO6G][kShell][sAzimuthal][1] = 4.936149294e-2; 
      exponents[STO6G][kShell][sAzimuthal][2] = 1.185056519e00;   coefficients[STO6G][kShell][sAzimuthal][2] = 1.685383049e-1; 
      exponents[STO6G][kShell][sAzimuthal][3] = 4.070988982e-1;   coefficients[STO6G][kShell][sAzimuthal][3] = 3.705627997e-1; 
      exponents[STO6G][kShell][sAzimuthal][4] = 1.580884151e-1;   coefficients[STO6G][kShell][sAzimuthal][4] = 4.164915298e-1; 
      exponents[STO6G][kShell][sAzimuthal][5] = 6.510953954e-2;   coefficients[STO6G][kShell][sAzimuthal][5] = 1.303340841e-1; 
      // 2s
      exponents[STO6G][lShell][sAzimuthal][0] = 2.768496241e01;   coefficients[STO6G][lShell][sAzimuthal][0] =-4.151277819e-3; 
      exponents[STO6G][lShell][sAzimuthal][1] = 5.077140627e00;   coefficients[STO6G][lShell][sAzimuthal][1] =-2.067024148e-2; 
      exponents[STO6G][lShell][sAzimuthal][2] = 1.426786050e00;   coefficients[STO6G][lShell][sAzimuthal][2] =-5.150303337e-2; 
      exponents[STO6G][lShell][sAzimuthal][3] = 2.040335729e-1;   coefficients[STO6G][lShell][sAzimuthal][3] = 3.346271174e-1; 
      exponents[STO6G][lShell][sAzimuthal][4] = 9.260298399e-2;   coefficients[STO6G][lShell][sAzimuthal][4] = 5.621061301e-1; 
      exponents[STO6G][lShell][sAzimuthal][5] = 4.416183978e-2;   coefficients[STO6G][lShell][sAzimuthal][5] = 1.712994697e-1; 
      // 2p
      exponents[STO6G][lShell][pAzimuthal][0] = 5.868285913e00;   coefficients[STO6G][lShell][pAzimuthal][0] = 7.924233646e-3; 
      exponents[STO6G][lShell][pAzimuthal][1] = 1.530329631e00;   coefficients[STO6G][lShell][pAzimuthal][1] = 5.144104825e-2; 
      exponents[STO6G][lShell][pAzimuthal][2] = 5.475665231e-1;   coefficients[STO6G][lShell][pAzimuthal][2] = 1.898400060e-1; 
      exponents[STO6G][lShell][pAzimuthal][3] = 2.288932733e-1;   coefficients[STO6G][lShell][pAzimuthal][3] = 4.049863191e-1; 
      exponents[STO6G][lShell][pAzimuthal][4] = 1.046655969e-1;   coefficients[STO6G][lShell][pAzimuthal][4] = 4.012362861e-1; 
      exponents[STO6G][lShell][pAzimuthal][5] = 4.948220127e-2;   coefficients[STO6G][lShell][pAzimuthal][5] = 1.051855189e-1; 
      // 3s
      exponents[STO6G][mShell][sAzimuthal][0] = 3.273031938e00;   coefficients[STO6G][mShell][sAzimuthal][0] =-6.775596947e-3; 
      exponents[STO6G][mShell][sAzimuthal][1] = 9.200611311e-1;   coefficients[STO6G][mShell][sAzimuthal][1] =-5.639325779e-2; 
      exponents[STO6G][mShell][sAzimuthal][2] = 3.593349765e-1;   coefficients[STO6G][mShell][sAzimuthal][2] =-1.587656086e-1; 
      exponents[STO6G][mShell][sAzimuthal][3] = 8.636686991e-2;   coefficients[STO6G][mShell][sAzimuthal][3] = 5.534527651e-1; 
      exponents[STO6G][mShell][sAzimuthal][4] = 4.797373812e-2;   coefficients[STO6G][mShell][sAzimuthal][4] = 5.015351020e-1; 
      exponents[STO6G][mShell][sAzimuthal][5] = 2.724741144e-2;   coefficients[STO6G][mShell][sAzimuthal][5] = 7.223633674e-2; 
      // 3p
      exponents[STO6G][mShell][pAzimuthal][0] = 5.077973607e00;   coefficients[STO6G][mShell][pAzimuthal][0] =-3.329929840e-3; 
      exponents[STO6G][mShell][pAzimuthal][1] = 1.340786940e00;   coefficients[STO6G][mShell][pAzimuthal][1] =-1.419488340e-2; 
      exponents[STO6G][mShell][pAzimuthal][2] = 2.248434849e-1;   coefficients[STO6G][mShell][pAzimuthal][2] = 1.639395770e-1; 
      exponents[STO6G][mShell][pAzimuthal][3] = 1.131741848e-1;   coefficients[STO6G][mShell][pAzimuthal][3] = 4.485358256e-1; 
      exponents[STO6G][mShell][pAzimuthal][4] = 6.076408893e-2;   coefficients[STO6G][mShell][pAzimuthal][4] = 3.908813050e-1; 
      exponents[STO6G][mShell][pAzimuthal][5] = 3.315424265e-2;   coefficients[STO6G][mShell][pAzimuthal][5] = 7.411456232e-2; 
      // 3d
      exponents[STO6G][mShell][dAzimuthal][0] = 2.488296923e00;   coefficients[STO6G][mShell][dAzimuthal][0] = 7.283828112e-3; 
      exponents[STO6G][mShell][dAzimuthal][1] = 7.981487853e-1;   coefficients[STO6G][mShell][dAzimuthal][1] = 5.386799363e-2; 
      exponents[STO6G][mShell][dAzimuthal][2] = 3.311327490e-1;   coefficients[STO6G][mShell][dAzimuthal][2] = 2.072139149e-1; 
      exponents[STO6G][mShell][dAzimuthal][3] = 1.559114463e-1;   coefficients[STO6G][mShell][dAzimuthal][3] = 4.266269092e-1; 
      exponents[STO6G][mShell][dAzimuthal][4] = 7.877734732e-2;   coefficients[STO6G][mShell][dAzimuthal][4] = 3.843100204e-1; 
      exponents[STO6G][mShell][dAzimuthal][5] = 4.058484363e-2;   coefficients[STO6G][mShell][dAzimuthal][5] = 8.902827546e-2; 
      // 4s
      exponents[STO6G][nShell][sAzimuthal][0] = 3.232838646e00;   coefficients[STO6G][nShell][sAzimuthal][0] = 1.374817488e-3; 
      exponents[STO6G][nShell][sAzimuthal][1] = 3.605788802e-1;   coefficients[STO6G][nShell][sAzimuthal][1] =-8.666390043e-2; 
      exponents[STO6G][nShell][sAzimuthal][2] = 1.717905487e-1;   coefficients[STO6G][nShell][sAzimuthal][2] =-3.130627309e-1; 
      exponents[STO6G][nShell][sAzimuthal][3] = 5.277666487e-2;   coefficients[STO6G][nShell][sAzimuthal][3] = 7.812787397e-1; 
      exponents[STO6G][nShell][sAzimuthal][4] = 3.163400284e-2;   coefficients[STO6G][nShell][sAzimuthal][4] = 4.389247988e-1; 
      exponents[STO6G][nShell][sAzimuthal][5] = 1.874093091e-2;   coefficients[STO6G][nShell][sAzimuthal][5] = 2.487178756e-2; 
      // 4p                                                               
      exponents[STO6G][nShell][pAzimuthal][0] = 2.389722618e00;   coefficients[STO6G][nShell][pAzimuthal][0] =-1.665913575e-3; 
      exponents[STO6G][nShell][pAzimuthal][1] = 7.960947826e-1;   coefficients[STO6G][nShell][pAzimuthal][1] =-1.657464971e-2; 
      exponents[STO6G][nShell][pAzimuthal][2] = 3.415541380e-1;   coefficients[STO6G][nShell][pAzimuthal][2] =-5.958513378e-2; 
      exponents[STO6G][nShell][pAzimuthal][3] = 8.847434525e-2;   coefficients[STO6G][nShell][pAzimuthal][3] = 4.053115554e-1; 
      exponents[STO6G][nShell][pAzimuthal][4] = 4.958248334e-2;   coefficients[STO6G][nShell][pAzimuthal][4] = 5.433958189e-1; 
      exponents[STO6G][nShell][pAzimuthal][5] = 2.816929784e-2;   coefficients[STO6G][nShell][pAzimuthal][5] = 1.204970491e-1; 
      // 4d                                                               
      exponents[STO6G][nShell][dAzimuthal][0] = 4.634239420e00;   coefficients[STO6G][nShell][dAzimuthal][0] =-4.749842876e-4; 
      exponents[STO6G][nShell][dAzimuthal][1] = 1.341648295e00;   coefficients[STO6G][nShell][dAzimuthal][1] =-3.566777891e-3; 
      exponents[STO6G][nShell][dAzimuthal][2] = 2.209593028e-1;   coefficients[STO6G][nShell][dAzimuthal][2] = 1.108670481e-1; 
      exponents[STO6G][nShell][dAzimuthal][3] = 1.101467943e-1;   coefficients[STO6G][nShell][dAzimuthal][3] = 4.159646930e-1; 
      exponents[STO6G][nShell][dAzimuthal][4] = 5.904190370e-2;   coefficients[STO6G][nShell][dAzimuthal][4] = 4.621672517e-1; 
      exponents[STO6G][nShell][dAzimuthal][5] = 3.232628887e-2;   coefficients[STO6G][nShell][dAzimuthal][5] = 1.081250196e-1; 
   }

}






double Cpp_Cndo2GetGaussianOverlapAOsSASB(double gaussianExponentA, 
                                        double gaussianExponentB,
                                        double rAB) {
   double value;
   double temp1 = 0.0;
   double temp2 = 0.0;
   double gauPlusAB = gaussianExponentA+gaussianExponentB;
   double gauMultAB = gaussianExponentA*gaussianExponentB;
   temp1 = 2.0*sqrt(gauMultAB)/gauPlusAB;
   temp2 = -1.0* gauMultAB/gauPlusAB;
   value = temp1*sqrt(temp1)*exp(temp2*rAB*rAB);
   return value;
}

double Cpp_Cndo2GetGaussianOverlapAOs(
                                    int valenceOrbitalA, 
                                    double gaussianExponentA, 
                                    int valenceOrbitalB, 
                                    double gaussianExponentB,
                                    double dx, double dy, double dz, 
                                    double rAB,
                                    double overlapSASB) {

   double value = 0.0;
   double gauPlusAB = gaussianExponentA+gaussianExponentB;
   double gauMultAB = gaussianExponentA*gaussianExponentB;
   if(valenceOrbitalA == s && valenceOrbitalB == s){
      value = 1.0;
   }

   else if(valenceOrbitalA == s && valenceOrbitalB == px){
      value = 2.0*gaussianExponentA*sqrt(gaussianExponentB)*dx;
      value /= gauPlusAB;
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == py){
      value = 2.0*gaussianExponentA*sqrt(gaussianExponentB)*dy;
      value /= gauPlusAB;
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == pz){
      value = 2.0*gaussianExponentA*sqrt(gaussianExponentB)*dz;
      value /= gauPlusAB;
   }

   else if(valenceOrbitalA == px && valenceOrbitalB == s){
      value = -2.0*sqrt(gaussianExponentA)*gaussianExponentB*dx;
      value /= gauPlusAB;
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == s){
      value = -2.0*sqrt(gaussianExponentA)*gaussianExponentB*dy;
      value /= gauPlusAB;
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == s){
      value = -2.0*sqrt(gaussianExponentA)*gaussianExponentB*dz;
      value /= gauPlusAB;
   }

   else if(valenceOrbitalA == px && valenceOrbitalB == px){
      double temp = 0.0;
      temp = -1.0*(dx*dx)*gauMultAB;
      temp /= gauPlusAB;
      temp += 0.5;
      value = 4.0*sqrt(gauMultAB);
      value /= gauPlusAB;
      value *= temp;
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == py){
      value = -4.0*gauMultAB*sqrt(gauMultAB);
      value /= gauPlusAB*gauPlusAB;
      value *= dx*dy;
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == pz){
      value = -4.0*gauMultAB*sqrt(gauMultAB);
      value /= gauPlusAB*gauPlusAB;
      value *= dx*dz;
   }

   else if(valenceOrbitalA == py && valenceOrbitalB == px){
      value = -4.0*gauMultAB*sqrt(gauMultAB);
      value /= gauPlusAB*gauPlusAB;
      value *= dy*dx;
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == py){
      double temp = 0.0;
      temp = -1.0*(dy*dy)*gauMultAB;
      temp /= gauPlusAB;
      temp += 0.5;
      value = 4.0*sqrt(gauMultAB);
      value /= gauPlusAB;
      value *= temp;
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == pz){
      value = -4.0*gauMultAB*sqrt(gauMultAB);
      value /= gauPlusAB*gauPlusAB;
      value *= dy*dz;
   }

   else if(valenceOrbitalA == pz && valenceOrbitalB == px){
      value = -4.0*gauMultAB*sqrt(gauMultAB);
      value /= gauPlusAB*gauPlusAB;
      value *= dz*dx;
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == py){
      value = -4.0*gauMultAB*sqrt(gauMultAB);
      value /= gauPlusAB*gauPlusAB;
      value *= dz*dy;
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == pz){
      double temp = 0.0;
      temp = -1.0*(dz*dz)*gauMultAB;
      temp /= gauPlusAB;
      temp += 0.5;
      value = 4.0*sqrt(gauMultAB);
      value /= gauPlusAB;
      value *= temp;
   }

   else if(valenceOrbitalA == dxy && valenceOrbitalB == s){
      value = 4.0*gaussianExponentA;
      value /= gauPlusAB*gauPlusAB;
      value *= gaussianExponentB*dx;
      value *= gaussianExponentB*dy;
   }
   else if(valenceOrbitalA == dyz && valenceOrbitalB == s){
      value = 4.0*gaussianExponentA;
      value /= gauPlusAB*gauPlusAB;
      value *= gaussianExponentB*dy;
      value *= gaussianExponentB*dz;
   }
   else if(valenceOrbitalA == dzx && valenceOrbitalB == s){
      value = 4.0*gaussianExponentA;
      value /= gauPlusAB*gauPlusAB;
      value *= gaussianExponentB*dz;
      value *= gaussianExponentB*dx;
   }

   else if(valenceOrbitalA == s && valenceOrbitalB == dxy){
      value = 4.0*gaussianExponentB;
      value /= gauPlusAB*gauPlusAB;
      value *= gaussianExponentA*dx;
      value *= gaussianExponentA*dy;
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == dyz){
      value = 4.0*gaussianExponentB;
      value /= gauPlusAB*gauPlusAB;
      value *= gaussianExponentA*dy;
      value *= gaussianExponentA*dz;
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == dzx){
      value = 4.0*gaussianExponentB;
      value /= gauPlusAB*gauPlusAB;
      value *= gaussianExponentA*dz;
      value *= gaussianExponentA*dx;
   }

   else if(valenceOrbitalA == dxy && valenceOrbitalB == px){
      double temp1 = -0.5*gaussianExponentB*dy;
      double temp2 = gaussianExponentB*dx*gaussianExponentA*dx*gaussianExponentB*dy;
      temp2 /= gauPlusAB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentA*sqrt(gaussianExponentB);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == dxy && valenceOrbitalB == py){
      double temp1 = -0.5*gaussianExponentB*dx;
      double temp2 = gaussianExponentB*dy*gaussianExponentA*dy*gaussianExponentB*dx;
      temp2 /= gauPlusAB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentA*sqrt(gaussianExponentB);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == dyz && valenceOrbitalB == py){
      double temp1 = -0.5*gaussianExponentB*dz;
      double temp2 = gaussianExponentB*dy*gaussianExponentA*dy*gaussianExponentB*dz;
      temp2 /= gauPlusAB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentA*sqrt(gaussianExponentB);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == dyz && valenceOrbitalB == pz){
      double temp1 = -0.5*gaussianExponentB*dy;
      double temp2 = gaussianExponentB*dz*gaussianExponentA*dz*gaussianExponentB*dy;
      temp2 /= gauPlusAB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentA*sqrt(gaussianExponentB);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == dzx && valenceOrbitalB == pz){
      double temp1 = -0.5*gaussianExponentB*dx;
      double temp2 = gaussianExponentB*dz*gaussianExponentA*dz*gaussianExponentB*dx;
      temp2 /= gauPlusAB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentA*sqrt(gaussianExponentB);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == dzx && valenceOrbitalB == px){
      double temp1 = -0.5*gaussianExponentB*dz;
      double temp2 = gaussianExponentB*dx*gaussianExponentA*dx*gaussianExponentB*dz;
      temp2 /= gauPlusAB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentA*sqrt(gaussianExponentB);
      value /= (gauPlusAB*gauPlusAB);
   }

   else if(valenceOrbitalA == px && valenceOrbitalB == dxy){
      double temp1 = 0.5*gaussianExponentA*dy;
      double temp2 = -1.0*gaussianExponentA*dx*gaussianExponentB*dx*gaussianExponentA*dy;
      temp2 /= gauPlusAB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentB*sqrt(gaussianExponentA);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == dxy){
      double temp1 = 0.5*gaussianExponentA*dx;
      double temp2 = -1.0*gaussianExponentA*dy*gaussianExponentB*dy*gaussianExponentA*dx;
      temp2 /= gauPlusAB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentB*sqrt(gaussianExponentA);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == dyz){
      double temp1 = 0.5*gaussianExponentA*dz;
      double temp2 = -1.0*gaussianExponentA*dy*gaussianExponentB*dy*gaussianExponentA*dz;
      temp2 /= gauPlusAB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentB*sqrt(gaussianExponentA);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == dyz){
      double temp1 = 0.5*gaussianExponentA*dy;
      double temp2 = -1.0*gaussianExponentA*dz*gaussianExponentB*dz*gaussianExponentA*dy;
      temp2 /= gauPlusAB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentB*sqrt(gaussianExponentA);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == dzx){
      double temp1 = 0.5*gaussianExponentA*dx;
      double temp2 = -1.0*gaussianExponentA*dz*gaussianExponentB*dz*gaussianExponentA*dx;
      temp2 /= gauPlusAB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentB*sqrt(gaussianExponentA);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == dzx){
      double temp1 = 0.5*gaussianExponentA*dz;
      double temp2 = -1.0*gaussianExponentA*dx*gaussianExponentB*dx*gaussianExponentA*dz;
      temp2 /= gauPlusAB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentB*sqrt(gaussianExponentA);
      value /= (gauPlusAB*gauPlusAB);
   }

   else if(valenceOrbitalA == dxy && valenceOrbitalB == pz){
      value = 8.0*gaussianExponentA*sqrt(gaussianExponentB);
      value /= (gauPlusAB*gauPlusAB*gauPlusAB);
      value *= gaussianExponentB*dx*gaussianExponentB*dy*gaussianExponentA*dz;
   }
   else if(valenceOrbitalA == dyz && valenceOrbitalB == px){
      value = 8.0*gaussianExponentA*sqrt(gaussianExponentB);
      value /= (gauPlusAB*gauPlusAB*gauPlusAB);
      value *= gaussianExponentB*dy*gaussianExponentB*dz*gaussianExponentA*dx;
   }
   else if(valenceOrbitalA == dzx && valenceOrbitalB == py){
      value = 8.0*gaussianExponentA*sqrt(gaussianExponentB);
      value /= (gauPlusAB*gauPlusAB*gauPlusAB);
      value *= gaussianExponentB*dz*gaussianExponentB*dx*gaussianExponentA*dy;
   }

   else if(valenceOrbitalA == pz && valenceOrbitalB == dxy){
      value = -8.0*gaussianExponentB*sqrt(gaussianExponentA);
      value /= (gauPlusAB*gauPlusAB*gauPlusAB);
      value *= gaussianExponentA*dx*gaussianExponentA*dy*gaussianExponentB*dz;
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == dyz){
      value = -8.0*gaussianExponentB*sqrt(gaussianExponentA);
      value /= (gauPlusAB*gauPlusAB*gauPlusAB);
      value *= gaussianExponentA*dy*gaussianExponentA*dz*gaussianExponentB*dx;
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == dzx){
      value = -8.0*gaussianExponentB*sqrt(gaussianExponentA);
      value /= (gauPlusAB*gauPlusAB*gauPlusAB);
      value *= gaussianExponentA*dz*gaussianExponentA*dx*gaussianExponentB*dy;
   }

   else if(valenceOrbitalA == dxxyy && valenceOrbitalB == s){
      value = 2.0*gaussianExponentA;
      value /= gauPlusAB*gauPlusAB;
      value *= (gaussianExponentB*gaussianExponentB*dx*dx) - (gaussianExponentB*gaussianExponentB*dy*dy);
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == dxxyy){
      value = 2.0*gaussianExponentB;
      value /= gauPlusAB*gauPlusAB;
      value *= (gaussianExponentA*gaussianExponentA*dx*dx) - (gaussianExponentA*gaussianExponentA*dy*dy);
   }
   else if(valenceOrbitalA == dxxyy && valenceOrbitalB == px){
      value = gaussianExponentB*dx;
      value -= (gaussianExponentB*gaussianExponentB*dx*dx)*gaussianExponentA*dx/gauPlusAB;
      value += (gaussianExponentB*gaussianExponentB*dy*dy)*gaussianExponentA*dx/gauPlusAB;
      value *= -4.0*gaussianExponentA*sqrt(gaussianExponentB);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == dxxyy){
      value = gaussianExponentA*dx;
      value -= (gaussianExponentA*gaussianExponentA*dx*dx)*gaussianExponentB*dx/gauPlusAB;
      value += (gaussianExponentA*gaussianExponentA*dy*dy)*gaussianExponentB*dx/gauPlusAB;
      value *= 4.0*gaussianExponentB*sqrt(gaussianExponentA);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == dxxyy && valenceOrbitalB == py){
      value = gaussianExponentB*dy;
      value += (gaussianExponentB*gaussianExponentB*dx*dx)*gaussianExponentA*dy/gauPlusAB;
      value -= (gaussianExponentB*gaussianExponentB*dy*dy)*gaussianExponentA*dy/gauPlusAB;
      value *= 4.0*gaussianExponentA*sqrt(gaussianExponentB);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == dxxyy){
      value = gaussianExponentA*dy;
      value += (gaussianExponentA*gaussianExponentA*dx*dx)*gaussianExponentB*dy/gauPlusAB;
      value -= (gaussianExponentA*gaussianExponentA*dy*dy)*gaussianExponentB*dy/gauPlusAB;
      value *= -4.0*gaussianExponentB*sqrt(gaussianExponentA);
      value /= (gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == dxxyy && valenceOrbitalB == pz){
      value = (gaussianExponentB*gaussianExponentB*dx*dx) - (gaussianExponentB*gaussianExponentB*dy*dy);
      value *= gaussianExponentA*dz;
      value *= 4.0*gaussianExponentA*sqrt(gaussianExponentB);
      value /= (gauPlusAB*gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == dxxyy){
      value = (gaussianExponentA*gaussianExponentA*dx*dx) - (gaussianExponentA*gaussianExponentA*dy*dy);
      value *= gaussianExponentB*dz;
      value *= -4.0*gaussianExponentB*sqrt(gaussianExponentA);
      value /= (gauPlusAB*gauPlusAB*gauPlusAB);
   }

   else if(valenceOrbitalA == dzz && valenceOrbitalB == s){
      double temp = 0.0;
      temp = 2.0*(gaussianExponentB*gaussianExponentB*dz*dz) 
            -    (gaussianExponentB*gaussianExponentB*dx*dx) 
            -    (gaussianExponentB*gaussianExponentB*dy*dy);
      value = 2.0*gaussianExponentA/sqrt(3.0);
      value /= gauPlusAB*gauPlusAB;
      value *= temp;
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == dzz){
      double temp = 0.0;
      temp = 2.0*(gaussianExponentA*gaussianExponentA*dz*dz) 
            -    (gaussianExponentA*gaussianExponentA*dx*dx) 
            -    (gaussianExponentA*gaussianExponentA*dy*dy);
      value = 2.0*gaussianExponentB/sqrt(3.0);
      value /= gauPlusAB*gauPlusAB;
      value *= temp;
   }
   else if(valenceOrbitalA == dzz && valenceOrbitalB == px){
      double temp = 0.0;
      temp = gaussianExponentB*dx;
      temp += 2.0*(gaussianExponentB*gaussianExponentB*dz*dz)*gaussianExponentA*dx/gauPlusAB;
      temp -=     (gaussianExponentB*gaussianExponentB*dx*dx)*gaussianExponentA*dx/gauPlusAB;
      temp -=     (gaussianExponentB*gaussianExponentB*dy*dy)*gaussianExponentA*dx/gauPlusAB;
      value = temp;
      value *= 4.0*gaussianExponentA*sqrt(gaussianExponentB)/sqrt(3.0);
      value /= gauPlusAB*gauPlusAB;
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == dzz){
      double temp = 0.0;
      temp = gaussianExponentA*dx;
      temp += 2.0*(gaussianExponentA*gaussianExponentA*dz*dz)*gaussianExponentB*dx/gauPlusAB;
      temp -=     (gaussianExponentA*gaussianExponentA*dx*dx)*gaussianExponentB*dx/gauPlusAB;
      temp -=     (gaussianExponentA*gaussianExponentA*dy*dy)*gaussianExponentB*dx/gauPlusAB;
      value = temp;
      value *= -4.0*gaussianExponentB*sqrt(gaussianExponentA)/sqrt(3.0);
      value /= gauPlusAB*gauPlusAB;
   }
   else if(valenceOrbitalA == dzz && valenceOrbitalB == py){
      double temp = 0.0;
      temp = gaussianExponentB*dy;
      temp += 2.0*(gaussianExponentB*gaussianExponentB*dz*dz)*gaussianExponentA*dy/gauPlusAB;
      temp -=     (gaussianExponentB*gaussianExponentB*dx*dx)*gaussianExponentA*dy/gauPlusAB;
      temp -=     (gaussianExponentB*gaussianExponentB*dy*dy)*gaussianExponentA*dy/gauPlusAB;
      value = temp;
      value *= 4.0*gaussianExponentA*sqrt(gaussianExponentB)/sqrt(3.0);
      value /= gauPlusAB*gauPlusAB;
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == dzz){
      double temp = 0.0;
      temp = gaussianExponentA*dy;
      temp += 2.0*(gaussianExponentA*gaussianExponentA*dz*dz)*gaussianExponentB*dy/gauPlusAB;
      temp -=     (gaussianExponentA*gaussianExponentA*dx*dx)*gaussianExponentB*dy/gauPlusAB;
      temp -=     (gaussianExponentA*gaussianExponentA*dy*dy)*gaussianExponentB*dy/gauPlusAB;
      value = temp;
      value *= -4.0*gaussianExponentB*sqrt(gaussianExponentA)/sqrt(3.0);
      value /= gauPlusAB*gauPlusAB;
   }
   else if(valenceOrbitalA == dzz && valenceOrbitalB == pz){
      double temp = 0.0;
      temp = -2.0*gaussianExponentB*dz;
      temp += 2.0*(gaussianExponentB*gaussianExponentB*dz*dz)*gaussianExponentA*dz/gauPlusAB;
      temp -=     (gaussianExponentB*gaussianExponentB*dx*dx)*gaussianExponentA*dz/gauPlusAB;
      temp -=     (gaussianExponentB*gaussianExponentB*dy*dy)*gaussianExponentA*dz/gauPlusAB;
      value = temp;
      value *= 4.0*gaussianExponentA*sqrt(gaussianExponentB)/sqrt(3.0);
      value /= gauPlusAB*gauPlusAB;
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == dzz){
      double temp = 0.0;
      temp = -2.0*gaussianExponentA*dz;
      temp += 2.0*(gaussianExponentA*gaussianExponentA*dz*dz)*gaussianExponentB*dz/gauPlusAB;
      temp -=     (gaussianExponentA*gaussianExponentA*dx*dx)*gaussianExponentB*dz/gauPlusAB;
      temp -=     (gaussianExponentA*gaussianExponentA*dy*dy)*gaussianExponentB*dz/gauPlusAB;
      value = temp;
      value *= -4.0*gaussianExponentB*sqrt(gaussianExponentA)/sqrt(3.0);
      value /= gauPlusAB*gauPlusAB;
   }

   else if(valenceOrbitalA == dxy && valenceOrbitalB == dxy){
      double temp = 0.25;
      temp -= 0.5*gaussianExponentB*dx*gaussianExponentA*dx/gauPlusAB;
      temp -= 0.5*gaussianExponentB*dy*gaussianExponentA*dy/gauPlusAB;
      temp += gaussianExponentB*dx*gaussianExponentA*dx
             *gaussianExponentB*dy*gaussianExponentA*dy
             /(gauPlusAB*gauPlusAB);
      value = 16.0*temp*gauMultAB
             /(gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == dyz && valenceOrbitalB == dyz){
      double temp = 0.25;
      temp -= 0.5*gaussianExponentB*dy*gaussianExponentA*dy/gauPlusAB;
      temp -= 0.5*gaussianExponentB*dz*gaussianExponentA*dz/gauPlusAB;
      temp += gaussianExponentB*dy*gaussianExponentA*dy
             *gaussianExponentB*dz*gaussianExponentA*dz
             /(gauPlusAB*gauPlusAB);
      value = 16.0*temp*gauMultAB
             /(gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == dzx && valenceOrbitalB == dzx){
      double temp = 0.25;
      temp -= 0.5*gaussianExponentB*dz*gaussianExponentA*dz/gauPlusAB;
      temp -= 0.5*gaussianExponentB*dx*gaussianExponentA*dx/gauPlusAB;
      temp += gaussianExponentB*dz*gaussianExponentA*dz
             *gaussianExponentB*dx*gaussianExponentA*dx
             /(gauPlusAB*gauPlusAB);
      value = 16.0*temp*gauMultAB
             /(gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == dxxyy && valenceOrbitalB == dxxyy){
      double temp1 = 1.0;
      temp1 -= 2.0*gaussianExponentB*dx*gaussianExponentA*dx/gauPlusAB;
      temp1 -= 2.0*gaussianExponentB*dy*gaussianExponentA*dy/gauPlusAB;
      double temp2 = gauMultAB*((dx*dx)-(dy*dy))
             /gauPlusAB;
      temp1 += temp2*temp2;
      value = 4.0*temp1*gauMultAB
             /(gauPlusAB*gauPlusAB);
   }
   else if(valenceOrbitalA == dzz && valenceOrbitalB == dzz){
      double temp = 3.0;
      temp -= gauMultAB
             /gauPlusAB
             *(8.0*(dz*dz)+2.0*(dx*dx)+2.0*(dy*dy));
      temp += (gaussianExponentA*gauMultAB*gaussianExponentB)
             *(4.0*(dz*dz*dz*dz)
                  +(dx*dx*dx*dx)
                  +(dy*dy*dy*dy)
              -4.0*(dx*dx*dz*dz)
              -4.0*(dy*dy*dz*dz)
              +2.0*(dx*dx*dy*dy))
             /(gauPlusAB*gauPlusAB);
      value = 4.0*temp*gauMultAB
             /(3.0*gauPlusAB*gauPlusAB);
   }

   else if((valenceOrbitalA == dxy && valenceOrbitalB == dyz) ||
           (valenceOrbitalA == dyz && valenceOrbitalB == dxy)){
      double temp = 0.5;
      temp -= gauMultAB*(dy*dy)/gauPlusAB;
      value = -16.0*(gaussianExponentA*gauMultAB*gaussianExponentB)
             *dx*dz*temp/(gauPlusAB*gauPlusAB*gauPlusAB);
   }
   else if((valenceOrbitalA == dyz && valenceOrbitalB == dzx) ||
           (valenceOrbitalA == dzx && valenceOrbitalB == dyz)){
      double temp = 0.5;
      temp -= gauMultAB*(dz*dz)/gauPlusAB;
      value = -16.0*(gaussianExponentA*gauMultAB*gaussianExponentB)
             *dy*dx*temp/(gauPlusAB*gauPlusAB*gauPlusAB);
   }
   else if((valenceOrbitalA == dzx && valenceOrbitalB == dxy) ||
           (valenceOrbitalA == dxy && valenceOrbitalB == dzx)){
      double temp = 0.5;
      temp -= gauMultAB*(dx*dx)/gauPlusAB;
      value = -16.0*(gaussianExponentA*gauMultAB*gaussianExponentB)
             *dz*dy*temp/(gauPlusAB*gauPlusAB*gauPlusAB);
   }

   else if((valenceOrbitalA == dxxyy && valenceOrbitalB == dxy) ||
           (valenceOrbitalA == dxy && valenceOrbitalB == dxxyy)){
      double temp = 2.0*gauMultAB;
      value = (temp*temp*temp)*(dy*(dx*dx*dx)-dx*(dy*dy*dy))
             /(gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB);
   }
   else if((valenceOrbitalA == dxxyy && valenceOrbitalB == dyz) ||
           (valenceOrbitalA == dyz && valenceOrbitalB == dxxyy)){
      double temp = 2.0*gauMultAB;
      value = (temp*temp*temp)*(dy*dz*gauPlusAB
                            /gauMultAB
                             +((dx*dx)*dy*dz - (dy*dy*dy)*dz))
             /(gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB);
   }
   else if((valenceOrbitalA == dxxyy && valenceOrbitalB == dzx) ||
           (valenceOrbitalA == dzx && valenceOrbitalB == dxxyy)){
      double temp = 2.0*gauMultAB;
      value = -1.0*(temp*temp*temp)*(dx*dz*gauPlusAB
                            /gauMultAB
                             +((dy*dy)*dx*dz - (dx*dx*dx)*dz))
             /(gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB);
   }

   else if((valenceOrbitalA == dzz && valenceOrbitalB == dxy) ||
           (valenceOrbitalA == dxy && valenceOrbitalB == dzz)){
      double temp = 2.0*dx*dy*(dz*dz) - (dx*dx*dx)*dy - dx*(dy*dy*dy);
      temp *= gauMultAB/gauPlusAB;
      temp += 2.0*dx*dy;
      value = 8.0*(gaussianExponentA*gauMultAB*gaussianExponentB)*temp;
      value /= sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB);
   }
   else if((valenceOrbitalA == dzz && valenceOrbitalB == dyz) ||
           (valenceOrbitalA == dyz && valenceOrbitalB == dzz)){
      double temp1 = -1.0*dy*dz;
      double temp2 = 2.0*dy*(dz*dz*dz) - (dy*dy*dy)*dz - (dx*dx)*dy*dz;
      temp2 *= gauMultAB/gauPlusAB;
      temp1 += temp2;
      value = 8.0*(gaussianExponentA*gauMultAB*gaussianExponentB)*temp1;
      value /= sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB);
   }
   else if((valenceOrbitalA == dzz && valenceOrbitalB == dzx) ||
           (valenceOrbitalA == dzx && valenceOrbitalB == dzz)){
      double temp1 = -1.0*dx*dz;
      double temp2 = 2.0*dx*(dz*dz*dz) - (dx*dx*dx)*dz - (dy*dy)*dx*dz;
      temp2 *= gauMultAB/gauPlusAB;
      temp1 += temp2;
      value = 8.0*(gaussianExponentA*gauMultAB*gaussianExponentB)*temp1;
      value /= sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB);
   }
   else if((valenceOrbitalA == dxxyy && valenceOrbitalB == dzz) ||
           (valenceOrbitalA == dzz && valenceOrbitalB == dxxyy)){
      double temp = 2.0*(dz*dz)-(dx*dx)-(dy*dy);
      temp *= gauMultAB/gauPlusAB;
      temp += 2.0;
      value = 4.0*(gaussianExponentA*gauMultAB*gaussianExponentB);
      value /= sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= ((dx*dx)-(dy*dy))*temp;
   }

   else{
      mexErrMsgTxt("Cpp_Cndo2GetGaussianOverlapAOs: Orbital type wrong.");
      value = 0.0;
   }
   value *= overlapSASB;
   return value;
}

double Cpp_Cndo2GetGaussianCartesianMatrix(int valenceOrbitalA, 
                                         double gaussianExponentA, 
                                         double const* xyzA,
                                         int valenceOrbitalB, 
                                         double gaussianExponentB,
                                         double const* xyzB,
                                         double rAB,
                                         double overlapSASB,
                                         int axis) {

   double value = 0.0;
   double gauPlusAB = gaussianExponentA+gaussianExponentB;
   double gauMultAB = gaussianExponentA*gaussianExponentB;
   double dxyz[3] = {xyzA[0] - xyzB[0], 
                     xyzA[1] - xyzB[1], 
                     xyzA[2] - xyzB[2]};
   if(valenceOrbitalA == s && valenceOrbitalB == s){
      value = gaussianExponentA*xyzA[axis] + gaussianExponentB*xyzB[axis];
      value /= gauPlusAB;
      value *= overlapSASB;
      return value;
   }
   else if( (valenceOrbitalA == s && axis == 0 && valenceOrbitalB == px) || 
            (valenceOrbitalA == s && axis == 0 && valenceOrbitalB == py) || 
            (valenceOrbitalA == s && axis == 0 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == s && axis == 0 && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == s && axis == 0 && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == s && axis == 0 && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == s && axis == 0 && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == s && axis == 0 && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == s && axis == 1 && valenceOrbitalB == px) || 
            (valenceOrbitalA == s && axis == 1 && valenceOrbitalB == py) || 
            (valenceOrbitalA == s && axis == 1 && valenceOrbitalB == pz) ||
            (valenceOrbitalA == s && axis == 1 && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == s && axis == 1 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == s && axis == 1 && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == s && axis == 1 && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == s && axis == 1 && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == s && axis == 2 && valenceOrbitalB == px) || 
            (valenceOrbitalA == s && axis == 2 && valenceOrbitalB == py) || 
            (valenceOrbitalA == s && axis == 2 && valenceOrbitalB == pz) ||
            (valenceOrbitalA == s && axis == 2 && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == s && axis == 2 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == s && axis == 2 && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == s && axis == 2 && valenceOrbitalB == dxxyy) ||
            (valenceOrbitalA == s && axis == 2 && valenceOrbitalB == dzz) ){
      int pOrbital;
      if(axis == 0){
         pOrbital = px;
      }
      else if(axis == 1){
         pOrbital = py;
      }
      else if(axis == 2){
         pOrbital = pz;
      }
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      double overlapAOs2 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       pOrbital, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = overlapAOs2/(2.0*sqrt(gaussianExponentA))+xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == px    && axis == 0 && valenceOrbitalB == s) || 
            (valenceOrbitalA == py    && axis == 0 && valenceOrbitalB == s) || 
            (valenceOrbitalA == pz    && axis == 0 && valenceOrbitalB == s) || 
            (valenceOrbitalA == dxy   && axis == 0 && valenceOrbitalB == s) || 
            (valenceOrbitalA == dyz   && axis == 0 && valenceOrbitalB == s) || 
            (valenceOrbitalA == dzx   && axis == 0 && valenceOrbitalB == s) ||
            (valenceOrbitalA == dxxyy && axis == 0 && valenceOrbitalB == s) ||
            (valenceOrbitalA == dzz   && axis == 0 && valenceOrbitalB == s) ||
            (valenceOrbitalA == px    && axis == 1 && valenceOrbitalB == s) || 
            (valenceOrbitalA == py    && axis == 1 && valenceOrbitalB == s) || 
            (valenceOrbitalA == pz    && axis == 1 && valenceOrbitalB == s) ||
            (valenceOrbitalA == dxy   && axis == 1 && valenceOrbitalB == s) ||
            (valenceOrbitalA == dyz   && axis == 1 && valenceOrbitalB == s) ||
            (valenceOrbitalA == dzx   && axis == 1 && valenceOrbitalB == s) ||
            (valenceOrbitalA == dxxyy && axis == 1 && valenceOrbitalB == s) ||
            (valenceOrbitalA == dzz   && axis == 1 && valenceOrbitalB == s) ||
            (valenceOrbitalA == px    && axis == 2 && valenceOrbitalB == s) || 
            (valenceOrbitalA == py    && axis == 2 && valenceOrbitalB == s) || 
            (valenceOrbitalA == pz    && axis == 2 && valenceOrbitalB == s) ||
            (valenceOrbitalA == dxy   && axis == 2 && valenceOrbitalB == s) ||
            (valenceOrbitalA == dyz   && axis == 2 && valenceOrbitalB == s) ||
            (valenceOrbitalA == dzx   && axis == 2 && valenceOrbitalB == s) ||
            (valenceOrbitalA == dxxyy && axis == 2 && valenceOrbitalB == s) ||
            (valenceOrbitalA == dzz   && axis == 2 && valenceOrbitalB == s) ){
      int pOrbital;
      if(axis == 0){
         pOrbital = px;
      }
      else if(axis == 1){
         pOrbital = py;
      }
      else if(axis == 2){
         pOrbital = pz;
      }
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      double overlapAOs2 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       pOrbital, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = overlapAOs2/(2.0*sqrt(gaussianExponentB))+xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == 0 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == py && axis == 1 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == pz && axis == 2 && valenceOrbitalB == valenceOrbitalA) ){
      double temp1 = gaussianExponentA*xyzA[axis] + gaussianExponentB*xyzB[axis];
      double temp2 = gaussianExponentA*xyzA[axis] - gaussianExponentA*xyzB[axis];
      double temp3 = gaussianExponentB*xyzA[axis] - gaussianExponentB*xyzB[axis];
      value = 0.5*(temp1+temp2-temp3);
      value -= temp1*temp2*temp3/gauPlusAB;
      value *= 4.0*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == 1 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == px && axis == 2 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == py && axis == 0 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == py && axis == 2 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == pz && axis == 0 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == pz && axis == 1 && valenceOrbitalB == valenceOrbitalA) ){
      int piDirection;
      if(valenceOrbitalA == px){
         piDirection = 0;
      }
      else if(valenceOrbitalA == py){
         piDirection = 1;
      }
      else if(valenceOrbitalA == pz){
         piDirection = 2;
      }
      double temp1 = gaussianExponentA*xyzA[piDirection] - gaussianExponentA*xyzB[piDirection];
      double temp2 = gaussianExponentB*xyzA[piDirection] - gaussianExponentB*xyzB[piDirection];
      value = 0.5 - temp1*temp2/gauPlusAB;
      value *= gaussianExponentA*xyzA[axis] + gaussianExponentB*xyzB[axis];
      value *= 4.0*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == 1 && valenceOrbitalB == py) || 
            (valenceOrbitalA == px && axis == 2 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == py && axis == 0 && valenceOrbitalB == px) || 
            (valenceOrbitalA == py && axis == 2 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == pz && axis == 0 && valenceOrbitalB == px) || 
            (valenceOrbitalA == pz && axis == 1 && valenceOrbitalB == py) ){
      int piDirectionA;
      if(valenceOrbitalA == px){
         piDirectionA = 0;
      }
      else if(valenceOrbitalA == py){
         piDirectionA = 1;
      }
      else if(valenceOrbitalA == pz){
         piDirectionA = 2;
      }
      double temp1 = gaussianExponentA*xyzA[axis] + gaussianExponentB*xyzB[axis];
      double temp2 = gaussianExponentA*xyzA[axis] - gaussianExponentA*xyzB[axis];
      value = 0.5 + temp1*temp2/gauPlusAB;
      value *= gaussianExponentB*xyzA[piDirectionA] - gaussianExponentB*xyzB[piDirectionA];
      value *= -4.0*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      return value;
   }
   else if( (valenceOrbitalA == py && axis == 1 && valenceOrbitalB == px) || 
            (valenceOrbitalA == py && axis == 1 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == px && axis == 0 && valenceOrbitalB == py) || 
            (valenceOrbitalA == px && axis == 0 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == pz && axis == 2 && valenceOrbitalB == px) || 
            (valenceOrbitalA == pz && axis == 2 && valenceOrbitalB == py) ){
      int piDirectionB;
      if(valenceOrbitalB == px){
         piDirectionB = 0;
      }
      else if(valenceOrbitalB == py){
         piDirectionB = 1;
      }
      else if(valenceOrbitalB == pz){
         piDirectionB = 2;
      }
      double temp1 = gaussianExponentA*xyzA[axis] + gaussianExponentB*xyzB[axis];
      double temp2 = gaussianExponentB*xyzA[axis] - gaussianExponentB*xyzB[axis];
      value = 0.5 - temp1*temp2/gauPlusAB;
      value *= gaussianExponentA*xyzA[piDirectionB] - gaussianExponentA*xyzB[piDirectionB];
      value *= 4.0*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == 1 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == py && axis == 2 && valenceOrbitalB == px) || 
            (valenceOrbitalA == pz && axis == 0 && valenceOrbitalB == py) ||
            (valenceOrbitalA == px && axis == 2 && valenceOrbitalB == py) ||
            (valenceOrbitalA == pz && axis == 1 && valenceOrbitalB == px) ||
            (valenceOrbitalA == py && axis == 0 && valenceOrbitalB == pz) ){
      int piDirectionA;
      int piDirectionB;
      if(valenceOrbitalA == px){
         piDirectionA = 0;
      }
      else if(valenceOrbitalA == py){
         piDirectionA = 1;
      }
      else if(valenceOrbitalA == pz){
         piDirectionA = 2;
      }
      if(valenceOrbitalB == px){
         piDirectionB = 0;
      }
      else if(valenceOrbitalB == py){
         piDirectionB = 1;
      }
      else if(valenceOrbitalB == pz){
         piDirectionB = 2;
      }
      double temp1 = gaussianExponentB*xyzA[piDirectionA] - gaussianExponentB*xyzB[piDirectionA];
      double temp2 = gaussianExponentA*xyzA[axis]         + gaussianExponentB*xyzB[axis];
      double temp3 = gaussianExponentA*xyzA[piDirectionB] - gaussianExponentA*xyzB[piDirectionB];
      value = -4.0*sqrt(gauMultAB)/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= temp1*temp2*temp3;
      value *= overlapSASB;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == 0 && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == py && axis == 1 && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == py && axis == 1 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == pz && axis == 2 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == pz && axis == 2 && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == px && axis == 0 && valenceOrbitalB == dzx) ){
      int anotherAxis;
      if(valenceOrbitalB == dxy){
         if(axis == 0){
            anotherAxis = 1;
         }
         else{
            anotherAxis = 0;
         }
      }
      else if(valenceOrbitalB == dyz){
         if(axis == 1){
            anotherAxis = 2;
         }
         else{
            anotherAxis = 1;
         }
      }
      else if(valenceOrbitalB == dzx){
         if(axis == 2){
            anotherAxis = 0;
         }
         else{
            anotherAxis = 2;
         }
      }
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*gaussianExponentA*dxyz[axis]
             -gaussianExponentB*dxyz[axis]
             +(gaussianExponentB*gaussianExponentB)*gaussianExponentA
             *dxyz[axis]*dxyz[axis]*dxyz[axis]/gauPlusAB;
      value *= 8.0*gaussianExponentA*sqrt(gaussianExponentA)
              *gaussianExponentB*dxyz[anotherAxis]*overlapSASB/(gauPlusAB*gauPlusAB*gauPlusAB);
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy && axis == 0 && valenceOrbitalB == px) || 
            (valenceOrbitalA == dxy && axis == 1 && valenceOrbitalB == py) ||
            (valenceOrbitalA == dyz && axis == 1 && valenceOrbitalB == py) ||
            (valenceOrbitalA == dyz && axis == 2 && valenceOrbitalB == pz) ||
            (valenceOrbitalA == dzx && axis == 2 && valenceOrbitalB == pz) ||
            (valenceOrbitalA == dzx && axis == 0 && valenceOrbitalB == px) ){
      int anotherAxis;
      if(valenceOrbitalA == dxy){
         if(axis == 0){
            anotherAxis = 1;
         }
         else{
            anotherAxis = 0;
         }
      }
      else if(valenceOrbitalA == dyz){
         if(axis == 1){
            anotherAxis = 2;
         }
         else{
            anotherAxis = 1;
         }
      }
      else if(valenceOrbitalA == dzx){
         if(axis == 2){
            anotherAxis = 0;
         }
         else{
            anotherAxis = 2;
         }
      }
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*gaussianExponentB*dxyz[axis]
             -gaussianExponentA*dxyz[axis]
             +(gaussianExponentA*gaussianExponentA)*gaussianExponentB
             *dxyz[axis]*dxyz[axis]*dxyz[axis]/gauPlusAB;
      value *= 8.0*gaussianExponentB*sqrt(gaussianExponentB)*gaussianExponentA
              *dxyz[anotherAxis]*overlapSASB/(gauPlusAB*gauPlusAB*gauPlusAB);
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == pz && axis == 2 && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == px && axis == 0 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == py && axis == 1 && valenceOrbitalB == dzx) ){
      int anotherAxis1;
      int anotherAxis2;
      if(axis == 0){
         anotherAxis1 = 1;
         anotherAxis2 = 2;
      }
      else if(axis == 1){
         anotherAxis1 = 0;
         anotherAxis2 = 2;
      }
      else if(axis == 2){
         anotherAxis1 = 0;
         anotherAxis2 = 1;
      }
      double overlapAOs1=0.0; 
      overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                valenceOrbitalA, 
                                                gaussianExponentA, 
                                                 
                                                valenceOrbitalB, 
                                                gaussianExponentB, 
                                                dxyz[0], 
                                                dxyz[1], 
                                                dxyz[2], 
                                                rAB,
                                                overlapSASB);
      value = 0.5+gaussianExponentB*gaussianExponentB*dxyz[axis]*dxyz[axis]/gauPlusAB;
      value *= 8.0*gaussianExponentA*gaussianExponentA*sqrt(gaussianExponentA)
              *gaussianExponentB*overlapSASB/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= dxyz[anotherAxis1]*dxyz[anotherAxis2];
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy && axis == 2 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dyz && axis == 0 && valenceOrbitalB == px) ||
            (valenceOrbitalA == dzx && axis == 1 && valenceOrbitalB == py) ){
      int anotherAxis1;
      int anotherAxis2;
      if(axis == 0){
         anotherAxis1 = 1;
         anotherAxis2 = 2;
      }
      else if(axis == 1){
         anotherAxis1 = 0;
         anotherAxis2 = 2;
      }
      else if(axis == 2){
         anotherAxis1 = 0;
         anotherAxis2 = 1;
      }
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB, 
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5+gaussianExponentA*gaussianExponentA*dxyz[axis]*dxyz[axis]/gauPlusAB;
      value *= 8.0*gaussianExponentB*gaussianExponentB*sqrt(gaussianExponentB)
              *gaussianExponentA*overlapSASB/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= dxyz[anotherAxis1]*dxyz[anotherAxis2];
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == 1 && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == px && axis == 2 && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == py && axis == 0 && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == py && axis == 2 && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == pz && axis == 0 && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == pz && axis == 1 && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == px && axis == 1 && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == px && axis == 2 && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == py && axis == 0 && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == py && axis == 2 && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == pz && axis == 0 && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == pz && axis == 1 && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == px && axis == 1 && valenceOrbitalB == dzx) || 
            (valenceOrbitalA == px && axis == 2 && valenceOrbitalB == dzx) || 
            (valenceOrbitalA == py && axis == 0 && valenceOrbitalB == dzx) || 
            (valenceOrbitalA == py && axis == 2 && valenceOrbitalB == dzx) || 
            (valenceOrbitalA == pz && axis == 0 && valenceOrbitalB == dzx) || 
            (valenceOrbitalA == pz && axis == 1 && valenceOrbitalB == dzx) || 
            (valenceOrbitalA == px && axis == 1 && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == px && axis == 2 && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == py && axis == 0 && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == py && axis == 2 && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == pz && axis == 0 && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == pz && axis == 1 && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == px && axis == 1 && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == px && axis == 2 && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == py && axis == 0 && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == py && axis == 2 && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == pz && axis == 0 && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == pz && axis == 1 && valenceOrbitalB == dzz) ){
      int dOrbital;
      if( (valenceOrbitalA == py && axis == 0) || 
          (valenceOrbitalA == px && axis == 1) ){
         dOrbital = dxy;
      }
      else if( (valenceOrbitalA == py && axis == 2) || 
               (valenceOrbitalA == pz && axis == 1) ){
         dOrbital = dyz;
      }
      else if( (valenceOrbitalA == px && axis == 2) || 
               (valenceOrbitalA == pz && axis == 0) ){
         dOrbital = dzx;
      }
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      double overlapAOs2 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       dOrbital, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = overlapAOs2/(2.0*sqrt(gaussianExponentA))+xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy   && axis == 1 && valenceOrbitalB == px) || 
            (valenceOrbitalA == dxy   && axis == 2 && valenceOrbitalB == px) || 
            (valenceOrbitalA == dxy   && axis == 0 && valenceOrbitalB == py) || 
            (valenceOrbitalA == dxy   && axis == 2 && valenceOrbitalB == py) || 
            (valenceOrbitalA == dxy   && axis == 0 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dxy   && axis == 1 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dyz   && axis == 1 && valenceOrbitalB == px) || 
            (valenceOrbitalA == dyz   && axis == 2 && valenceOrbitalB == px) || 
            (valenceOrbitalA == dyz   && axis == 0 && valenceOrbitalB == py) || 
            (valenceOrbitalA == dyz   && axis == 2 && valenceOrbitalB == py) || 
            (valenceOrbitalA == dyz   && axis == 0 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dyz   && axis == 1 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dzx   && axis == 1 && valenceOrbitalB == px) || 
            (valenceOrbitalA == dzx   && axis == 2 && valenceOrbitalB == px) || 
            (valenceOrbitalA == dzx   && axis == 0 && valenceOrbitalB == py) || 
            (valenceOrbitalA == dzx   && axis == 2 && valenceOrbitalB == py) || 
            (valenceOrbitalA == dzx   && axis == 0 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dzx   && axis == 1 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dxxyy && axis == 1 && valenceOrbitalB == px) || 
            (valenceOrbitalA == dxxyy && axis == 2 && valenceOrbitalB == px) || 
            (valenceOrbitalA == dxxyy && axis == 0 && valenceOrbitalB == py) || 
            (valenceOrbitalA == dxxyy && axis == 2 && valenceOrbitalB == py) || 
            (valenceOrbitalA == dxxyy && axis == 0 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dxxyy && axis == 1 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dzz   && axis == 1 && valenceOrbitalB == px) || 
            (valenceOrbitalA == dzz   && axis == 2 && valenceOrbitalB == px) || 
            (valenceOrbitalA == dzz   && axis == 0 && valenceOrbitalB == py) || 
            (valenceOrbitalA == dzz   && axis == 2 && valenceOrbitalB == py) || 
            (valenceOrbitalA == dzz   && axis == 0 && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dzz   && axis == 1 && valenceOrbitalB == pz) ){
      int dOrbital;
      if( (valenceOrbitalB == py && axis == 0) || 
          (valenceOrbitalB == px && axis == 1) ){
         dOrbital = dxy;
      }
      else if( (valenceOrbitalB == py && axis == 2) || 
               (valenceOrbitalB == pz && axis == 1) ){
         dOrbital = dyz;
      }
      else if( (valenceOrbitalB == px && axis == 2) || 
               (valenceOrbitalB == pz && axis == 0) ){
         dOrbital = dzx;
      }
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      double overlapAOs2 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       dOrbital, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = overlapAOs2/(2.0*sqrt(gaussianExponentB))+xyzB[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == px && axis == 0 && valenceOrbitalB == dxxyy){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5-2.0*gauMultAB*(dxyz[0]*dxyz[0])/gauPlusAB;
      value += 0.5*(gaussianExponentA*gaussianExponentA)*((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/gauPlusAB;
      value += (gauMultAB*dxyz[0]
               *gauMultAB*dxyz[0])
              *((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/(gauPlusAB*gauPlusAB);
      value *= 4.0*sqrt(gaussianExponentA)*gaussianExponentB/(gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == 0 && valenceOrbitalB == px){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5-2.0*gauMultAB*(dxyz[0]*dxyz[0])/gauPlusAB;
      value += 0.5*(gaussianExponentB*gaussianExponentB)*((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/gauPlusAB;
      value += (gauMultAB*dxyz[0]
               *gauMultAB*dxyz[0])
              *((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/(gauPlusAB*gauPlusAB);
      value *= 4.0*sqrt(gaussianExponentB)*gaussianExponentA/(gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == py && axis == 1 && valenceOrbitalB == dxxyy){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5-2.0*gauMultAB*(dxyz[1]*dxyz[1])/gauPlusAB;
      value += 0.5*(gaussianExponentA*gaussianExponentA)*((dxyz[1]*dxyz[1])-(dxyz[0]*dxyz[0]))/gauPlusAB;
      value += (gauMultAB*dxyz[1]
               *gauMultAB*dxyz[1])
              *((dxyz[1]*dxyz[1])-(dxyz[0]*dxyz[0]))/(gauPlusAB*gauPlusAB);
      value *= -4.0*sqrt(gaussianExponentA)*gaussianExponentB/(gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == 1 && valenceOrbitalB == py){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5-2.0*gauMultAB*(dxyz[1]*dxyz[1])/gauPlusAB;
      value += 0.5*(gaussianExponentB*gaussianExponentB)*((dxyz[1]*dxyz[1])-(dxyz[0]*dxyz[0]))/gauPlusAB;
      value += (gauMultAB*dxyz[1]
               *gauMultAB*dxyz[1])
              *((dxyz[1]*dxyz[1])-(dxyz[0]*dxyz[0]))/(gauPlusAB*gauPlusAB);
      value *= -4.0*sqrt(gaussianExponentB)*gaussianExponentA/(gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == pz && axis == 2 && valenceOrbitalB == dxxyy){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]));
      value += gaussianExponentB*gaussianExponentB*dxyz[2]*dxyz[2]
              *((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))
              /gauPlusAB;
      value *= 4.0*gaussianExponentA*gaussianExponentA*sqrt(gaussianExponentA)
              *gaussianExponentB/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == 2 && valenceOrbitalB == pz){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]));
      value += gaussianExponentA*gaussianExponentA*dxyz[2]*dxyz[2]
              *((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))
              /gauPlusAB;
      value *= 4.0*gaussianExponentB*gaussianExponentB*sqrt(gaussianExponentB)
              *gaussianExponentA/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == 0 && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == py && axis == 1 && valenceOrbitalB == dzz) ){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5
             +2.0*gauMultAB*(dxyz[axis]*dxyz[axis])/gauPlusAB
             +0.5*(gaussianExponentA*gaussianExponentA)
                 *(2.0*(dxyz[2]*dxyz[2])-(dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/gauPlusAB
             +(gauMultAB*gauMultAB*dxyz[axis]*dxyz[axis]/(gauPlusAB*gauPlusAB))
              *(2.0*(dxyz[2]*dxyz[2])-(dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]));
      value *= 4.0*sqrt(gaussianExponentA)*gaussianExponentB/(gauPlusAB*gauPlusAB*sqrt(3.0));
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dzz && axis == 0 && valenceOrbitalB == px) || 
            (valenceOrbitalA == dzz && axis == 1 && valenceOrbitalB == py) ){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5
             +2.0*gauMultAB*(dxyz[axis]*dxyz[axis])/gauPlusAB
             +0.5*(gaussianExponentB*gaussianExponentB)
                 *(2.0*(dxyz[2]*dxyz[2])-(dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/gauPlusAB
             +(gauMultAB*gauMultAB*dxyz[axis]*dxyz[axis]/(gauPlusAB*gauPlusAB))
              *(2.0*(dxyz[2]*dxyz[2])-(dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]));
      value *= 4.0*sqrt(gaussianExponentB)*gaussianExponentA/(gauPlusAB*gauPlusAB*sqrt(3.0));
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == pz && axis == 2 && valenceOrbitalB == dzz){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 1.0
             -4.0*gauMultAB*(dxyz[axis]*dxyz[axis])/gauPlusAB
             +0.5*(gaussianExponentA*gaussianExponentA)
                 *(2.0*(dxyz[2]*dxyz[2])-(dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/gauPlusAB
             +(gauMultAB*gauMultAB*dxyz[axis]*dxyz[axis]/(gauPlusAB*gauPlusAB))
              *(2.0*(dxyz[2]*dxyz[2])-(dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]));
      value *= 4.0*sqrt(gaussianExponentA)*gaussianExponentB/(gauPlusAB*gauPlusAB*sqrt(3.0));
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dzz && axis == 2 && valenceOrbitalB == pz){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 1.0
             -4.0*gauMultAB*(dxyz[axis]*dxyz[axis])/gauPlusAB
             +0.5*(gaussianExponentB*gaussianExponentB)
                 *(2.0*(dxyz[2]*dxyz[2])-(dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/gauPlusAB
             +(gauMultAB*gauMultAB*dxyz[axis]*dxyz[axis]/(gauPlusAB*gauPlusAB))
              *(2.0*(dxyz[2]*dxyz[2])-(dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]));
      value *= 4.0*sqrt(gaussianExponentB)*gaussianExponentA/(gauPlusAB*gauPlusAB*sqrt(3.0));
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy && axis == 0 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dxy && axis == 1 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dyz && axis == 1 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dyz && axis == 2 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dzx && axis == 2 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dzx && axis == 0 && valenceOrbitalB == valenceOrbitalA) ){
      int anotherAxis;
      if(valenceOrbitalB == dxy){
         if(axis == 0){
            anotherAxis = 1;
         }
         else{
            anotherAxis = 0;
         }
      }
      else if(valenceOrbitalB == dyz){
         if(axis == 1){
            anotherAxis = 2;
         }
         else{
            anotherAxis = 1;
         }
      }
      else if(valenceOrbitalB == dzx){
         if(axis == 2){
            anotherAxis = 0;
         }
         else{
            anotherAxis = 2;
         }
      }
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      double gauMinsAB = gaussianExponentA-gaussianExponentB;
      value = 0.5*gauMinsAB*dxyz[axis]
             +gauMinsAB*gauMultAB
              *dxyz[axis]*(dxyz[anotherAxis]*dxyz[anotherAxis])/gauPlusAB;
      value *= 8.0*gauMultAB/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy   && axis == 2 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dyz   && axis == 0 && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dzx   && axis == 1 && valenceOrbitalB == valenceOrbitalA) ||
            (valenceOrbitalA == dxxyy && axis == 2 && valenceOrbitalB == valenceOrbitalA) ){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxxyy && axis == 0 && valenceOrbitalB == valenceOrbitalA) ||
            (valenceOrbitalA == dxxyy && axis == 1 && valenceOrbitalB == valenceOrbitalA) ){
      int anotherAxis;
      if(axis == 0){
         anotherAxis = 1;
      }
      else{
         anotherAxis = 0;
      }
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      double gauMinsAB = gaussianExponentA-gaussianExponentB;
      value = gauMinsAB*dxyz[axis]
             -gauMinsAB*gauMultAB
              *((dxyz[axis]*dxyz[axis]) - (dxyz[anotherAxis]*dxyz[anotherAxis]))*dxyz[axis]/gauPlusAB;
      value *= 4.0*gauMultAB/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dzz && axis == 0 && valenceOrbitalB == valenceOrbitalA) ||
            (valenceOrbitalA == dzz && axis == 1 && valenceOrbitalB == valenceOrbitalA) ){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      double gauMinsAB = gaussianExponentA-gaussianExponentB;
      value = gauMinsAB*dxyz[axis]
             -gauMinsAB*gauMultAB
              *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))*dxyz[axis]/gauPlusAB;
      value *= 4.0*gauMultAB/(3.0*gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( valenceOrbitalA == dzz && axis == 2 && valenceOrbitalB == valenceOrbitalA) {
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      double gauMinsAB = gaussianExponentA-gaussianExponentB;
      value = 2.0*gauMinsAB*dxyz[axis]
             -gauMinsAB*gauMultAB
              *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))*dxyz[axis]/gauPlusAB;
      value *= 8.0*gauMultAB/(3.0*gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy && axis == 1 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dyz && axis == 2 && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == dzx && axis == 0 && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == dyz && axis == 1 && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == dzx && axis == 2 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dxy && axis == 0 && valenceOrbitalB == dzx) ){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      double gauMinsAB = gaussianExponentA-gaussianExponentB;
      value = -8.0*(gaussianExponentA*gauMultAB*gaussianExponentB)*dxyz[0]*dxyz[1]*dxyz[2]
             *gauMinsAB/(gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy && axis == 0 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dyz && axis == 1 && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == dzx && axis == 2 && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == dyz && axis == 2 && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == dzx && axis == 0 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dxy && axis == 1 && valenceOrbitalB == dzx) ){
      int anotherAxis1;
      int anotherAxis2;
      if(valenceOrbitalA == dxy && valenceOrbitalB == dyz){
         anotherAxis1 = 1;
         anotherAxis2 = 2;
      }
      else if(valenceOrbitalA == dyz && valenceOrbitalB == dzx){
         anotherAxis1 = 2;
         anotherAxis2 = 0;
      }
      else if(valenceOrbitalA == dzx && valenceOrbitalB == dxy){
         anotherAxis1 = 0;
         anotherAxis2 = 1;
      }
      else if(valenceOrbitalA == dyz && valenceOrbitalB == dxy){
         anotherAxis1 = 1;
         anotherAxis2 = 0;
      }
      else if(valenceOrbitalA == dzx && valenceOrbitalB == dyz){
         anotherAxis1 = 2;
         anotherAxis2 = 1;
      }
      else if(valenceOrbitalA == dxy && valenceOrbitalB == dzx){
         anotherAxis1 = 0;
         anotherAxis2 = 2;
      }
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5-gauMultAB*dxyz[anotherAxis1]*dxyz[anotherAxis1]/gauPlusAB;
      value *= 8.0*(gaussianExponentA*gaussianExponentA)*gaussianExponentB*dxyz[anotherAxis2]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dyz && axis == 0 && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == dzx && axis == 1 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dxy && axis == 2 && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == dxy && axis == 2 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dyz && axis == 0 && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == dzx && axis == 1 && valenceOrbitalB == dxy) ){
      int anotherAxis1;
      int anotherAxis2;
      if(valenceOrbitalA == dyz && valenceOrbitalB == dxy){
         anotherAxis1 = 1;
         anotherAxis2 = 2;
      }
      else if(valenceOrbitalA == dzx && valenceOrbitalB == dyz){
         anotherAxis1 = 2;
         anotherAxis2 = 0;
      }
      else if(valenceOrbitalA == dxy && valenceOrbitalB == dzx){
         anotherAxis1 = 0;
         anotherAxis2 = 1;
      }
      else if(valenceOrbitalA == dxy && valenceOrbitalB == dyz){
         anotherAxis1 = 1;
         anotherAxis2 = 0;
      }
      else if(valenceOrbitalA == dyz && valenceOrbitalB == dzx){
         anotherAxis1 = 2;
         anotherAxis2 = 1;
      }
      else if(valenceOrbitalA == dzx && valenceOrbitalB == dxy){
         anotherAxis1 = 0;
         anotherAxis2 = 2;
      }

      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5-gauMultAB*dxyz[anotherAxis1]*dxyz[anotherAxis1]/gauPlusAB;
      value *= -8.0*(gaussianExponentB*gaussianExponentB)*gaussianExponentA*dxyz[anotherAxis2]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == 0 && valenceOrbitalB == dxy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*gauPlusAB
             -(gaussianExponentA*gaussianExponentA)*gaussianExponentB*(dxyz[0]*dxyz[0])/gauPlusAB
             +(gaussianExponentB*gaussianExponentB)*gaussianExponentA
             *((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/(2.0*gauPlusAB);
      value *= 8.0*gauMultAB*dxyz[1]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxy && axis == 0 && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*gauPlusAB
             -(gaussianExponentB*gaussianExponentB)*gaussianExponentA*(dxyz[0]*dxyz[0])/gauPlusAB
             +(gaussianExponentA*gaussianExponentA)*gaussianExponentB
             *((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/(2.0*gauPlusAB);
      value *= -8.0*gauMultAB*dxyz[1]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   } 
   else if(valenceOrbitalA == dxxyy && axis == 1 && valenceOrbitalB == dxy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5*gauPlusAB
             +(gaussianExponentA*gaussianExponentA)*gaussianExponentB*(dxyz[1]*dxyz[1])/gauPlusAB
             +(gaussianExponentB*gaussianExponentB)*gaussianExponentA
             *((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/(2.0*gauPlusAB);
      value *= 8.0*gauMultAB*dxyz[0]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxy && axis == 1 && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5*gauPlusAB
             +(gaussianExponentB*gaussianExponentB)*gaussianExponentA*(dxyz[1]*dxyz[1])/gauPlusAB
             +(gaussianExponentA*gaussianExponentA)*gaussianExponentB
             *((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/(2.0*gauPlusAB);
      value *= -8.0*gauMultAB*dxyz[0]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxxyy && axis == 0 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dxxyy && axis == 1 && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == dxxyy && axis == 2 && valenceOrbitalB == dxy) ){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 8.0*(gaussianExponentA*gaussianExponentA*gaussianExponentA*gaussianExponentA)
             *(gaussianExponentB*gaussianExponentB*gaussianExponentB)
             *dxyz[0]*dxyz[1]*dxyz[2]
             *((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))
             /(gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dyz && axis == 0 && valenceOrbitalB == dxxyy) ||
            (valenceOrbitalA == dzx && axis == 1 && valenceOrbitalB == dxxyy) ||
            (valenceOrbitalA == dxy && axis == 2 && valenceOrbitalB == dxxyy) ){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = -8.0*(gaussianExponentA*gaussianExponentA*gaussianExponentA)
             *(gaussianExponentB*gaussianExponentB*gaussianExponentB*gaussianExponentB)
             *dxyz[0]*dxyz[1]*dxyz[2]
             *((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))
             /(gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == 1 && valenceOrbitalB == dyz){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5*gaussianExponentA
             +(gaussianExponentA*gaussianExponentA)*gaussianExponentB*(dxyz[1]*dxyz[1])/gauPlusAB
             +(gaussianExponentB*gaussianExponentB)*gaussianExponentA
             *((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/(2.0*gauPlusAB);
      value *= 8.0*gauMultAB*dxyz[2]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dyz && axis == 1 && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5*gaussianExponentB
             +(gaussianExponentB*gaussianExponentB)*gaussianExponentA*(dxyz[1]*dxyz[1])/gauPlusAB
             +(gaussianExponentA*gaussianExponentA)*gaussianExponentB
             *((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/(2.0*gauPlusAB);
      value *= -8.0*gauMultAB*dxyz[2]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == 0 && valenceOrbitalB == dzx){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5*gaussianExponentA
             +(gaussianExponentA*gaussianExponentA)*gaussianExponentB*(dxyz[0]*dxyz[0])/gauPlusAB
             +(gaussianExponentB*gaussianExponentB)*gaussianExponentA
             *((dxyz[1]*dxyz[1])-(dxyz[0]*dxyz[0]))/(2.0*gauPlusAB);
      value *= -8.0*gauMultAB*dxyz[2]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dzx && axis == 0 && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5*gaussianExponentB
             +(gaussianExponentB*gaussianExponentB)*gaussianExponentA*(dxyz[1]*dxyz[1])/gauPlusAB
             +(gaussianExponentA*gaussianExponentA)*gaussianExponentB
             *((dxyz[1]*dxyz[1])-(dxyz[0]*dxyz[0]))/(2.0*gauPlusAB);
      value *= 8.0*gauMultAB*dxyz[2]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == 2 && valenceOrbitalB == dyz){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = gauMultAB*((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/gauPlusAB+1.0;
      value *= 4.0*gaussianExponentA*(gaussianExponentB*gaussianExponentB)*dxyz[1]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dyz && axis == 2 && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = gauMultAB*((dxyz[0]*dxyz[0])-(dxyz[1]*dxyz[1]))/gauPlusAB+1.0;
      value *= -4.0*gaussianExponentB*(gaussianExponentA*gaussianExponentA)*dxyz[1]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == 2 && valenceOrbitalB == dzx){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = gauMultAB*((dxyz[1]*dxyz[1])-(dxyz[0]*dxyz[0]))/gauPlusAB+1.0;
      value *= -4.0*gaussianExponentA*(gaussianExponentB*gaussianExponentB)*dxyz[0]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dzx && axis == 2 && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = gauMultAB*((dxyz[1]*dxyz[1])-(dxyz[0]*dxyz[0]))/gauPlusAB+1.0;
      value *= 4.0*gaussianExponentB*(gaussianExponentA*gaussianExponentA)*dxyz[0]/(gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dzz && axis == 0 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dzz && axis == 1 && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == dzz && axis == 2 && valenceOrbitalB == dxy) ){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]);
      value *= dxyz[0]*dxyz[1]*dxyz[2];
      value *= 8.0*(gaussianExponentA*gaussianExponentA*gaussianExponentA*gaussianExponentA)
              *(gaussianExponentB*gaussianExponentB*gaussianExponentB);
      value /= sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dyz && axis == 0 && valenceOrbitalB == dzz) ||
            (valenceOrbitalA == dzx && axis == 1 && valenceOrbitalB == dzz) ||
            (valenceOrbitalA == dxy && axis == 2 && valenceOrbitalB == dzz) ){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]);
      value *= dxyz[0]*dxyz[1]*dxyz[2];
      value *= -8.0*(gaussianExponentB*gaussianExponentB*gaussianExponentB*gaussianExponentB)*(gaussianExponentA*gaussianExponentA*gaussianExponentA);
      value /= sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dzz && axis == 0 && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == dzz && axis == 1 && valenceOrbitalB == dxy) ){
      int anotherAxis;
      if(axis == 0){
         anotherAxis = 1;
      }
      else if(axis == 1){
         anotherAxis = 0;
      }
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*(gaussianExponentB-gaussianExponentA)
             +3.0*(gaussianExponentA*gaussianExponentA)*gaussianExponentB
             *(dxyz[axis]*dxyz[axis])/(gauPlusAB*gauPlusAB)
             +gaussianExponentA*(gaussianExponentB*gaussianExponentB)
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/(2.0*gauPlusAB)
             +(gaussianExponentA*gaussianExponentA*gaussianExponentA)*(gaussianExponentB*gaussianExponentB)*(dxyz[axis]*dxyz[axis])
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))
             /(gauPlusAB*gauPlusAB);
      value *= 8.0*gauMultAB*dxyz[anotherAxis]/(sqrt(3.0)*gauPlusAB*gauPlusAB*gauPlusAB);
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy && axis == 0 && valenceOrbitalB == dzz) ||
            (valenceOrbitalA == dxy && axis == 1 && valenceOrbitalB == dzz) ){
      int anotherAxis;
      if(axis == 0){
         anotherAxis = 1;
      }
      else if(axis == 1){
         anotherAxis = 0;
      }
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      double gauMinsAB = gaussianExponentA-gaussianExponentB;
      value = 0.5*gauMinsAB
             +3.0*(gaussianExponentB*gaussianExponentB)*gaussianExponentA
             *(dxyz[axis]*dxyz[axis])/(gauPlusAB*gauPlusAB)
             +gaussianExponentB*(gaussianExponentA*gaussianExponentA)
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/(2.0*gauPlusAB)
             +(gaussianExponentB*gaussianExponentB*gaussianExponentB)*(gaussianExponentA*gaussianExponentA)*(dxyz[axis]*dxyz[axis])
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/(gauPlusAB*gauPlusAB);
      value *= -8.0*gauMultAB*dxyz[anotherAxis]/(sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB));
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dzz && axis == 1 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dzz && axis == 0 && valenceOrbitalB == dzx) ){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentB-0.5*gaussianExponentA
             +gaussianExponentA*(gaussianExponentB*gaussianExponentB)
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/(2.0*gauPlusAB)
             +(gaussianExponentA*gaussianExponentA*gaussianExponentA)
             *(gaussianExponentB*gaussianExponentB*dxyz[axis]*dxyz[axis])
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))
             /(gauPlusAB*gauPlusAB);
      value *= 8.0*gauMultAB*dxyz[2]
              /(sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB));
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dyz && axis == 1 && valenceOrbitalB == dzz) ||
            (valenceOrbitalA == dzx && axis == 0 && valenceOrbitalB == dzz) ){
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentA-0.5*gaussianExponentB
             +gaussianExponentB*(gaussianExponentA*gaussianExponentA)
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/(2.0*gauPlusAB)
             +(gaussianExponentB*gaussianExponentB*gaussianExponentB)*(gaussianExponentA*gaussianExponentA*dxyz[axis]*dxyz[axis])
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))
             /(gauPlusAB*gauPlusAB);
      value *= -8.0*gauMultAB*dxyz[2]
              /(sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB));
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dzz && axis == 2 && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dzz && axis == 2 && valenceOrbitalB == dzx) ){
      int anotherAxis;
      if(valenceOrbitalB == dyz){
         anotherAxis = 1;
      }
      else if(valenceOrbitalB == dzx){
         anotherAxis = 0;
      }
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentA-0.5*gaussianExponentB
             -3.0*(gaussianExponentA*gaussianExponentA)*gaussianExponentB
             *dxyz[axis]*dxyz[axis]/(gauPlusAB*gauPlusAB)
             +gaussianExponentA*(gaussianExponentB*gaussianExponentB)
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/(2.0*gauPlusAB)
             +(gaussianExponentA*gaussianExponentA*gaussianExponentA)
             *(gaussianExponentB*gaussianExponentB*dxyz[axis]*dxyz[axis])
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))
             /(gauPlusAB*gauPlusAB);
      value *= 8.0*gauMultAB*dxyz[anotherAxis]
              /(sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB));
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dyz && axis == 2 && valenceOrbitalB == dzz) ||
            (valenceOrbitalA == dzx && axis == 2 && valenceOrbitalB == dzz) ){
      int anotherAxis;
      if(valenceOrbitalA == dyz){
         anotherAxis = 1;
      }
      else if(valenceOrbitalA == dzx){
         anotherAxis = 0;
      }
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentB-0.5*gaussianExponentA
             -3.0*(gaussianExponentB*gaussianExponentB)*gaussianExponentA*dxyz[axis]*dxyz[axis]/(gauPlusAB*gauPlusAB)
             +gaussianExponentB*(gaussianExponentA*gaussianExponentA)
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/(2.0*gauPlusAB)
             +(gaussianExponentB*gaussianExponentB*gaussianExponentB)
             *(gaussianExponentA*gaussianExponentA*dxyz[axis]*dxyz[axis])
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))
             /(gauPlusAB*gauPlusAB);
      value *= -8.0*gauMultAB*dxyz[anotherAxis]
              /(sqrt(3.0)*(gauPlusAB*gauPlusAB*gauPlusAB));
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dzz && axis == 0 && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentB - gaussianExponentA
             +gaussianExponentA*(gaussianExponentB*gaussianExponentB)
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/gauPlusAB
             +(gaussianExponentA*gaussianExponentA)*gaussianExponentB
             *((dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/gauPlusAB;
      value *= 4.0*gauMultAB*dxyz[0]/(gauPlusAB*gauPlusAB*gauPlusAB*sqrt(3.0));
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == 0 && valenceOrbitalB == dzz){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentA - gaussianExponentB
             +gaussianExponentB*(gaussianExponentA*gaussianExponentA)
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/gauPlusAB
             +(gaussianExponentB*gaussianExponentB)*gaussianExponentA
             *((dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/gauPlusAB;
      value *= -4.0*gauMultAB*dxyz[0]/(gauPlusAB*gauPlusAB*gauPlusAB*sqrt(3.0));
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dzz && axis == 1 && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentB - gaussianExponentA
             +gaussianExponentA*(gaussianExponentB*gaussianExponentB)
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/gauPlusAB
             +(gaussianExponentA*gaussianExponentA)*gaussianExponentB
             *((dxyz[1]*dxyz[1]) - (dxyz[0]*dxyz[0]))/gauPlusAB;
      value *= -4.0*gauMultAB*dxyz[1]/(gauPlusAB*gauPlusAB*gauPlusAB*sqrt(3.0));
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == 1 && valenceOrbitalB == dzz){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentA - gaussianExponentB
             +gaussianExponentB*(gaussianExponentA*gaussianExponentA)
             *(2.0*(dxyz[2]*dxyz[2]) - (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]))/gauPlusAB
             +(gaussianExponentB*gaussianExponentB)*gaussianExponentA
             *((dxyz[1]*dxyz[1]) - (dxyz[0]*dxyz[0]))/gauPlusAB;
      value *= 4.0*gauMultAB*dxyz[1]/(gauPlusAB*gauPlusAB*gauPlusAB*sqrt(3.0));
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dzz && axis == 2 && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]);
      value *= -8.0*(gaussianExponentA*gaussianExponentA*gaussianExponentA)*(gaussianExponentB*gaussianExponentB)*dxyz[2];
      value /= sqrt(3.0)*gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB;
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == 2 && valenceOrbitalB == dzz){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/gauPlusAB;
      double overlapAOs1 = Cpp_Cndo2GetGaussianOverlapAOs(
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                        
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[0], 
                                                       dxyz[1], 
                                                       dxyz[2], 
                                                       rAB,
                                                       overlapSASB);
      value = (dxyz[0]*dxyz[0]) - (dxyz[1]*dxyz[1]);
      value *= 8.0*(gaussianExponentA*gaussianExponentA)*(gaussianExponentB*gaussianExponentB*gaussianExponentB)*dxyz[2];
      value /= sqrt(3.0)*gauPlusAB*gauPlusAB*gauPlusAB*gauPlusAB;
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else{
      mexErrMsgTxt("Cpp_Cndo2GetGaussianCartesianMatrix: Orbital type wrong.");
   }

   return value;
}

void Cpp_Cndo2CalcCartesianMatrixElementsByGTOExpansionLoop( 
    double& xComponent, double& yComponent, double& zComponent, 
    int shellTypeA, int valenceOrbitalA, double orbitalExponentA, double const* xyzA, 
    int shellTypeB, int valenceOrbitalB, double orbitalExponentB, double const* xyzB, 
    double rAB, int stonG) {
        
    xComponent=0.0;
    yComponent=0.0;
    zComponent=0.0;
    double gaussianExponentA = 0.0;
    double gaussianExponentB = 0.0;
   double overlapSASB = 0.0;
   double temp  = 0.0;
   double tempX = 0.0;
   double tempY = 0.0;
   double tempZ = 0.0;
   for(int i=0; i<=stonG; i++){
      for(int j=0; j<=stonG; j++){
         temp = GTOExpansionSTO::GetInstance()->GetCoefficient(stonG, 
                                                               shellTypeA, 
                                                               valenceOrbitalA, 
                                                               i); 
         temp *= GTOExpansionSTO::GetInstance()->GetCoefficient(stonG, 
                                                                shellTypeB, 
                                                                valenceOrbitalB, 
                                                                j); 
         gaussianExponentA = (orbitalExponentA*orbitalExponentA) *
                             GTOExpansionSTO::GetInstance()->GetExponent(stonG, 
                                                                         shellTypeA, 
                                                                         valenceOrbitalA, 
                                                                         i);
         gaussianExponentB = (orbitalExponentB*orbitalExponentB) *
                             GTOExpansionSTO::GetInstance()->GetExponent(stonG, 
                                                                         shellTypeB, 
                                                                         valenceOrbitalB, 
                                                                         j);
         overlapSASB = Cpp_Cndo2GetGaussianOverlapAOsSASB(gaussianExponentA, gaussianExponentB, rAB);
         tempX = Cpp_Cndo2GetGaussianCartesianMatrix(valenceOrbitalA, gaussianExponentA, xyzA,
                                                  valenceOrbitalB, gaussianExponentB, xyzB, 
                                                  rAB, overlapSASB,
                                                  0);
         tempY = Cpp_Cndo2GetGaussianCartesianMatrix(valenceOrbitalA, gaussianExponentA, xyzA,
                                                  valenceOrbitalB, gaussianExponentB, xyzB, 
                                                  rAB, overlapSASB, 
                                                  1);
         tempZ = Cpp_Cndo2GetGaussianCartesianMatrix(valenceOrbitalA, gaussianExponentA, xyzA,
                                                  valenceOrbitalB, gaussianExponentB, xyzB, 
                                                  rAB, overlapSASB,
                                                  2);
         xComponent += temp*tempX;
         yComponent += temp*tempY;
         zComponent += temp*tempZ;
      }
   }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nlhs != 1)
        mexErrMsgTxt("Cpp_Cndo2CalcCartesianMatrixElementsByGTOExpansionLoop: 1 output expected.");
    if (nrhs != 10 || mxGetM(prhs[3]) != 3 || mxGetN(prhs[3]) != 1 || mxGetM(prhs[7]) != 3 || mxGetN(prhs[7]) != 1)
        mexErrMsgTxt("Cpp_Cndo2CalcCartesianMatrixElementsByGTOExpansionLoop: 8 input expected.");
    int shellTypeA = (int)mxGetScalar(prhs[0]) - 1;
    int valenceOrbitalA = (int)mxGetScalar(prhs[1]) - 1;
    double orbitalExponentA = mxGetScalar(prhs[2]);
    double const* xyzA = mxGetPr(prhs[3]);
    int shellTypeB = (int)mxGetScalar(prhs[4]) - 1;
    int valenceOrbitalB = (int)mxGetScalar(prhs[5]) - 1;
    double orbitalExponentB = mxGetScalar(prhs[6]);
    double const* xyzB = mxGetPr(prhs[7]);
    double rAB = mxGetScalar(prhs[8]);
    int stonG = (int)mxGetScalar(prhs[9]) - 1;
    
    double xComponent;
    double yComponent;
    double zComponent;
    Cpp_Cndo2CalcCartesianMatrixElementsByGTOExpansionLoop(
       xComponent, yComponent, zComponent, 
       shellTypeA, valenceOrbitalA, orbitalExponentA, xyzA, 
       shellTypeB, valenceOrbitalB, orbitalExponentB, xyzB, 
       rAB, stonG);
    plhs[0] = mxCreateDoubleMatrix( 3, 1, mxREAL);
    double* outputPtr = mxGetPr(plhs[0]);
    outputPtr[0] = xComponent;
    outputPtr[1] = yComponent;
    outputPtr[2] = zComponent;
}
