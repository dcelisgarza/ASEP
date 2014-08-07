#ifndef ADD_H
#define ADD_H
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
float ran3(long *idum); 
 
#endif 

float ran3(long *idum)
//Returns a uniform random deviate between 0.0 and 1.0. Set idum to any negative value to
//initialize or reinitialize the sequence.
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  
  if (*idum < 0 || iff == 0) {  //Initialization.
    iff=1;
    mj=MSEED- (*idum < 0 ?-*idum : *idum);
    mj%= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++){    //Now initialize the rest of the table, in a slightly random order, with numbers that are espacially random.
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk< MZ) mk +=MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)  // We randomize them by "warming up the generator."
      for (i=1;i<=55;i++){ 
	ma[i]-= ma[1+(i+30) % 55];
	if (ma[i]< MZ) ma[i] += MBIG;
      }
    inext=0;       // Prepare indices for our first generated number.
    inextp=31;      // The constant 31 is special; see Knuth.
    *idum=1;
  }
  // Here is where we start, except on initialization.
  if (++inext== 56) inext=1;  // Increment inext and inextp, wrapping around 56 to 1.
  if (++inextp== 56) inextp=1;  // Generate a random number subtractively.
  mj=ma[inext]-ma[inextp]; // Be sure that it is in range.
  if (mj< MZ) mj+=MBIG;  
  ma[inext]=mj;                 // store it;
  return mj*FAC;        // and output the derived uniform deviate.
} 
