/*
	Programmer: Chris Mirkes
	Code optimized by Ali Aghaeifar and modified to support high dimention data
*/


#include "mex.h"
#include <stdio.h>
#include "matrix.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>


#include <random>

// uncomment following 3 lines for windows OS 
/*
#include <windows.h>
#include <ppl.h>
using namespace concurrency;
 */
using namespace std;

#define pi 3.14159265358979323846
/*
int round( double r ) {
    return (r > 0.0) ? (r + 0.5) : (r - 0.5); 
}
*/

void linspace(double *Array, double d1, double d2, long n){
    double Increment;
    Increment = (d2-d1)/((double)(n-1)); 
	
    for (int i = 0; i < n-1; i++)
        Array[i] = d1+ i*Increment;
    
    Array[n-1] = d2;
}

void linspace(double *Array, long d1, long d2, long n){
    double Increment;
    Increment = (d2-d1)/((double)(n-1));
    
    for (int i = 0; i < n-1; i++)
        Array[i] = d1+ i*Increment;

    Array[n-1] = d2;
}



void InterpolateToSTRFoV(double *dImage,double *X, double *Y, double *Z,double dFOVST, long lNST, double dPixelSizeOrig, long NX, long NY, long NZ, double *dImageST){
    
    // Input check
//    if(!(lNST%2))
//        mexErrMsgTxt("PixelSizeST must be odd!\n");   
    
    // Square of pixel number of standard image
    long lNSTsq = lNST*lNST;
    long lNSTqu = lNSTsq*lNST;
    long lNOrigIndices = NX*NY*NZ;
    
    double dPixelSize = dFOVST/lNST; // dFOVST/(lNST-1.0); //aa:
    long lCentreSTR   = (lNST+1)/2;  // (lNST+1)/2-1;
    
    // The kernel is a box in which we look for the nearest given point on the pre-gridded image
    long lKernelsize   = (long)ceil(dPixelSizeOrig/dPixelSize) + ((long)ceil(dPixelSizeOrig/dPixelSize))%2+1;
    long lKernelExtend = (lKernelsize-1)/2;
    long lKernelsizesq = lKernelsize*lKernelsize;
    
	// mexPrintf("Size = %f %f %d\n",dPixelSizeOrig, dPixelSize, lKernelsize);
    // Tells if the pixel has been assigned a value in the pre-gridded image
    long *lBinaryMask = new long[lNSTqu];    
    double *dImageSTtemp = new double[lNSTqu];
        
    // Initialize matrices
	memset(lBinaryMask	, 0, lNSTqu*sizeof(long));
	memset(dImageST		, 0, lNSTqu*sizeof(double));
	memset(dImageSTtemp	, 0, lNSTqu*sizeof(double));
	
   // Create search cube to look for lBinaryMask[Index] == 1 around a given voxel
   // The cubes XYZ store the relative coordinates
   long lNCube = lKernelsize*lKernelsize*lKernelsize;
   double *cubeX = new double [lNCube];
   double *cubeY = new double [lNCube];
   double *cubeZ = new double [lNCube];
   double *cubeR = new double [lNCube]; // contains distance from center of cube to all voxels in the cube
   int *lIndicesOrder = new int [lNCube];

    // Initialise indices
    for ( int i = 0; i < lNCube; i++ )
        lIndicesOrder[i] = i;
      
   double *temp;
   temp = new double [lKernelsize];
   
   linspace(temp, -lKernelExtend, lKernelExtend, lKernelsize);
   /*
   for (int o=0; o<lKernelsize; o++)
	   mexPrintf("%f \n",temp[o]);
*/
   for (int  l = 0; l < lKernelsize; l++ )
        for (int  m = 0; m < lKernelsize; m++ )
            for (int n = 0; n < lKernelsize; n++ )
			{
                long Index = l + m*lKernelsize + n*lKernelsizesq;
                cubeX[Index] = temp[l];
                cubeY[Index] = temp[m];
                cubeZ[Index] = temp[n];
                cubeR[Index] = sqrt(cubeX[Index]*cubeX[Index] + cubeY[Index]*cubeY[Index] + cubeZ[Index]*cubeZ[Index]);
				// mexPrintf("Size = %d %d %d %f %f\n",lKernelExtend, lKernelsize, Index, temp[n], cubeR[Index]);
            }        
     
    // Sort by distance
    for (long nStartIndex = 0; nStartIndex < lNCube; nStartIndex++) // lNCube = lKernelsize*lKernelsize*lKernelsize;
    {
        // nSmallestIndex is the index of the smallest element we've encountered so far.
        long nSmallestIndex = nStartIndex;
        // Search through every element starting at nStartIndex+1
        for (long nCurrentIndex = nStartIndex + 1; nCurrentIndex < lNCube; nCurrentIndex++)
            // If the current element is smaller than our previously found smallest
            if (cubeR[nCurrentIndex] < cubeR[nSmallestIndex])
                // Store the index in nSmallestIndex
                nSmallestIndex = nCurrentIndex;
              
        // Swap our start element with our smallest element
        double tempA = cubeR[nStartIndex];
        cubeR[nStartIndex] = cubeR[nSmallestIndex];
        cubeR[nSmallestIndex] = tempA;
		// Save index of the our smallest element
        long tempB = lIndicesOrder[nStartIndex];
        lIndicesOrder[nStartIndex] = lIndicesOrder[nSmallestIndex];
        lIndicesOrder[nSmallestIndex] = tempB;
    }

   /*
   for ( int z = 0; z < lNCube; z++ ){
         mexPrintf("cubeR[%d] = %f  \n",z,cubeR[z]);  
   }
   */
   
    // Do pregridding by finding the pixel in the standard image closest to the given pixel in the original image
	// lNOrigIndices = NX*NY*NZ
	// lNCube = lKernelsize*lKernelsize*lKernelsize
	// size(lBinaryMask) = lNST*lNST*lNST
    // parallel_for (size_t(0),  (size_t)lNOrigIndices, [&](size_t i) {  
    for ( int i = 0; i < lNOrigIndices; i++ ){     
        // Find nearest pixel in new image    
        double Xnearest = X[i]/dPixelSize + lCentreSTR; // lCentreSTR = (201+1)/2 - 1 = 100;    
        double Ynearest = Y[i]/dPixelSize + lCentreSTR;
        double Znearest = Z[i]/dPixelSize + lCentreSTR;      
		//if (Xnearest>lNST || Ynearest>lNST || Znearest>lNST)
		//	mexPrintf("%d  %f %f %f\n", i, Xnearest, Ynearest, Znearest);
      
		// a pixel in the old image can support up to lNCube voxels in neighborhood in the new image. Here we label all voxels in the new image which are supported.
        for ( int cSortedIndex = 1; cSortedIndex < lNCube; cSortedIndex++ ){
            long cubeIndex = lIndicesOrder[cSortedIndex];
            long cIndexSTImageX = Xnearest+cubeX[cubeIndex];
            long cIndexSTImageY = Ynearest+cubeY[cubeIndex];
            long cIndexSTImageZ = Znearest+cubeZ[cubeIndex];
		
            if( cIndexSTImageX>=0 && cIndexSTImageY>=0  && cIndexSTImageZ>=0 && cIndexSTImageX<lNST && cIndexSTImageY<lNST && cIndexSTImageZ<lNST ){		
				long cIndexSTImage = cIndexSTImageX + cIndexSTImageY*lNST + cIndexSTImageZ*lNSTsq;			
                lBinaryMask[cIndexSTImage] = 1;     
			}				
        }
    } 
     //);

   // lNOrigIndices = NX*NY*NZ;   
   // parallel_for (size_t(0),  (size_t)lNOrigIndices, [&](size_t i) {     
    for ( int i = 0; i < lNOrigIndices; i++ ){
        // Find nearest pixel in new image    
        double Xnearest = X[i]/dPixelSize + lCentreSTR;        
        double Ynearest = Y[i]/dPixelSize + lCentreSTR;
        double Znearest = Z[i]/dPixelSize + lCentreSTR;

        long cubeIndex = lIndicesOrder[0];

        long cIndexSTImageX = Xnearest+cubeX[cubeIndex];
        long cIndexSTImageY = Ynearest+cubeY[cubeIndex];
        long cIndexSTImageZ = Znearest+cubeZ[cubeIndex];
		
		// the value of the closest points between two images is copied from old image to the new image, then the binary mask is filled by 2
		// sometimes the closest point is not uniqe! 
        if( cIndexSTImageX>=0 && cIndexSTImageY>=0  && cIndexSTImageZ>=0 && cIndexSTImageX<lNST && cIndexSTImageY<lNST && cIndexSTImageZ<lNST ){
            long cIndexSTImage = cIndexSTImageX + cIndexSTImageY*lNST + cIndexSTImageZ*lNSTsq;					
            lBinaryMask[cIndexSTImage]  = 2;            
            dImageSTtemp[cIndexSTImage] = dImage[i];
			dImageST[cIndexSTImage]     = dImage[i];  				  
        }
    }
	// mexPrintf("%f   %f\n",lCentreSTR, dPixelSize);
   //);

   
   // Loop over all pixel of standard image
   //parallel_for (size_t(0),  (size_t)lNST, [&](size_t x)   {
    for ( int x = 0; x < lNST; x++ ){   
        for ( int y = 0; y < lNST; y++ ){
            for ( int z = 0; z < lNST; z++ ){    
			   long lMainIndex = x + y*lNST + z*lNSTsq;       
			   if (lBinaryMask[ lMainIndex ] == 1){ // is in neighborhood of one of the old image voxels		   
					double dWeightTot = 0.0;                  
					int SamplesFoundinVicinity = 0;		
					for ( int cSortedIndex = 0; cSortedIndex < lNCube; cSortedIndex++ ){					
						if(SamplesFoundinVicinity == 0){
							long cubeIndex = lIndicesOrder[cSortedIndex];
							double cIndexSTImageX = x+cubeX[cubeIndex];
							double cIndexSTImageY = y+cubeY[cubeIndex];
							double cIndexSTImageZ = z+cubeZ[cubeIndex];

							if( cIndexSTImageX>=0 && cIndexSTImageY>=0  && cIndexSTImageZ>=0 && cIndexSTImageX<lNST && cIndexSTImageY<lNST && cIndexSTImageZ<lNST ){
								long cIndexSTImage = cIndexSTImageX + cIndexSTImageY*lNST + cIndexSTImageZ*lNSTsq;
								if (lBinaryMask[ cIndexSTImage ] == 2){																		
									double dDist   = cubeR[cSortedIndex]/cubeR[lNCube-1];
									double dWeight =  pow(cos(dDist*pi/2.0),2);//1.0-dDist;//						
									dImageST[lMainIndex] += dImageSTtemp[cIndexSTImage]*dWeight;																								
									dWeightTot += dWeight;
									SamplesFoundinVicinity++;
								}
							}
						}
					}
					if(dWeightTot != 0.0)
						dImageST[lMainIndex] = dImageST[lMainIndex]/dWeightTot;				 
                }
            }
        }
	}
	//);
    
   delete []temp;

   delete []cubeX;
   delete []cubeY;
   delete []cubeZ;
   delete []cubeR;
   
   delete []lIndicesOrder;
   delete []lBinaryMask;
   delete []dImageSTtemp;
}

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    mexEvalString("drawnow");
    ///////////////////////////////////////////////////////////////////////
    // Variables
    ///////////////////////////////////////////////////////////////////////
    double *dtempPointer;
    double *dImage;
    
    ///////////////////////////////////////////////////////////////////////
    // Check input
    ///////////////////////////////////////////////////////////////////////
    
    if (nrhs<7){
        mexErrMsgTxt("At least 7 arguments have to be passed!\n");
    }
    
    ///////////////////////////////////////////////////////////////////////
    // Get Image
    ///////////////////////////////////////////////////////////////////////
    if (!mxIsDouble(prhs[0]))mexErrMsgTxt("First argument must be of double precision\n");
    
    dImage = mxGetPr(prhs[0]);
         
    const mwSize *dims_data = mxGetDimensions(prhs[0]);
    mwSize number_of_dimensions_data = mxGetNumberOfDimensions(prhs[0]);
  
    int NX = dims_data[0], NY = 1, NZ = 1;
	if (number_of_dimensions_data>1)
		NY = dims_data[1];
	if (number_of_dimensions_data>2)
        NZ = dims_data[2];
	if (number_of_dimensions_data>3)
        mexErrMsgTxt("Only up to 3D data supported \n"); 
               
    ///////////////////////////////////////////////////////////////////////
    // Get coordinates
    ///////////////////////////////////////////////////////////////////////
    // X	
    if (!mxIsDouble(prhs[1]))mexErrMsgTxt("Xcoords argument must be of double precision\n");
    if (mxIsComplex(prhs[1]))mexErrMsgTxt("Xcoords argument must be real\n"); 
    double *dXCoords = mxGetPr(prhs[1]);
    
	int i = 0;
    const mwSize *dims = mxGetDimensions(prhs[1]);
    mwSize number_of_dimensions = mxGetNumberOfDimensions(prhs[1]); 
	for(i=0; i<number_of_dimensions; i++)
		if (dims_data[i] != dims[i] || number_of_dimensions != number_of_dimensions_data)
			mexErrMsgTxt("Xcoords: Dimension mismatch\n"); 
    
    // Y
    if (!mxIsDouble(prhs[2]))mexErrMsgTxt("Ycoords argument must be of double precision\n");
    if (mxIsComplex(prhs[2]))mexErrMsgTxt("Ycoords argument must be real\n"); 
    double *dYCoords = mxGetPr(prhs[2]);
    
    dims = mxGetDimensions(prhs[2]);
    number_of_dimensions = mxGetNumberOfDimensions(prhs[2]);
	for(i=0; i<number_of_dimensions; i++)
		if (dims_data[i] != dims[i] || number_of_dimensions != number_of_dimensions_data)
			mexErrMsgTxt("Ycoords: Dimension mismatch\n"); 
      
    // Z
    if (!mxIsDouble(prhs[3]))mexErrMsgTxt("Zcoords argument must be of double precision\n");
    if (mxIsComplex(prhs[3]))mexErrMsgTxt("Zcorrds argument must be real\n"); 
    double *dZCoords = mxGetPr(prhs[3]);
    
    dims = mxGetDimensions(prhs[3]);
    number_of_dimensions = mxGetNumberOfDimensions(prhs[3]);
	for(i=0; i<number_of_dimensions; i++)
		if (dims_data[i] != dims[i] || number_of_dimensions != number_of_dimensions_data)
			mexErrMsgTxt("Zcorrds: Dimension mismatch\n"); 
    
    
    ///////////////////////////////////////////////////////////////////////
    // Get params
    ///////////////////////////////////////////////////////////////////////
    if (!mxIsDouble(prhs[4]))mexErrMsgTxt("5th argument must be of double precision\n");
    if (mxIsComplex(prhs[4]))mexErrMsgTxt("5th argument must be real\n"); 
    
    if (!mxIsDouble(prhs[5]))mexErrMsgTxt("6th argument must be of double precision\n");
    if (mxIsComplex(prhs[5]))mexErrMsgTxt("6th argument must be real\n"); 
    
    if (!mxIsDouble(prhs[6]))mexErrMsgTxt("7th argument must be of double precision\n");
    if (mxIsComplex(prhs[6]))mexErrMsgTxt("7th argument must be real\n"); 
    
    dtempPointer = mxGetPr(prhs[4]);
    double dFOVST        = dtempPointer[0];
    //mexPrintf("FOVST: %f  \n",dFOVST);
    
    dtempPointer = mxGetPr(prhs[5]);
    long lNST       = (long)dtempPointer[0];
    //mexPrintf("NST: %d  \n",lNST);
        
    dtempPointer = mxGetPr(prhs[6]);
    double dPixelSizeOrig   = dtempPointer[0];
    //mexPrintf("Kernelsize: %d  \n",lKernelsize);
            
    // ######################################################################
	// Output variables
	// ######################################################################
    mwSize ndim = 3;
    mwSize *Dims  = new mwSize[ndim];
    
    Dims[0] = lNST;
    Dims[1] = lNST;
    Dims[2] = lNST;
    
    double *dImageST;
    plhs[0] = mxCreateNumericArray(ndim,Dims , mxDOUBLE_CLASS, mxREAL);
    dImageST   = (double*)mxGetPr(plhs[0]);

    InterpolateToSTRFoV(dImage,dXCoords,dYCoords,dZCoords,dFOVST,lNST,dPixelSizeOrig, NX, NY, NZ, dImageST);
    delete[] Dims;  
}