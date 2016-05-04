/**
	@file
	matffttrig~: a simple audio object for Max
	original by: matthew aidekman
	@ingroup examples
*/

#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "z_dsp.h"			// required for MSP objects
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "fft.h"

// struct to represent the object's state
typedef struct _matffttrig {
	t_pxobject		ob;			// the object itself (t_pxobject in MSP instead of t_object)
    void            *out;
    void            *out2;
    void            *out3;
    void            *out4;
    
    int             barkBins[60];
    int             barkFreqsLength;
    int             barkCenters[60];
    
    double          rollingBuffer[2048];  // keep track of last 2048 samples
    long            rollingBufCounter; // where are we in rolling Buffer
    long            paddingSamples;  //analyzes this many samples prior to trigger.  (theoretically allowing quicker times between pulse and analysis.
    int             fftLength;  //512
    long            fftingModeCounter;  //count down from when click is heard to when to perform analysis.
    
} t_matffttrig;


// method prototypes
void *matffttrig_new(t_symbol *s, long argc, t_atom *argv);
void matffttrig_free(t_matffttrig *x);
void matffttrig_assist(t_matffttrig *x, void *b, long m, long a, char *s);
void matffttrig_float(t_matffttrig *x, double f);
void matffttrig_dsp64(t_matffttrig *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void matffttrig_perform64(t_matffttrig *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);


// global class pointer variable
static t_class *matffttrig_class = NULL;


//***********************************************************************************************

void ext_main(void *r)
{
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	t_class *c = class_new("mat.ffttrig~", (method)matffttrig_new, (method)dsp_free, (long)sizeof(t_matffttrig), 0L, A_GIMME, 0);

	class_addmethod(c, (method)matffttrig_float,		"float",	A_FLOAT, 0);
	class_addmethod(c, (method)matffttrig_dsp64,		"dsp64",	A_CANT, 0);
	class_addmethod(c, (method)matffttrig_assist,       "assist",	A_CANT, 0);

	class_dspinit(c);
	class_register(CLASS_BOX, c);
	matffttrig_class = c;
}


void *matffttrig_new(t_symbol *s, long argc, t_atom *argv)
{
	t_matffttrig *x = (t_matffttrig *)object_alloc(matffttrig_class);

	if (x) {
		dsp_setup((t_pxobject *)x, 2);	// MSP inlets: arg is # of inlets and is REQUIRED!
		// use 0 if you don't need inlets
        x->rollingBufCounter = 0;
        x->paddingSamples = 30;
        x->fftLength=512;
        x->fftingModeCounter = -1;  //dont count right now
		x->out=outlet_new(x, NULL); 		// signal outlet (note "signal" rather than NULL)
        x->out2=outlet_new(x, NULL); 		// signal outlet (note "signal" rather than NULL)
        x->out3=outlet_new(x, NULL); 		// signal outlet (note "signal" rather than NULL)
        x->out4=outlet_new(x, NULL); 		// signal outlet (note "signal" rather than NULL)
	}
	return (x);
}


// NOT CALLED!, we use dsp_free for a generic free function
void matffttrig_free(t_matffttrig *x)
{
	;
}


void matffttrig_assist(t_matffttrig *x, void *b, long m, long a, char *s)
{
	if (m == ASSIST_INLET) { //inlet
		sprintf(s, "I am inlet %ld", a);
	}
	else {	// outlet
		sprintf(s, "I am outlet %ld", a);
	}
}


void matffttrig_float(t_matffttrig *x, double f)
{
}


// registers a function for the signal chain in Max
void matffttrig_dsp64(t_matffttrig *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
	post("my sample rate is: %f", samplerate);

    
    //These are the bark edge frequencies
    float barkFreqs[] = {100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500};
    int lastBinVal = 0;
    int lastBinIdx =-1;
    
    //these would be the bin indexes for each bark frequency, exchept we're going to fill in a few more.
    int tempBarkBins[24];
 
    // calculate bin indexes, leave out repeats;
    for(int i = 0 ; i < 24;i++){
        int calcThisBin=ceil(barkFreqs[i]*x->fftLength/samplerate);
        if(lastBinVal!=calcThisBin){
            lastBinVal = calcThisBin;
            lastBinIdx = lastBinIdx+1;
            tempBarkBins[lastBinIdx] = lastBinVal;
        }
        
    }

    int onesUntil =5;  // we're going to ABSOLUTELY have a bark band every bin until bin number onesUntil;
    int indexBeforeOnesUntil = 0;
    for(int i = 0 ; i < 24;i++){
        if(tempBarkBins[i]<onesUntil){
            indexBeforeOnesUntil = i;
        }
        else{
            i = 25;
        }
    }
    
    int addingIndex = 0;
    for(int i = 1 ; i < onesUntil;i++){
        x->barkBins[addingIndex] = i;
        addingIndex++;
    }
    
    for(int i = indexBeforeOnesUntil+1;i<24;i++){
        x->barkBins[addingIndex] = tempBarkBins[i];
        addingIndex++;
    }
    
    x->barkFreqsLength = addingIndex;

    
    //calculate the center of each band for centroid analysis.
    x->barkCenters[0] = 0;
    for(int i = 1 ; i < x->barkFreqsLength;i++){
        x->barkCenters[i] = x->barkBins[i]*samplerate/x->fftLength;  // temporarily fill with frequency of each band edge.
    }
    
    for(int i = 0 ; i < x->barkFreqsLength;i++){
        x->barkCenters[i] = (x->barkCenters[i]+x->barkCenters[i+1])/2;  //average out band edges to get band center
    }
    
    // the last one is between itself and the nyquist frequency
    x->barkCenters[x->barkFreqsLength] = (x->barkCenters[x->barkFreqsLength-1]+samplerate/2)/2;
    
	// instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
	// the arguments passed are:
	// 1: the dsp64 object passed-in by the calling function
	// 2: the symbol of the "dsp_add64" message we are sending
	// 3: a pointer to your object
	// 4: a pointer to your 64-bit perform method
	// 5: flags to alter how the signal chain handles your object -- just pass 0
	// 6: a generic pointer that you can use to pass any additional data to your perform method

	object_method(dsp64, gensym("dsp_add64"), x, matffttrig_perform64, 0, NULL);
}


// this is the 64-bit perform method audio vectors
void matffttrig_perform64(t_matffttrig *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
//	t_double *inL = ins[0];		// we get audio for each inlet of the object from the **ins argument
//	t_double *outL = outs[0];	// we get audio for each outlet of the object from the **outs argument
	int n = sampleframes;

	// this perform method simply copies the input to the output, offsetting the value

    if(numins!=2){
        return;
    }
    
    for(int i =  0 ; i < n ; i++){
        
        
        //// when we get a click in the right inlet
        // we set a count down till we can start FFTing again
        if(ins[1][i]>.5){
            x->fftingModeCounter = x->fftLength - x->paddingSamples;
        }
        
        x->rollingBuffer[x->rollingBufCounter] = ins[0][i];
        
        //increment and wrap rolling buffer counter
        x->rollingBufCounter = (x->rollingBufCounter+ 1) % 2048;
        
        if(x->fftingModeCounter>=0){
            x->fftingModeCounter--;
        }
        
        
        
        
        //Perform all the analysis
        if(x->fftingModeCounter==0){
        
            //FFT IN rollingBufCounter to (fftLength-1), wraps around to 0 to (rolingBufCounter-1)
            double fftinputR[x->fftLength];
            double fftinputI[x->fftLength];
            int fftStartIndex = (2048+x->rollingBufCounter - x->fftLength)%2048;
            int readIndex;
            
            
            //fill the arrays with values to FFT
            for(int j = 0 ; j < x->fftLength;j++){
                readIndex = (fftStartIndex+j)%2048;
                
                
                //terrible windowing
                double window = 1;
                if(j<10)
                    window = (double)j/10.0;
                else if(j>x->fftLength-10)
                    window = abs((double)j-(double)x->fftLength)/10.0;
                
                fftinputR[j] = x->rollingBuffer[readIndex];
                fftinputI[j] = 0;
            }
            
            
            //perform actual FFT
            transform_radix2(fftinputR,fftinputI,x->fftLength);
            
            
            //Send out FFT
            t_atom myAtom[x->fftLength];
            for(int j = 0; j < x->fftLength;j++){
                atom_setfloat(&myAtom[j],fftinputR[j]);
            }
            outlet_list(x->out,NULL,x->fftLength,&myAtom);

            t_atom myAtom2[x->fftLength];
            for(int j = 0; j < x->fftLength;j++){
                atom_setfloat(&myAtom2[j],fftinputI[j]);
            }
            outlet_list(x->out2,NULL,x->fftLength,&myAtom2);
            
            
            
            //do the pseudobark analysis
            double runningBandMagnitude = 0.0;
            int    numRunningBins = 0;
            int    nextBandIndex = 0;
            double loopMagnitude;
            double barkAnalysis[x->barkFreqsLength+2];
            
            for(int j = 0 ; j < x->fftLength/2;j++){
                if(j==x->barkBins[nextBandIndex]){
                    barkAnalysis[nextBandIndex] = runningBandMagnitude/(double)numRunningBins;
                    nextBandIndex++;
                }
                
                loopMagnitude = sqrt(fftinputI[j]*fftinputI[j]+fftinputR[j]*fftinputR[j]);
                runningBandMagnitude+=loopMagnitude;
                numRunningBins++;
            }
            barkAnalysis[nextBandIndex] = runningBandMagnitude/(double)numRunningBins;
            
            t_atom myAtom3[x->barkFreqsLength+2];
            for(int j = 0; j < x->barkFreqsLength+2;j++){
                atom_setfloat(&myAtom3[j],barkAnalysis[j]);
            }
            outlet_list(x->out3,NULL, x->barkFreqsLength+2,&myAtom3);
            
            //perform centroid analysis on bark spectrum
            double weightedFreqSum = 0;
            double magsum = 0;
            for(int j = 0; j < x->barkFreqsLength+2;j++){
                weightedFreqSum+=x->barkCenters[j]*barkAnalysis[j];
                magsum+=barkAnalysis[j];
            }
            
            outlet_float(x->out4, weightedFreqSum/magsum);
            
            return;
            
        }
    }
}

