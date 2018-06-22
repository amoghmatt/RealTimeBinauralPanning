#include <Bela.h>
#include <cmath>
#include "filter.h"

//#define DELAY_BUFFER_SIZE 1100
#define LPF_FS 4000
#define BANDPASSLOWER_FS 4000
#define BANDPASSHIGHER_FS 9000
#define DELAY_BUFFER_SIZE 2205

filter gLPF;
filter gBandpassHPF;
filter gBandpassLPF;

// Enumeration to set the orientation states
enum {
	flat = 0,
    upsideDown,
};

//This variables indicates whether the effect should be triggered or not
int gEffectApply =  0;

// Varilables for button 1
int gButtonPin= P8_07;
int gButtonState = 0;
int gButtonPreviousState = 0;

// Varilables for LED
int gLedPin = P8_09;
int gLedState = 0;

// Varilables for Accelerometer
int gAccelPinX = 1;
int gAccelPinY = 2;
int gAccelPinZ = 3;
int gNumAudioFramesPerAnalog;

// Varilables for Orientation
int gOrientation;
float gXState = 0;
float gYState = 0;
float gZState = 0;

// Varilables for Hysterisis while calculating orientation
float gHysterisisUpperThreshold = 0.07;
float gHysterisisLowerThreshold = 0.04;

// Variables for Spatial Audio Calculation
float gTheta = (M_PI)/4;	// Angle of the speakers
float gPhi = 0; 	// Angle of the source
float gLeftGain;
float gRightGain;

// Buffer holding previous samples per channel
float gBackDelayBuffer_l[DELAY_BUFFER_SIZE] = {0};
float gBackDelayBuffer_r[DELAY_BUFFER_SIZE] = {0};

float gFrontDelayBuffer_l[DELAY_BUFFER_SIZE] = {0};
float gFrontDelayBuffer_r[DELAY_BUFFER_SIZE] = {0};

float gDelayBuffer_l[DELAY_BUFFER_SIZE] = {0};
float gDelayBuffer_r[DELAY_BUFFER_SIZE] = {0};

// Write pointer
int gDelayBufWritePtr = 0;
// Amount of delay
float gDelayAmount = 0.75;
// Amount of feedback
float gDelayFeedbackAmount = 0.25;
// Level of pre-delay input
float gDelayAmountPre = 0.75;
// Amount of delay in samples (needs to be smaller than or equal to the buffer size defined above)
int gDelayInSamples = DELAY_BUFFER_SIZE;

// Variables for Head Related Transfer Functions
float gBeta = 0;
float gAlpha = 0;
float gSamplePeriod;

float gYLeft = 0;
float gYPreviousLeft = 0;
float gYRight = 0;
float gYPreviousRight = 0;

float gXLeft = 0;
float gXPreviousLeft = 0;
float gXRight = 0;
float gXPreviousRight = 0;

// Method to calculate the orientation of the Accelerometer and set the state
void calculateOrientation(float x, float y, float z)
{
    if(gXState == 0 && fabs(x) > gHysterisisUpperThreshold)
        gXState = x/fabs(x);
    else if (gXState != 0 && fabs(x) < gHysterisisLowerThreshold)
        gXState = 0;
    
    if(gYState == 0 && fabs(y) > gHysterisisUpperThreshold)
        gYState = y/fabs(y);
    else if (gYState != 0 && fabs(y) < gHysterisisLowerThreshold)
        gYState = 0;
    
    if(gZState == 0 && fabs(z) > gHysterisisUpperThreshold)
        gZState = z/fabs(z);
    else if (gZState != 0 && fabs(z) < gHysterisisLowerThreshold)
        gZState = 0;
    
    if (gXState == 0 && gYState == 0 && gZState == 1){
        gOrientation = flat;
    }
    if (gXState == 0 && gYState == 0 && gZState == -1){
        gOrientation = upsideDown;
    }
}

bool setup(BelaContext *context, void *userData)
{	
    pinMode(context, 0, gButtonPin, INPUT);
    pinMode(context, 0, gLedPin, OUTPUT);

	// For this example we need the same amount of audio and analog input and output channels
	if(context->audioInChannels != context->audioOutChannels ||
			context->analogInChannels != context-> analogOutChannels){
		printf("Error: for this project, you need the same number of input and output channels.\n");
		return false;
	}
	// Compute the coefficients based ont the frequency and the sampling rate
	gLPF.getCoefficients(LPF_FS, context->audioSampleRate);
	gBandpassLPF.getCoefficients(BANDPASSLOWER_FS, context->audioSampleRate);
	gBandpassHPF.getCoefficients(BANDPASSHIGHER_FS, context->audioSampleRate);
	
    // Clear any previous values	
	gLPF.resetFilter();
	gBandpassLPF.resetFilter();
	gBandpassHPF.resetFilter();
	
    // Set the calculated coefficient values.
	gLPF.setCoefficients(gLPF.getLPFB(), gLPF.getLPFA());
	gBandpassLPF.setCoefficients(gBandpassLPF.getLPFB(), gBandpassLPF.getLPFA());
	gBandpassHPF.setCoefficients(gBandpassHPF.getLPFB(), gBandpassHPF.getLPFA());

	gNumAudioFramesPerAnalog = context->audioFrames / context->analogFrames;
	gSamplePeriod = 1 / context->audioSampleRate;
	// Beta = 2 * speed of sound / head radius
	gBeta = (2 * 343) / 0.9;
	
	return true;
}

void render(BelaContext *context, void *userData)
{
	for(unsigned int n = 0; n < context->audioFrames; n++) {
	
		// Obtain the incoming audio
		float leftInput = audioRead(context, n, 0);
		float rightInput = audioRead(context, n, 1);
		
		// Declare the output audio variables and make it equal to the input
        float output_l = leftInput;
        float output_r = rightInput;
	
		// Read the button state
		gButtonState = digitalRead(context, n, gButtonPin);
		
        // Obtain the x,y,z axis values from the Accelerometer and calibrate it
        float x = analogRead(context, n/gNumAudioFramesPerAnalog, gAccelPinX) - 0.415;
        float y = analogRead(context, n/gNumAudioFramesPerAnalog, gAccelPinY) - 0.436;
        float z = analogRead(context, n/gNumAudioFramesPerAnalog, gAccelPinZ) - 0.415;
        
        // Obtain the orientation of the Accelerometer
        calculateOrientation(x,y,z);
		
		// Map the obtained Accelerometer values
		gPhi = map(x, -0.201, 0.201, gTheta, -gTheta);
		float lrDepth = map((0.201 - fabs(x)), 0, 0.201, 0, 1);
        float fbDepth = map((0.202 - fabs(y)), 0, 0.202, 0, 1);
        float goingBack = map(y, -0.202, 0.202, 1, 0);
        float goingFront = 1 - goingBack;

		// Calculate the Left and Right Gain using ILD
		gLeftGain = (cos(gPhi) * sin(gTheta) + sin(gPhi) * cos(gTheta)) / (2 * cos(gTheta) * sin(gTheta));
		gRightGain = (cos(gPhi) * sin(gTheta) - sin(gPhi) * cos(gTheta)) / (2 * cos(gTheta) * sin(gTheta));
	
		// Head Model Filter
		gAlpha = sin(1 - goingBack);
		
		gXLeft = leftInput;
		gXRight = rightInput;
		
		// Filter coefficient calculations
		float b1 = (2 * gAlpha) + (gBeta * gSamplePeriod);
		float b2 = (gBeta * gSamplePeriod) - (2 * gAlpha);
		float a1 = (2 + (gBeta * gSamplePeriod));
		float a2 = ((gBeta * gSamplePeriod) - 2);
		
		// Discrete Time Domain Transfer Function
		gYLeft = ((b1 * gXLeft + b2 * gXPreviousLeft - a2 * gYPreviousLeft) / a1);
		gYRight = ((b1 * gXRight + b2 * gXPreviousRight - a2 * gYPreviousRight) / a1);
		 
		float gHeadFilteredLeftOutput = gYLeft;
		float gHeadFilteredRightOutput = gYRight;
		
		// Update previous values
		gYPreviousLeft = gYLeft;
		gYPreviousRight = gYRight;
		gXPreviousLeft = gXLeft;
		gXPreviousRight = gXRight;
		
		// Low Pass Filtering, f < 4kHz
		float lpfLeft = gLPF.processFilter(leftInput);
		float lpfRight = gLPF.processFilter(rightInput);
		
		// High Pass Filtering, f > 4kHz 
       	float hpfLeft = gBandpassHPF.processFilter(leftInput);
       	float hpfRight = gBandpassHPF.processFilter(rightInput);
       	
       	// Band Pass Filtering, 4kHz < f > 7kHz
       	float bandpassedLeft = gBandpassLPF.processFilter(hpfLeft);
       	float bandpassedRight = gBandpassLPF.processFilter(hpfRight);
       	
       	// Delay buffer with the bandpass filtered input for when the audio transitions front
       	float delayedFrontLeftInput =  (gDelayAmountPre * bandpassedLeft + gFrontDelayBuffer_l[(gDelayBufWritePtr - gDelayInSamples + DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE] * gDelayFeedbackAmount);
       	float delayedFrontRightInput = (gDelayAmountPre * bandpassedRight + gFrontDelayBuffer_r[(gDelayBufWritePtr - gDelayInSamples + DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE] * gDelayFeedbackAmount);
        
        // Delay buffer with the lowpass filtered input for when the audio transitions back
		float delayedBackLeftInput = (gDelayAmountPre * lpfLeft + gBackDelayBuffer_l[(gDelayBufWritePtr - gDelayInSamples + DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE] * gDelayFeedbackAmount);
		float delayedBackRightInput = (gDelayAmountPre * lpfRight + gBackDelayBuffer_r[(gDelayBufWritePtr - gDelayInSamples + DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE] * gDelayFeedbackAmount);
		
		// Delay buffer with the lowpass filtered input for when the audio transitions to the sides
   		gDelayBuffer_l[gDelayBufWritePtr] = (gDelayAmountPre * lpfLeft  + gDelayBuffer_l[(gDelayBufWritePtr - gDelayInSamples + DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE] * gDelayFeedbackAmount);
		gDelayBuffer_r[gDelayBufWritePtr] = (gDelayAmountPre * lpfRight + gDelayBuffer_r[(gDelayBufWritePtr - gDelayInSamples + DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE] * gDelayFeedbackAmount);
		
        gFrontDelayBuffer_l[gDelayBufWritePtr] = delayedFrontLeftInput;
        gFrontDelayBuffer_r[gDelayBufWritePtr] = delayedFrontRightInput;
		
        gBackDelayBuffer_l[gDelayBufWritePtr] = delayedBackLeftInput;
        gBackDelayBuffer_r[gDelayBufWritePtr] = delayedBackRightInput;
        
		// 2 Pinna Reflection
		float pinnaLeft = 0;
		float pinnaRight = 0;
		
		// tau_k (theta,phi) = A_k cos(phi/2)sin(D_k(pi/2 - theta)) + B_k
		// phi = elevaltion = 30 degrees
		float tau = 1 * cos((gPhi)/2) * sin(1 * (M_PI/2 - (gTheta))) + 2;
		
		// s_out(n) = s_in(n) + sum(for all k (p_k * s_in(n - tau)))
		float leftSum = 0;
		float rightSum = 0;
		for(int i = 0; i < 2; i++) {
			leftSum +=  0.5 * gBackDelayBuffer_l[gDelayBufWritePtr - int(ceil(tau) + DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE];
			rightSum += 0.5 * gBackDelayBuffer_r[gDelayBufWritePtr - int(ceil(tau) + DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE];
		}
		pinnaLeft = gHeadFilteredLeftOutput + leftSum;
		pinnaRight = gHeadFilteredRightOutput + rightSum;
			
        // Read the button and toggle effect state
        if(gButtonState == 1 && gButtonPreviousState == 0) {
        	gEffectApply = !gEffectApply;
        }
    	gButtonPreviousState = gButtonState;
    	
    	// Configuring the various panning states
    	// Centre of the space
    	float midLeft = (fbDepth * lrDepth * leftInput);
    	float midRight = (fbDepth * lrDepth * rightInput);
    	
    	// Transitioning to the front
    	float movingFrontLeft = 0.5 * goingFront * (leftInput + gFrontDelayBuffer_l[(gDelayBufWritePtr-gDelayInSamples+DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE] * gDelayAmount);
    	float movingFrontRight = 0.5 * goingFront * (rightInput + gFrontDelayBuffer_r[(gDelayBufWritePtr-gDelayInSamples+DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE] * gDelayAmount);
    	
    	// Transitioning to the back
    	float movingBackLeft = (goingBack) * (1 - fbDepth) * (pinnaLeft + gBackDelayBuffer_l[(gDelayBufWritePtr-gDelayInSamples+DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE] * gDelayAmount);
    	float movingBackRight = (goingBack) * (1 - fbDepth) * (pinnaRight + gBackDelayBuffer_r[(gDelayBufWritePtr-gDelayInSamples+DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE] * gDelayAmount);
    	
    	// Transitioning to the left / right 
    	float leftDepth = (fbDepth * (1 - lrDepth) * (leftInput + gDelayBuffer_l[(gDelayBufWritePtr-gDelayInSamples+DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE] * gDelayAmount));
    	float rightDepth = (fbDepth * (1 - lrDepth) * (rightInput + gDelayBuffer_r[(gDelayBufWritePtr-gDelayInSamples+DELAY_BUFFER_SIZE)%DELAY_BUFFER_SIZE] * gDelayAmount));
    	
    	// Transitioning up / down
    	float leftDown = (fbDepth * lrDepth * lpfLeft);
    	float rightDown = (fbDepth * lrDepth * lpfRight);
    	
    	// Alter the output if the effect has to be applied (button pressed)
    	if(gEffectApply){
    		if(gOrientation == upsideDown) {
    			output_l = gLeftGain * (leftDown + leftDepth + movingBackLeft + movingFrontLeft);
    			output_r = gRightGain * (rightDown + rightDepth + movingBackRight + movingFrontRight);
    		}
    		else{
	    		output_l = gLeftGain * (midLeft + leftDepth + movingBackLeft + movingFrontLeft);
    			output_r = gRightGain * (midRight + rightDepth + movingBackRight + movingFrontRight);
    		}
    	}
    	// Turn on/off LED to indicate effect applied status
    	digitalWrite(context, n, gLedPin, gEffectApply);
    	
        // Write the sample into the output buffer -- done!
        audioWrite(context, n, 0, output_l);
        audioWrite(context, n, 1, output_r);
        
        // Increment delay buffer write pointer
        if(++gDelayBufWritePtr>DELAY_BUFFER_SIZE)
            gDelayBufWritePtr = 0;
	}
}

void cleanup(BelaContext *context, void *userData)
{

}