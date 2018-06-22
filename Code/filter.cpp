#include "filter.h"
#include "stdio.h"
#include "math.h"

void filter::resetFilter()
{
    // Clear the delayed input values
    xLast = 0;
    xSecondLast = 0;
    
    // Clear the delayed output values
    yLast = 0;
    ySecondLast = 0;
}

// Method to invoke the caluclation of the filter coefficients.
void filter::getCoefficients(double cutoff, double fs)
{
    // Calculate the Filters
    calculateLPF(cutoff, fs);
    calculateHPF(cutoff, fs);
}

// Method to calculate the coefficients of the Low Pass Filter
void filter::calculateLPF(double cutoff, double fs)
{	
    const double wa = tan((M_PI * cutoff)/fs);	// analog cutoff frequency
    const double c = 1/wa;
    
    lpfA[0] = 1 / (1 + (sqrt(2)*c) + (c*c));
    lpfA[1] = 2 * lpfA[0];
    lpfA[2] = lpfA[0];

    lpfB[1] = 2 * lpfA[0] * (1 - (c*c));
    lpfB[2] = lpfA[0] * (1 - (sqrt(2)*c) + (c*c));
}

// Method to calculate the coefficients of the High Pass Filter
void filter::calculateHPF(double cutoff, double fs)
{
    const double wa = tan((M_PI * cutoff)/fs);	// analog cutoff frequency
    const double c = wa;
    
    hpfA[0] = 1 / (1 + (sqrt(2)*c) + (c*c));
    hpfA[1] = -2 * hpfA[0];
    hpfA[2] = hpfA[0];
    
    hpfB[1] = 2 * hpfA[0] * ((c*c) - 1);
    hpfB[2] = hpfA[0] * (1 - (sqrt(2)*c) + (c*c));
}

// Method to set the filter coefficients calculated earlier
void filter::setCoefficients(const std::vector<double>& bUpdated, const std::vector<double>& aUpdated)
{
    a[0] = aUpdated[0];
    a[1] = aUpdated[1];
    a[2] = aUpdated[2];
    
    b[1] = bUpdated[1];
    b[2] = bUpdated[2];
}

// Method to apply the filter to the signal
float filter::processFilter(float y)
{
    const double x = y;
    
    y = a[0]*x + a[1]*xLast + a[2]*xSecondLast - b[1]*yLast - b[2]*ySecondLast;
    
    xSecondLast = xLast;
    xLast = x;
    ySecondLast = yLast;
    yLast = y;
    
    return y;
}


