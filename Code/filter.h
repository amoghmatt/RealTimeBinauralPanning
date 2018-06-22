#ifndef filter_h
#define filter_h

#include <vector>

class filter
{
public:
    filter() {};
    ~filter() {};
    
    void resetFilter();
    void setCoefficients(const std::vector<double>& b, const std::vector<double>& a);
    float processFilter(float x);
    
    void getCoefficients(double cutoff, double fs);

    const std::vector<double>& getLPFA() {return lpfA; }
    const std::vector<double>& getLPFB() {return lpfB; }

    const std::vector<double>& getHPFA() {return hpfA; }
    const std::vector<double>& getHPFB() {return hpfB; }
    
private:
    // Numerator and Denominator Coefficients for the Low Pass Filter
    std::vector<double> lpfA {1.0, 0.0, 0.0};
    std::vector<double> lpfB {1.0, 0.0, 0.0};

    // Numerator and Denominator Coefficients for the Low Pass Filter
    std::vector<double> hpfA {1.0, 0.0, 0.0};
    std::vector<double> hpfB {1.0, 0.0, 0.0};

    void calculateLPF(double cutoff, double fs);
    void calculateHPF(double cutoff, double fs);
	
    std::vector<double> a = {1.0, 0.0, 0.0}; // Setting Coefficients for Low Pass as Default
    std::vector<double> b = {1.0, 0.0, 0.0};
    
    double xLast {0.0};
    double xSecondLast {0.0};
    double yLast {0.0};
    double ySecondLast {0.0};

};

#endif /* filter_h */
