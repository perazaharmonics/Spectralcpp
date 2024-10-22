#include "IIR_Filter.h"

IIRFilter::IIRFilter(const std::vector<StkFloat>& numerator, const std::vector<StkFloat>& denominator)
    : filter_(numerator, denominator),
    output_(20, 1)   // These can be adjsuted
{
    output_[0] = 1.0;
}

void IIRFilter::applyFilter(StkFrames& frames) {
    filter_.tick(frames);
}

void IIRFilter::printOutput(const StkFrames& frames) const {
    for (unsigned int i = 0; i < frames.size(); i++) {
        std::cout << "i = " << i << " : output = " << frames[i] << std::endl;
    }
}
