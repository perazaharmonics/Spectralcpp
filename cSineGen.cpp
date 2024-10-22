#include <cmath>
#include <vector>
#include <algorithm>
#include "FileLoop.h"
#include "FileWvOut.h"

// cDelayline.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

using namespace stk;
using namespace std;

int main()
{
    // Set the global sample rate before creating class instances.
    Stk::setSampleRate(44100.0);

    FileLoop input;
    FileWvOut output;

    // Load the sine wave file.
    try {
		input.openFile("rawwaves/sinewave.raw");
	}
    catch (StkError&) {
		exit(1);
	}

    // Open a 16-bit per-channel, one-channel WAV formatted output file
    try {
        output.openFile("hellosine.wav", 1, FileWrite::FILE_WAV, Stk::STK_SINT16);

    }
    catch (StkError&) {
			exit(1);
		}
    // Run the oscillator for 40000 samples, writing to the output file
    for (int i = 0; i < 40000; i++) {
		output.tick(input.tick());
	}
    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
