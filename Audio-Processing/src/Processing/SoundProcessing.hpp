//
//  SoundProcessing.hpp
//  Visualiser
//
//  Created by Rufus Vijayaratnam on 22/10/2020.
//  Copyright Â© 2020 Rufus Vijayaratnam. All rights reserved.
//

#ifndef SoundProcessing_hpp
#define SoundProcessing_hpp
#include <vector>
#include <complex>
#include <iostream>
#include <string>
#include <Parse/AudioFile.h>
#include <cmath>
#include <numeric>
#define BLOCK_SIZE 1024

namespace Audio {
std::vector<std::complex<double>> FFT(std::vector<std::complex<double>> &samples);

void TrimAudioSample(std::vector<std::complex<double>> &samples);

std::vector<std::vector<std::complex<double>>> FrameSamples(std::vector<std::complex<double>> &samples);

std::vector<double> ExtractFrequencies(int sampleRate);

std::vector<std::vector<std::complex<double>>> TrimFrames(std::vector<std::vector<std::complex<double>>> &frames);

std::vector<std::vector<double>> ExtractMagnitudes(std::vector<std::vector<std::complex<double>>> &frames);

std::vector<std::vector<std::complex<double>>> FramedFFT(std::vector<std::vector<std::complex<double>>> &frames);

std::vector<std::complex<double>> ChannelAveraged(AudioFile<double> &song);

std::vector<std::complex<double>> DFT(std::vector<std::complex<double>> &samples);

void NormaliseAmplitude(std::vector<std::vector<double>> &magnitudes);

std::vector<double> SpectrumFrequencies(std::vector<double> &frequencies, int N);

std::vector<std::vector<double>> MagnitudeToSpectrum(std::vector<std::vector<double>> &magnitudes, std::vector<double> &spectrumFrequencies, std::vector<double> &frequencies);

void LogNormaliseAmplitude(std::vector<std::vector<double>> &magnitudes);

std::vector<int> GetFrequencyIndexes(const std::vector<double> &frequencies, const int minimumFrequency, const int maxFrequency, const int bars, const int sampleRate);

void SpectrumMagnitudes(std::vector<std::vector<double>> &magnitudes, const std::vector<int> &frequencyIndexes);
}
#endif 
