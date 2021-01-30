#include <iostream>
#include <complex>
#include <vector>
#include <FFT/SoundProcessing.hpp>
#include <Parse/AudioFile.h>
#include <gsl/gsl_wavelet.h>
#include <stdio.h>

int main() {
    AudioFile<double> Song;
    std::string directory = "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Sounds/";
    std::string audioFile = "test-audio.wav";
    std::string filePath = directory + audioFile;
    Song.load(filePath);
   
    std::vector<std::vector<std::complex<double>>> frames, frequencyBins;
    std::vector<std::complex<double>> samples;
    samples = Audio::ChannelAveraged(Song);
    
    frames = Audio::FrameSamples(samples);
    frequencyBins = Audio::FramedFFT(frames);
    frames = Audio::TrimFrames(frequencyBins);
    std::vector<double> frequencies = Audio::ExtractFrequencies(Song.getSampleRate());
    std::vector<std::vector<double>> magnitudes = Audio::ExtractMagnitudes(frames);

    printf("Shape of bins: (%i, %i) \n", int(frequencyBins.size()), int(frequencyBins[3].size()));

    for (int i = 0; i < 10; i++) {
        std::cout << frequencyBins[0][i] << std::endl;
    }


    return 0;
}