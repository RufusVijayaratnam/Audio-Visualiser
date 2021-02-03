#include <iostream>
#include <complex>
#include <vector>
#include <MainGFX/Graphics.hpp>
#include <Processing/SoundProcessing.hpp>
#include <Parse/AudioFile.h>
#include <stdio.h>

GLFWwindow * window;

int main() {
    bool doGfx = false;
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
    magnitudes = Audio::NormaliseAmplitude(magnitudes);




    if(doGfx) {
        bool windowOpened;
        if (!gfx::InitialiseGLFW()) return -1;
        window = gfx::OpenWindow("Music Visualiser", windowOpened);
        if (!windowOpened) return -1;
        gfx::Main(window);
    }

    return 0;
}