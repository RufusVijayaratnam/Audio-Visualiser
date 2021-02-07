#include <iostream>
#include <complex>
#include <vector>
#include <MainGFX/Graphics.hpp>
#include <Processing/SoundProcessing.hpp>
#include <Parse/AudioFile.h>
#include <stdio.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_mixer.h>


GLFWwindow * window;

int main() {
    bool doGfx = true;
    bool doAudio = true;
    int spectrumBars = 20;
    AudioFile<double> Song;
    std::string directory = "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Sounds/";
    std::string audioFile = "test-audio.wav";
    std::string filePath = directory + audioFile;
    Song.load(filePath);
    int sampleRate = Song.getSampleRate();
    double frameTime = 1024.0 / double(sampleRate);
   
    std::vector<std::vector<std::complex<double>>> frames, frequencyBins;
    std::vector<std::complex<double>> samples;
    samples = Audio::ChannelAveraged(Song);
    Audio::TrimAudioSample(samples);
    
    frames = Audio::FrameSamples(samples);
    frequencyBins = Audio::FramedFFT(frames);
    frames = Audio::TrimFrames(frequencyBins);
    std::vector<double> frequencies = Audio::ExtractFrequencies(Song.getSampleRate());
    std::vector<std::vector<double>> magnitudes = Audio::ExtractMagnitudes(frames);
    Audio::NormaliseAmplitude(magnitudes);
    //std::vector<double> spectrumFrequencies = Audio::SpectrumFrequencies(frequencies, spectrumBars);
   // std::vector<std::vector<double>> spectrumMagnitudes = Audio::MagnitudeToSpectrum(magnitudes, spectrumFrequencies, frequencies);
    
    

    
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO) != 0) { 
    fprintf(stderr, "Unable to initialize SDL: %s\n", SDL_GetError()); return 1; 
    }

    Uint16 audio_format = AUDIO_S16SYS; 
    int audio_channels = 2; 
    int audio_buffers = 4096;

    if(Mix_OpenAudio(sampleRate, audio_format, audio_channels, audio_buffers) != 0) { 
        fprintf(stderr, "Unable to initialize audio: %s\n", Mix_GetError()); exit(1); 
    }


    Mix_Chunk* sound = NULL;

    sound = Mix_LoadWAV(filePath.c_str());
    if (sound == NULL) {
        fprintf(stderr, "unable to locate WAV file: %s\n", Mix_GetError());
    }

   
  

    if(doGfx) {
        bool windowOpened;
        if (!gfx::InitialiseGLFW()) return -1;
        window = gfx::OpenWindow("Music Visualiser", windowOpened);
        if (!windowOpened) return -1;
        gfx::Main(window, magnitudes, frameTime, sound);
    }

    return 0;
}