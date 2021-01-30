//
//  SoundProcessing.cpp
//  Visualiser
//
//  Created by Rufus Vijayaratnam on 22/10/2020.
//  Copyright Â© 2020 Rufus Vijayaratnam. All rights reserved.
//

#include "SoundProcessing.hpp"

std::vector<std::complex<double>> Audio::FFT(std::vector<std::complex<double>> &samples) {
    int N = samples.size();
	
	if (N == 1) {return samples;};

	int M = N / 2;

	std::vector<std::complex<double>> Xeven(M, 0);
	std::vector<std::complex<double>> Xodd(M, 0);

	for(int i = 0; i < M; i++) {
		Xeven[i] = samples[2 * i];
		Xodd[i] = samples[2 * i + 1];
	}

	std::vector<std::complex<double>> Feven(M, 0);
	Feven = FFT(Xeven);
	std::vector<std::complex<double>> Fodd(M, 0);
	Fodd = FFT(Xodd);
	double pi = 3.14159265359;
	std::vector<std::complex<double>> freqbins(N, 0);
	for(int k = 0; k < (N / 2); k++) {
		std::complex<double> cmplxexponential = std::polar(1.0, -2*pi*k/N) * Fodd[k];
		freqbins[k] = Feven[k] + cmplxexponential;
		freqbins[k+N/2] = Feven[k] - cmplxexponential;
	}
	return freqbins;
}

std::vector<std::vector<std::complex<double>>> Audio::FramedFFT(std::vector<std::vector<std::complex<double>>> &frames) {
	int nFrames = frames.size();
	std::vector<std::vector<std::complex<double>>> frequencyBins;
	for (int i = 0; i < nFrames; i++) {
		std::vector<std::complex<double>> frequencyBin;
		frequencyBin = Audio::FFT(frames[i]);
		frequencyBins.push_back(frequencyBin);
	}
	return frequencyBins;
}


void Audio::TrimAudioSample(std::vector<std::complex<double>> &samples) {
	int extraSamples = samples.size() % BLOCK_SIZE;
	samples = std::vector<std::complex<double>>(samples.begin(), samples.end() - extraSamples);
}

/*
Create a vector containing each frame to be processed subsequently.
If a lower time resolution is desired, this should be achieved by skipping frames
during the graphics phase, not by changing the BLOCK_SIZE.
Each frame is a vector of size BLOCK_SIZE
*/
std::vector<std::vector<std::complex<double>>> Audio::FrameSamples(std::vector<std::complex<double>> &samples) {
	std::vector<std::vector<std::complex<double>>> framedSamples;
	int nFrames = samples.size() / BLOCK_SIZE;
	for (int i = 0; i < nFrames; i++) {
		std::vector<std::complex<double>> frame;
		for (int j = 0; j < BLOCK_SIZE; j++) {
			std::complex<double> element;
			element = samples[i * BLOCK_SIZE + j];
			frame.push_back(element);
		}
		framedSamples.push_back(frame);
	}
	return framedSamples;
}

std::vector<double> Audio::ExtractFrequencies(int sampleRate) {
	std::vector<double> frequencies;
	int nyquistLimitIndex = BLOCK_SIZE / 2;
	for (int i = 0; i < nyquistLimitIndex; i++) {
		double frequency;
		frequency = sampleRate / BLOCK_SIZE * i;
		frequencies.push_back(frequency);
	}
	return frequencies;
}

std::vector<std::vector<std::complex<double>>> Audio::TrimFrames(std::vector<std::vector<std::complex<double>>> &frames) {
	int nyquistLimitIndex = BLOCK_SIZE / 2;
	int nFrames = frames.size();
	for (int i = 0; i < nFrames; i++) {
		frames[i] = std::vector<std::complex<double>>(frames[i].begin(), frames[i].end() - nyquistLimitIndex + 1);
	}
	return frames;
}

std::vector<std::vector<double>> Audio::ExtractMagnitudes(std::vector<std::vector<std::complex<double>>> &frames) {
	std::vector<std::vector<double>> magnitudes;
	int frameSize = BLOCK_SIZE / 2;
	int nFrames = frames.size();
	for (int i = 0; i < nFrames; i++) {
		std::vector<double> frameMagnitudes;
		for (int j = 0; j < frameSize; j++) {
			double magnitude;
			double real = std::real(frames[i][j]);
			double imag = std::imag(frames[i][j]);
			magnitude = std::sqrt(std::pow(real, 2) + std::pow(imag, 2));
			frameMagnitudes.push_back(magnitude);
		}
		magnitudes.push_back(frameMagnitudes);
	}
	return magnitudes;


}

std::vector<std::complex<double>> Audio::ChannelAveraged(AudioFile<double> &song) {
	std::vector<std::complex<double>> samples;
    int numSamples = song.getNumSamplesPerChannel();
    for (int i = 0; i < numSamples; i++) {
        double channelAvg = (song.samples[0][i] + song.samples[1][i]) / 2;
        samples.push_back(std::complex<double>(channelAvg));
    }
	return samples;
}


