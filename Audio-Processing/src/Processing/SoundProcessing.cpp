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
	
	if (N == 1) {
		return samples;
	} 
	else if (N <= 32) {
		return Audio::DFT(samples);
	}

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
	std::vector<std::complex<double>> freqbins(N, 0);
	for(int k = 0; k < (N / 2); k++) {
		std::complex<double> cmplxexponential = std::polar(1.0, -2*M_PI*k/N) * Fodd[k];
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


std::vector<std::complex<double>> Audio::DFT(std::vector<std::complex<double>> &samples) {
	int N = samples.size();
	std::vector<std::complex<double>> frequency_bin;
	for (int k = 0; k < N; k++) {
		std::vector<std::complex<double>> expnt;
		for (int n = 0; n < N; n++) {
		expnt.push_back(std::polar(1.0, -2 * M_PI * k * n / N));
		}
		std::vector<std::complex<double>> res;
		//What the fuck is res I may ask? Who knows, it's copied from my Python implementation.
		//This is why you should comment your code.
		for (int j = 0; j < N; j++) {
			res.push_back(samples[j] * expnt[j]);
		}
		std::complex<double> sum;
		for (int j = 0; j < N; j++) {
			sum += res[j];
		}
		frequency_bin.push_back(sum);
	}	
	return frequency_bin;
}

//OpenGL coordinates work between -1 and 1 so we must put values into this range with some desired baseline.
//For example, if we want 0 amplitude to appear just off the bottom of the window we can assign that as -0.8
//So that leaves 1.8 remainder for dynamic spectrum. So maximum magnitude must be equal to 1.8.
//But lets actually make it 1.6 so that there is a bit of a header and it doesn't look like the 
//spectrum bars are going off the screen.
std::vector<std::vector<double>> Audio::NormaliseAmplitude(std::vector<std::vector<double>> &frames) {
	int nFrames = frames.size();
	std::vector<std::vector<double>> normalisedFrames;
	for (int i = 0; i < nFrames; i++) {
		double maxAmplitude = *std::max_element(frames[i].begin(), frames[i].end());
		std::vector<double> normalisedFrame;
		for (int j = 0; j < (BLOCK_SIZE / 2); j++) {
			//Line below ensures that all amplitudes will range from 0 to 1.6
			double element = frames[i][j] / (1.0 / 1.6 * maxAmplitude);
			element -= 0.8; //Now all elements are between - 0.8 and 0.8
			normalisedFrame.push_back(element);
		}
			normalisedFrames.push_back(normalisedFrame);
	}
	return normalisedFrames;
}

//int N is the number of frequency representations desired.
std::vector<double> Audio::SpectrumFrequencies(std::vector<double> &frequencies, int N) {
	std::vector<double> spectrumFrequencies;
	const double minFrequency = 50;
	double maxFrequency;
	if (frequencies.back() > 20000) {
		maxFrequency = 20000;
	} else {
		maxFrequency = frequencies.back();
	}

	double span = (maxFrequency - minFrequency) / N;
	for (int i = 0; i <= N; i++) {
		double spectrumFrequency = 50 + span * i;
		spectrumFrequencies.push_back(spectrumFrequency);
	}
	return spectrumFrequencies;
}

//This function takes in the full magnitude vector and the vector of frequency ranges to be represented by each 
//segment of the spectrum.
//It divides the magnitudes into the correct respective ranges. The length of each frame will be N where N
//is the desired number of bars on the spectrum.
std::vector<std::vector<double>> Audio::MagnitudeToSpectrum(std::vector<std::vector<double>> &magnitudes, std::vector<double> &spectrumFrequencies, std::vector<double> &frequencies) {
	std::vector<std::vector<double>> dividedMagnitudes;
	int N = spectrumFrequencies.size() - 1;
	printf("here n is %i \n", N);
	int nFrames = magnitudes.size();
	for (int k = 0; k < nFrames; k++) {

	int prevMinIndex = 0;
	int currentIndex; //These are used to increase efficiency of the loop. To remember where 
		//For single frame
		std::vector<double> dividedMagnitudesFrame;
		for (int i = 0; i < (N - 1); i++) { //Iterate through each frequency range
			std::vector<double> frequencyBandMags;
			double max = spectrumFrequencies[i + 1];
			for (int j = prevMinIndex; j < (BLOCK_SIZE / 2); j++) {
				double mag = magnitudes[k][j];
				double magFrequency = frequencies[j];
				if (magFrequency < max) {
					frequencyBandMags.push_back(mag);
					currentIndex = j;
				} else {
					break;
				}
			}
			double bandAvg = std::accumulate(frequencyBandMags.begin(), frequencyBandMags.end(), 0.0) / ((currentIndex + 1) - prevMinIndex);
			dividedMagnitudesFrame.push_back(bandAvg);
			prevMinIndex = currentIndex;
		}

		dividedMagnitudes.push_back(dividedMagnitudesFrame);
		//Single frame end



	}
	return dividedMagnitudes;
}
