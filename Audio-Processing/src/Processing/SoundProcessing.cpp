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
		frames[i] = std::vector<std::complex<double>>(frames[i].begin(), frames[i].end() - nyquistLimitIndex);
	}
	return frames;
}

std::vector<std::vector<double>> Audio::ExtractMagnitudes(std::vector<std::vector<std::complex<double>>> &frames) {
	std::vector<std::vector<double>> magnitudes;
	int frameSize = frames[0].size();
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
	int numChannels = song.getNumChannels();
	if (numChannels == 1) {
		for (int i = 0; i < numSamples; i++) {
			samples.push_back(song.samples[0][i]);
		}
		return samples;
	}
    for (int i = 0; i < numSamples; i++) {
        double channelAvg = (song.samples[0][i] + song.samples[1][i]) / 2;
        samples.push_back(std::complex<double>(channelAvg) * 32767.0);
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
void Audio::NormaliseAmplitude(std::vector<std::vector<double>> &magnitudes) {
	int nFrames = magnitudes.size();
	int frameSize = magnitudes[0].size();
	std::vector<std::vector<double>> normalisedFrames;
	for (int i = 0; i < nFrames; i++) {
		double maxAmplitude = *std::max_element(magnitudes[i].begin(), magnitudes[i].end());
		std::vector<double> normalisedFrame;
		for (int j = 0; j < frameSize; j++) {
			//Line below ensures that all amplitudes will range from 0 to 1.6
			double element = magnitudes[i][j] / (1.0 / 1.6 * maxAmplitude);
			element -= 0.8; //Now all elements are between - 0.8 and 0.8
			normalisedFrame.push_back(element);
		}
			//normalisedFrames.push_back(normalisedFrame);
			magnitudes[i] = normalisedFrame;
	}
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

	double span = (maxFrequency - minFrequency) / (N);
	for (int i = 0; i <= (N); i++) {
		double spectrumFrequency = 50 + span * i;
		spectrumFrequencies.push_back(spectrumFrequency);
	}
	return spectrumFrequencies;
}


void Audio::LogNormaliseAmplitude(std::vector<std::vector<double>> &magnitudes) {
	int nFrames = magnitudes.size();
	int frameSize = magnitudes[0].size();
	std::vector<std::vector<double>> normalisedFrames;
	for (int i = 0; i < nFrames; i++) {
		double maxAmplitude = *std::max_element(magnitudes[i].begin(), magnitudes[i].end());
		double minAmplitude = *std::min_element(magnitudes[i].begin(), magnitudes[i].end());

		std::vector<double> normalisedFrame;
		for (int j = 0; j < frameSize; j++) {
			//Line below ensures that all amplitudes will range from 0 to 1.6
			double logNormalisedMag = std::log10(magnitudes[i][j] - minAmplitude) / std::log10(maxAmplitude - minAmplitude);
			//logNormalisedMag -= 0.8;
			normalisedFrame.push_back(logNormalisedMag);
		}
			magnitudes[i] = normalisedFrame;
	}
}

std::vector<int> Audio::GetFrequencyIndexes(const std::vector<double> &frequencies, const int minFrequency, const int maxFrequency, const int bars, const int sampleRate) {
	std::vector<int> frequencyIndexes;
	int minIndex = std::floor(double(minFrequency) / double(sampleRate) * BLOCK_SIZE);
	int maxIndex = std::ceil(double(maxFrequency) / double(sampleRate) * BLOCK_SIZE);
	
	double indexStep = double(maxIndex - minIndex) / double(bars);
	frequencyIndexes.push_back(minIndex);
	for (int i = 1; i <= bars; i++) {
		int index = std::round(double(minIndex) + double(i) * indexStep);
		frequencyIndexes.push_back(index);
	}
	if (frequencyIndexes.back() > frequencies.size()) frequencyIndexes.back() = frequencies.size(); 
	return frequencyIndexes;
}

void Audio::SpectrumMagnitudes(std::vector<std::vector<double>> &magnitudes, const std::vector<int> &frequencyIndexes, int mode) {
	std::vector<std::vector<double>> spectrumMagnitudes;
	const std::vector<std::vector<double>> tempMags = magnitudes;

	if (mode==0) {
		for (int j = 0; j < tempMags.size(); j++) {
			for (int i = 0; i < (frequencyIndexes.size() - 1); i++) {
				int begin = frequencyIndexes[i];
				int end = tempMags[j].size() - frequencyIndexes[i + 1];
				double sum = std::accumulate(tempMags[j].begin() + begin, tempMags[j].end() - end + 1, 0.0);
				double avg = sum / (frequencyIndexes[i + 1] + 1 - begin);
				magnitudes[j][i] = avg;
				//printf("j: %i i: %i, begin: %i, end: %i, sum: %f, avg %f \n", j, i, begin, frequencyIndexes[i + 1], sum, avg);
			}
				magnitudes[j].resize(frequencyIndexes.size() - 1);

		}
	} else if (mode==1) {
		for (int j = 0; j < tempMags.size(); j++) {
			for (int i = 0; i < (frequencyIndexes.size() - 1); i++) {
				int begin = frequencyIndexes[i];
				int end = tempMags[j].size() - frequencyIndexes[i + 1];
				double max = *std::max_element(tempMags[j].begin() + begin, tempMags[j].end() - end + 1);
				magnitudes[j][i] = max;
				}
				magnitudes[j].resize(frequencyIndexes.size() - 1);
			}
	}
	
}

void Audio::ThreshholdMagnitudes(std::vector<std::vector<double>> &magnitudes, double val, double ceil, double floor) {
    for (int j = 0; j < magnitudes.size(); j++) {
        for (int i = 0; i < magnitudes[j].size(); i++) {
            double &amp = magnitudes[j][i];
            if (amp >= val) {
                amp = amp;
            } else {
                amp = floor;
            }
        }
    }
}