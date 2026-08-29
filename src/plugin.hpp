#pragma once
#include <rack.hpp>


using namespace rack;

// Forward-declare the Plugin, defined in Template.cpp
extern Plugin *pluginInstance;

// Forward-declare each Model, defined in each module source file
extern Model *modelTfVCA;
extern Model *modelTfSlop;
extern Model *modelTfSlop4;
extern Model *modelTfVDPO;
extern Model *modelTf303VoiceCore;
extern Model *modelTf303Oscillator;
extern Model *modelTf4072VoiceCore;
extern Model *modelTfElectricPiano;
extern Model *modelTfRideCymbal;
extern Model *modelTfHiHat;
extern Model *modelTfWavefoldOscillator;
extern Model *modelTfUnisonOscillator;
extern Model *modelTfScenePack4;
extern Model *modelTfEventMerge2;
extern Model *modelTfReverb;
extern Model *modelTfTransport;
extern Model *modelTfProgSequencer;
