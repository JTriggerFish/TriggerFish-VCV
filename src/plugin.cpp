#include "plugin.hpp"
#include "tfdsp/minblep.hpp"
#include "tfdsp/wavefolder.hpp"

Plugin *pluginInstance;

void init(Plugin *p)
{
	pluginInstance = p;

	// Generate the shared oscillator reconstruction table while Rack loads the
	// plugin. MinBlepGenerator also prepares it on construction as a fallback,
	// but doing it here keeps all one-time FFT work away from audio processing.
	tfdsp::MinBlepGenerator<>::PrepareKernel();
	tfdsp::Wavefolder::PrepareTable();

	// Add all Models defined throughout the pluginInstance
	p->addModel(modelTfVCA);
	p->addModel(modelTfSlop);
	p->addModel(modelTfSlop4);
	p->addModel(modelTfVDPO);
	p->addModel(modelTf303VoiceCore);
	p->addModel(modelTf303Oscillator);
	p->addModel(modelTf4072VoiceCore);
	p->addModel(modelTfWavefoldOscillator);

	// Any other pluginInstance initialization may go here.
	// As an alternative, consider lazy-loading assets and lookup tables when your module is created to reduce startup times of Rack.
}
