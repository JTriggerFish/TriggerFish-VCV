#include "models/VdpOscillator.hpp"
#include <memory>
#include "plugin.hpp"
#include "components.hpp"
#include "tfdsp/noise.hpp"


// Analog style modulation of pitch for VCOs and filter cutoffs
struct TfVDPO : Module 
{
	enum ParamIds {
        DAMPING,
        FREQ,
        INPUT_GAIN,
		LEVEL,
		VOCT_SCALING,
		DAMPING_ATTENUVERT,
		NUM_PARAMS
	};
	enum InputIds {
		VOCT_INPUT,
		AUDIO_INPUT,
		DAMPING_INPUT,
		NUM_INPUTS
	};
	enum OutputIds {
        OUTPUT,
		NUM_OUTPUTS
	};
	enum LightIds {
		NUM_LIGHTS
	};

	//NOTE : The order 3 BDF integration is chosen as a decent compromise between stability and accuracy,
	// in particular harmonic tuning ( higher order methods have a tendency to have more inharmonic partials ).
	// in the future it might be worth investing other methods such as Implicit Runge-Kutta for instance Radau II.
	// Alternatively higher oversampling also provides better stability and tuning, which is why a HQ mode is offered
    //VdpOscillator<tfdsp::X2Resampler_Order7, 3> _vdp{tfdsp::CreateX2Resampler_Chebychev7};
    VdpOscillator<tfdsp::X4Resampler_Order7, 3> _vdpHq{tfdsp::CreateX4Resampler_Cheby7};

    //----------------------------------------------------------------

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	TfVDPO()
	{
    config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
    configParam(TfVDPO::FREQ, -5.0f, 5.0f, 0.0f, "");
    configParam(TfVDPO::DAMPING, 0.001f, 9.0f, 0.5f, "");
    configParam(TfVDPO::INPUT_GAIN, 0.0, 1.0f, 1.0f, "");
    configParam(TfVDPO::LEVEL, 0.0, 1.0f, 2.0f, "");
    configParam(TfVDPO::VOCT_SCALING, 0.0f, 1.0f, 2.0f, "");
    configParam(TfVDPO::DAMPING_ATTENUVERT, -1.0f, 1.0f, 1.0f, "");
    //configParam(TfVDPO::HQ_MODE, -1.0f, 1.0f, -1.0f, "");
		//_resampler = tfdsp::CreateX2Resampler_Butterworth5();
    float gSampleRate = APP->engine->getSampleRate();
    init(gSampleRate);
	}

	void process(const ProcessArgs& args) override;
	void init(float sampleRate);
	void onSampleRateChange() override;


	// For more advanced Module features, read Rack's engine.hpp header file
	// - dataToJson, dataFromJson: serialization of internal data
	// - onReset, onRandomize, onCreate, onDelete: implements special behavior when user clicks these from the context menu
};


void TfVDPO::init(float sampleRate)
{
    //_vdp.SetSampleRate(sampleRate);
    _vdpHq.SetSampleRate(sampleRate);
}

void TfVDPO::process(const ProcessArgs& args)
 {
	 float x = inputs[AUDIO_INPUT].getVoltage() * params[INPUT_GAIN].getValue();
	 float vOct = inputs[VOCT_INPUT].getVoltage() * params[VOCT_SCALING].getValue() + params[FREQ].getValue();
	 float mu = params[DAMPING].getValue() + params[DAMPING_ATTENUVERT].getValue() * inputs[DAMPING_INPUT].getVoltage();
	 float freq = 261.626f * powf(2.0f, vOct );

	 //TODO: menu item for low quality, leave high quality by default for now.

	 //auto y = params[HQ_MODE].getValue() > 0 ?
	 //_vdpHq.Step(double(x), double(mu), double(2 * PI * freq))
	 //: _vdp.Step(double(x), double(mu), double(2 * PI * freq));
	 auto y = _vdpHq.Step(double(x), double(mu), double(2 * PI * freq));

	 outputs[OUTPUT].setVoltage(y * params[LEVEL].getValue());

}
void TfVDPO::onSampleRateChange()
{
  float gSampleRate = APP->engine->getSampleRate();
  init(gSampleRate);
}


struct TfVDPOWidget : ModuleWidget {
	TfVDPOWidget(TfVDPO *module) {
		setModule(module);
		setPanel(APP->window->loadSvg(asset::plugin(pluginInstance, "res/TfVDPO.svg")));

		//Panel screws
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		//Knobs
        addParam(createParam<TfAudioKob>(Vec(14, 58), module, TfVDPO::FREQ));
        addParam(createParam<TfAudioKob>(Vec(14, 112), module, TfVDPO::DAMPING));
		addParam(createParam<TfCvKnob>(Vec(18, 175), module, TfVDPO::INPUT_GAIN));
		addParam(createParam<TfCvKnob>(Vec(76, 175), module, TfVDPO::LEVEL));

		addParam(createParam<TfTrimpot>(Vec(23, 256), module, TfVDPO::VOCT_SCALING));
		addParam(createParam<TfTrimpot>(Vec(81, 256), module, TfVDPO::DAMPING_ATTENUVERT));

		//High quality switch
		//addParam(createParam<ToggleSwitch>(Vec(50, 280), module, TfVDPO::HQ_MODE));

		//Jacks at the bottom
		addInput(createInput<PJ301MPort>(Vec(20, 280), module, TfVDPO::VOCT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(78, 280), module, TfVDPO::DAMPING_INPUT));
		addInput(createInput<PJ301MPort>(Vec(20, 324), module, TfVDPO::AUDIO_INPUT));
		addOutput(createOutput<PJ301MPort>(Vec(78, 324), module, TfVDPO::OUTPUT));

	}
};


// Specify the Module and ModuleWidget subclass, human-readable
// author name for categorization per pluginInstance, module slug (should never
// change), human-readable module name, and any number of tags
// (found in `include/tags.hpp`) separated by commas.
Model *modelTfVDPO = createModel<TfVDPO, TfVDPOWidget>("TfVDPO");
