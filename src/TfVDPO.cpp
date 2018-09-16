#include "models/VdpOscillator.hpp"
#include <memory>
#include "TfElements.hpp"
#include "components.hpp"
#include "dsp/noise.hpp"


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
    //VdpOscillator<dsp::X2Resampler_Order7, 3> _vdp{dsp::CreateX2Resampler_Chebychev7};
    VdpOscillator<dsp::X4Resampler_Order7, 3> _vdpHq{dsp::CreateX4Resampler_Cheby7};

    //----------------------------------------------------------------

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	TfVDPO() : Module(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS)
	{
		auto engineSampleRate = engineGetSampleRate();
		//_resampler = dsp::CreateX2Resampler_Butterworth5();
		init(engineSampleRate);
	}

	void step() override;
	void init(float sampleRate);
	void onSampleRateChange() override;


	// For more advanced Module features, read Rack's engine.hpp header file
	// - toJson, fromJson: serialization of internal data
	// - onReset, onRandomize, onCreate, onDelete: implements special behavior when user clicks these from the context menu
};


void TfVDPO::init(float sampleRate)
{
    //_vdp.SetSampleRate(sampleRate);
    _vdpHq.SetSampleRate(sampleRate);
}

void TfVDPO::step()
 {
	 float x = inputs[AUDIO_INPUT].value * params[INPUT_GAIN].value;
	 float vOct = inputs[VOCT_INPUT].value * params[VOCT_SCALING].value + params[FREQ].value;
	 float mu = params[DAMPING].value + params[DAMPING_ATTENUVERT].value * inputs[DAMPING_INPUT].value;
	 float freq = 261.626f * powf(2.0f, vOct );

	 //TODO: menu item for low quality, leave high quality by default for now.

	 //auto y = params[HQ_MODE].value > 0 ?
	 //_vdpHq.Step(double(x), double(mu), double(2 * PI * freq))
	 //: _vdp.Step(double(x), double(mu), double(2 * PI * freq));
	 auto y = _vdpHq.Step(double(x), double(mu), double(2 * PI * freq));

	 outputs[OUTPUT].value = y * params[LEVEL].value;

}
void TfVDPO::onSampleRateChange()
{
	float gSampleRate = engineGetSampleRate();
	init(gSampleRate);
}


struct TfVDPOWidget : ModuleWidget {
	TfVDPOWidget(TfVDPO *module) : ModuleWidget(module) {
		setPanel(SVG::load(assetPlugin(plugin, "res/TfVDPO.svg")));

		//Panel screws
		addChild(Widget::create<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(Widget::create<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(Widget::create<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(Widget::create<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		//Knobs
        addParam(ParamWidget::create<TfAudioKob>(Vec(14, 58), module, TfVDPO::FREQ, -5.0f, 5.0f, 0.0f));
        addParam(ParamWidget::create<TfAudioKob>(Vec(14, 112), module, TfVDPO::DAMPING, 0.001f, 9.0f, 0.5f));
		addParam(ParamWidget::create<TfCvKnob>(Vec(18, 175), module, TfVDPO::INPUT_GAIN, 0.0, 1.0f, 1.0f));
		addParam(ParamWidget::create<TfCvKnob>(Vec(76, 175), module, TfVDPO::LEVEL, 0.0, 1.0f, 2.0f));

		addParam(ParamWidget::create<TfTrimpot>(Vec(23, 256), module, TfVDPO::VOCT_SCALING, 0.0f, 1.0f, 2.0f));
		addParam(ParamWidget::create<TfTrimpot>(Vec(81, 256), module, TfVDPO::DAMPING_ATTENUVERT, -1.0f, 1.0f, 1.0f));

		//High quality switch
		//addParam(ParamWidget::create<ToggleSwitch>(Vec(50, 280), module, TfVDPO::HQ_MODE, -1.0f, 1.0f, -1.0f));

		//Jacks at the bottom
		addInput(Port::create<PJ301MPort>(Vec(20, 280), Port::INPUT, module, TfVDPO::VOCT_INPUT));
		addInput(Port::create<PJ301MPort>(Vec(78, 280), Port::INPUT, module, TfVDPO::DAMPING_INPUT));
		addInput(Port::create<PJ301MPort>(Vec(20, 324), Port::INPUT, module, TfVDPO::AUDIO_INPUT));
		addOutput(Port::create<PJ301MPort>(Vec(78, 324), Port::OUTPUT, module, TfVDPO::OUTPUT));

	}
};


// Specify the Module and ModuleWidget subclass, human-readable
// author name for categorization per plugin, module slug (should never
// change), human-readable module name, and any number of tags
// (found in `include/tags.hpp`) separated by commas.
Model *modelTfVDPO = Model::create<TfVDPO, TfVDPOWidget>("TriggerFish-Elements", "TfVDPO", "TriggerFish-VDPO", rack::ModelTag::OSCILLATOR_TAG);
