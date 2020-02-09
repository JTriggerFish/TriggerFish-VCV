#include <memory>
#include "TfElements.hpp"
#include "components.hpp"
#include "dsp/noise.hpp"


// Analog style modulation of pitch for VCOs and filter cutoffs
struct TfSlop : Module 
{
	enum ParamIds {
        HUM_LEVEL,
        DRIFT_LEVEL,
        TRACK_SCALING,
		DETUNE_MODE,
		NUM_PARAMS
	};
	enum InputIds {
		VOCT_INPUT,
		NUM_INPUTS
	};
	enum OutputIds {
		VOCT_OUTPUT,
		NUM_OUTPUTS
	};
	enum LightIds {
		NUM_LIGHTS
	};

    std::random_device _seed{};
    std::minstd_rand _rng;
    std::normal_distribution<double> _gaussian{};

    static constexpr float _maxHum{1.0e-2f};
    static constexpr float _humFreq{60.0f};
    
    float _humPhaseIncrement{};
    float _humPhase{};

    //Temperature drift is modelled as an OU process with a simple Euler discretization.
	// i.e an AR(1) process
	static constexpr double _tau{60.0}; //Time constant ( average decay time) in seconds
	static constexpr double _sigmaCents{0.2 / 12}; 
	static constexpr double _sigmaHz{2};
    double _ou{};
    double _phi;
	float _prevDetuneMode{};

    //----------------------------------------------------------------

	TfSlop() : Module(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS),  _rng(_seed())
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


void TfSlop::init(float sampleRate)
{
	double T = 1.0 / sampleRate;
	_humPhaseIncrement = _humFreq * T;
	_phi = 1.0 - T / _tau;
	_gaussian = std::normal_distribution<double>(0.0, std::sqrt(T));
}

void TfSlop::step()
 {
	 if(_prevDetuneMode != params[DETUNE_MODE].value)
	 {
		 _ou = 0.0;
		 _prevDetuneMode = params[DETUNE_MODE].value;
	 }
	float voct = inputs[VOCT_INPUT].value * params[TRACK_SCALING].value;

    _humPhase += _humPhaseIncrement;
    if(_humPhase >= 1.0)
        _humPhase -= 1.0;

    float hum = _maxHum * params[HUM_LEVEL].value * std::sin(2 * PI * _humPhase);

	double sigma =  params[DETUNE_MODE].value < 0 ? _sigmaHz : _sigmaCents;

	_ou = _phi * _ou + sigma * _gaussian(_rng);
	float drift = params[DRIFT_LEVEL].value * _ou;

	voct = voct + hum;

	
	if(params[DETUNE_MODE].value < 0 ) //Hz i.e linear detune mode
		outputs[VOCT_OUTPUT].value = dsp::detune::linear(voct, drift);
	else //Cents i.e proportional detune mode
		outputs[VOCT_OUTPUT].value = voct + drift;


}
void TfSlop::onSampleRateChange()
{
	float gSampleRate = engineGetSampleRate();
	init(gSampleRate);
}


struct TfSlopWidget : ModuleWidget {
	TfSlopWidget(TfSlop *module) : ModuleWidget(module) {
		setPanel(SVG::load(assetPlugin(plugin, "res/TfSlop.svg")));

		//Panel screws
		addChild(Widget::create<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(Widget::create<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(Widget::create<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(Widget::create<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		//Knobs
        addParam(ParamWidget::create<TfCvKnob>(Vec(30, 55), module, TfSlop::HUM_LEVEL, 0.0f, 1.0f, 0.25f));
        addParam(ParamWidget::create<TfCvKnob>(Vec(10, 127), module, TfSlop::DRIFT_LEVEL, 0.0f, 1.0f, 0.25f));
		addParam(ParamWidget::create<TfCvKnob>(Vec(30, 190), module, TfSlop::TRACK_SCALING, 1.0f - 0.2f / 12, 1.0f, 1.0f));

		//Drift mode switch
		addParam(ParamWidget::create<CKSS>(Vec(65, 135), module, TfSlop::DETUNE_MODE, -1.0f, 1.0f, -1.0f));

		//Jacks at the bottom
		addInput(Port::create<PJ301MPort>(Vec(13.5, 317), Port::INPUT, module, TfSlop::VOCT_INPUT));
		addOutput(Port::create<PJ301MPort>(Vec(55, 317), Port::OUTPUT, module, TfSlop::VOCT_OUTPUT));

	}
};


// Specify the Module and ModuleWidget subclass, human-readable
// author name for categorization per plugin, module slug (should never
// change), human-readable module name, and any number of tags
// (found in `include/tags.hpp`) separated by commas.
Model *modelTfSlop = Model::create<TfSlop, TfSlopWidget>("TfSlop");
