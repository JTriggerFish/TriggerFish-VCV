#pragma once
#include "componentlibrary.hpp"
#include <vector>
#include <jansson.h>
#include "rack.hpp"
#include <iostream>

using namespace std;

namespace rack
{

struct TfSlider : SVGSlider
{
	TfSlider()
	{
		Vec margin = Vec(4, 4);
		maxHandlePos = Vec(-1.5, -8).plus(margin);
		minHandlePos = Vec(-1.5, 104).plus(margin);
		background->svg = APP->window->loadSvg(asset::plugin(pluginInstance, "res/slider.svg"));
		background->wrap();
		background->box.pos = margin;
		box.size = background->box.size.plus(margin.mult(2));
		handle->svg = APP->window->loadSvg(asset::plugin(pluginInstance, "res/sliderHandle.svg"));
		handle->wrap();
	}
};
struct TfCvKnob : RoundBlackKnob
{
	TfCvKnob()
	{
		shadow->blurRadius = 2;
	}
};
struct TfLargeAudioKnob : Davies1900hLargeBlackKnob
{
	TfLargeAudioKnob()
	{
		shadow->blurRadius = 4;
	}
};
struct TfAudioKob : Davies1900hBlackKnob
{
	TfAudioKob()
	{
		shadow->blurRadius = 4;
	}
};
struct TfTrimpot : Trimpot
{
	TfTrimpot()
	{
		shadow->blurRadius = 1;
	}
};

struct TfSnapKnob : RoundBlackSnapKnob
{
	TfSnapKnob()
	{
		shadow->blurRadius = 2;
	}
};

struct TfRotarySwitchKnob : RoundBigBlackKnob
{
	TfRotarySwitchKnob()
	{
		snap = true;
		shadow->blurRadius = 3;
	}
};

struct TfSvgWatermark : SvgWidget
{
	Vec drawScale = Vec(1, 1);
	float opacity = 1.0f;

	void setScaledSvg(std::shared_ptr<window::Svg> svg, Vec size)
	{
		setSvg(svg);
		drawScale = Vec(size.x / box.size.x, size.y / box.size.y);
		box.size = size;
	}

	void draw(const DrawArgs& args) override
	{
		nvgGlobalAlpha(args.vg, opacity);
		nvgScale(args.vg, drawScale.x, drawScale.y);
		SvgWidget::draw(args);
	}
};

} // namespace rack
