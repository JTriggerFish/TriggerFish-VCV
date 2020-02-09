#pragma once
#include "componentlibrary.hpp"
#include <vector>
#include <jansson.h>
#include "rack0.hpp"
#include <iostream>

using namespace std;

namespace rack {

struct TfSlider : SVGSlider {
     TfSlider()
     {
		Vec margin = Vec(4, 4);
		maxHandlePos = Vec(-1.5, -8).plus(margin);
		minHandlePos = Vec(-1.5, 104).plus(margin);
		background->svg = SVG::load(assetPlugin(pluginInstance,"res/slider.svg"));
		background->wrap();
		background->box.pos = margin;
		box.size = background->box.size.plus(margin.mult(2));
		handle->svg = SVG::load(assetPlugin(pluginInstance,"res/sliderHandle.svg"));
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

}
