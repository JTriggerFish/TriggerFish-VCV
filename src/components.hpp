#pragma once
#include "componentlibrary.hpp"
#include <vector>
#include <jansson.h>
#include "rack.hpp"
#include <iostream>

using namespace std;

namespace rack
{

struct TfSlider : app::SvgSlider
{
	TfSlider()
	{
		constexpr float width = 20.0f;
		constexpr float height = 100.0f;
		constexpr float handleHeight = 13.0f;
		setBackgroundSvg(window::Svg::load(asset::plugin(
			pluginInstance, "res/TfSlider.svg")));
		setHandleSvg(window::Svg::load(asset::plugin(
			pluginInstance, "res/TfSliderHandle.svg")));
		setHandlePosCentered(
			math::Vec(width / 2.0f, height - handleHeight / 2.0f),
			math::Vec(width / 2.0f, handleHeight / 2.0f));
	}
};

struct TfEnvelopeSliderLight : RectangleLight<BlueLight>
{
	TfEnvelopeSliderLight()
	{
		box.size = Vec(5.0f, 7.0f);
		// BlueLight supplies Rack's standard #29b2ef LED colour. Keep the
		// unlit window neutral so it does not muddy the illuminated blue.
		bgColor = nvgRGBA(0x18, 0x18, 0x18, 0xff);
		borderColor = nvgRGBA(0, 0, 0, 0x50);
	}
};

struct TfEnvelopeSlider : LightSlider<TfSlider, TfEnvelopeSliderLight>
{
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
