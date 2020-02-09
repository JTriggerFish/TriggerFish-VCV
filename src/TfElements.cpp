#include "TfElements.hpp"


Plugin *plugin;


void init(Plugin *p) {
	plugin = p;

	// Add all Models defined throughout the plugin
	p->addModel(modelTfVCA);
	p->addModel(modelTfSlop);
	p->addModel(modelTfSlop4);
	p->addModel(modelTfVDPO);

	// Any other plugin initialization may go here.
	// As an alternative, consider lazy-loading assets and lookup tables when your module is created to reduce startup times of Rack.
}
