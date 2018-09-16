//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION


//Workaround for crappy numpy headers ??
#define PyTuple_GET_SIZE PyTuple_Size
#define PyTuple_GET_ITEM PyTuple_GetItem


#include <Windows.h>
#include <cmath>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <memory>
#include "../../../src/models/VCAcore.hpp"
#include "../../../src/models/VdpOscillator.hpp"
#include "../../../src/dsp/sampleRate.hpp"

const double e = 2.7182818284590452353602874713527;

double sinh_impl(double x) {
	return (1 - pow(e, (-2 * x))) / (2 * pow(e, -x));
}

double cosh_impl(double x) {
	return (1 + pow(e, (-2 * x))) / (2 * pow(e, -x));
}

double tanh_impl(double x) {
	return sinh_impl(x) / cosh_impl(x);
}

PyObject* tanh_impl(PyObject *, PyObject* o) {
	double x = PyFloat_AsDouble(o);
	double tanh_x = sinh_impl(x) / cosh_impl(x);
	return PyFloat_FromDouble(tanh_x);
}

double *pyvector_to_Carrayptrs(PyArrayObject *arrayin) 
{
	int n = arrayin->dimensions[0];
	return reinterpret_cast<double *>(arrayin->data); 
}
/* ==== Check that PyArrayObject is a double (Float) type and a vector ==============
	return 1 if an error and raise exception */
bool  not_doublevector(PyArrayObject *vec) {
	if (vec->descr->type_num != NPY_DOUBLE || vec->nd != 1) {
		PyErr_SetString(PyExc_ValueError,
			"In not_doublevector: array must be of type Float and 1 dimensional (n).");
		return true;
	}
	return false;
}

template<typename Model, typename Resampler>
PyObject* CallVca(PyObject* args, std::unique_ptr<VCACore<Model,Resampler>> vca)
{
	PyArrayObject *audioIn, *cvIn;  // The python objects to be extracted from the args
	double sampleRate;

	if (!PyArg_ParseTuple(args, "O!O!d", &PyArray_Type, &audioIn, &PyArray_Type, &cvIn, &sampleRate))
		return nullptr;
	if (audioIn == nullptr) return nullptr;
	if (cvIn == nullptr) return nullptr;

	if (not_doublevector(audioIn)) return nullptr;
	if (not_doublevector(cvIn)) return nullptr;

	auto audio = pyvector_to_Carrayptrs(audioIn);
	auto cv = pyvector_to_Carrayptrs(cvIn);
	vca->SetSampleRate(sampleRate);

	int dims[2]{};
	auto n = dims[0] = audioIn->dimensions[0];
	if (n > cvIn->dimensions[0])
		n = cvIn->dimensions[0];
	auto vecout = (PyArrayObject *)PyArray_FromDims(1, dims, NPY_DOUBLE);
	auto cout = pyvector_to_Carrayptrs(vecout);

	for (int i = 0; i < n; ++i)
	{
		cout[i] = vca->Step(audio[i], cv[i], 1.0);
	}

	return PyArray_Return(vecout);
}
template<typename Resampler, int order>
PyObject* CallVdpO(PyObject* args, VdpOscillator<Resampler, order>& vdp)
{
	PyArrayObject *audioIn, *muIn, *wIn;  // The python objects to be extracted from the args
	double sampleRate;

	if (!PyArg_ParseTuple(args, "O!O!O!d", &PyArray_Type, &audioIn, &PyArray_Type, &muIn, &PyArray_Type, &wIn, &sampleRate))
		return nullptr;
	if (audioIn == nullptr) return nullptr;
	if (muIn == nullptr) return nullptr;
	if (wIn == nullptr) return nullptr;

	if (not_doublevector(audioIn)) return nullptr;
	if (not_doublevector(muIn)) return nullptr;
	if (not_doublevector(wIn)) return nullptr;

	auto audio = pyvector_to_Carrayptrs(audioIn);
	auto mu = pyvector_to_Carrayptrs(muIn);
	auto w = pyvector_to_Carrayptrs(wIn);
	vdp.SetSampleRate(sampleRate);

	int dims[2]{};
	auto n = dims[0] = audioIn->dimensions[0];
	if (n > wIn->dimensions[0])
		n = wIn->dimensions[0];
	if (n > muIn->dimensions[0])
		n = muIn->dimensions[0];
	auto vecout = (PyArrayObject *)PyArray_FromDims(1, dims, NPY_DOUBLE);
	auto cout = pyvector_to_Carrayptrs(vecout);

	for (int i = 0; i < n; ++i)
	{
		auto y = vdp.Step(audio[i], mu[i], w[i]);
		cout[i] = y;
	}

	return PyArray_Return(vecout);
}

/// arguments = audio numpy array, cv numpy array
/// return numpy array
PyObject* vca_OTA_noOversampling(PyObject*, PyObject* args)
{
	auto vca = std::unique_ptr<VCA_OTACore<dsp::DummyResampler>>{ new VCA_OTACore<dsp::DummyResampler>(dsp::CreateDummyResampler) };
	return CallVca(args, std::move(vca));
}
/// arguments = audio numpy array, cv numpy array
/// return numpy array
PyObject* vca_OTA_butterworth5(PyObject*, PyObject* args)
{
	auto vca = std::unique_ptr<VCA_OTACore<dsp::X2Resampler_Order5>>{ new VCA_OTACore<dsp::X2Resampler_Order5>(dsp::CreateX2Resampler_Butterworth5) };
	return CallVca(args, std::move(vca));
}
/// arguments = audio numpy array, cv numpy array
/// return numpy array
PyObject* vca_OTA_cheby7(PyObject*, PyObject* args)
{
	auto vca = std::unique_ptr<VCA_OTACore<dsp::X2Resampler_Order7>>{ new VCA_OTACore<dsp::X2Resampler_Order7>(dsp::CreateX2Resampler_Chebychev7) };
	return CallVca(args, std::move(vca));
}
/// arguments = audio numpy array, cv numpy array
/// return numpy array
PyObject* vca_OTA_cheby9(PyObject*, PyObject* args)
{
	auto vca = std::unique_ptr<VCA_OTACore<dsp::X2Resampler_Order9>>{ new VCA_OTACore<dsp::X2Resampler_Order9>(dsp::CreateX2Resampler_Chebychev9) };
	return CallVca(args, std::move(vca));
}
/// arguments = audio numpy array, cv numpy array
/// return numpy array
PyObject* vca_OTA_x4_cheby7(PyObject*, PyObject* args)
{
	auto vca = std::unique_ptr<VCA_OTACore<dsp::X4Resampler_Order7>>{ new VCA_OTACore<dsp::X4Resampler_Order7>(dsp::CreateX4Resampler_Cheby7) };
	return CallVca(args, std::move(vca));
}
PyObject* vca_Transistor_noOversampling(PyObject*, PyObject* args)
{
	auto vca = std::unique_ptr<VCA_TransistorCore<dsp::DummyResampler>>{ new VCA_TransistorCore<dsp::DummyResampler>(dsp::CreateDummyResampler) };
	return CallVca(args, std::move(vca));
}
/// arguments = audio numpy array, cv numpy array
/// return numpy array
PyObject* vca_Transistor_butterworth5(PyObject*, PyObject* args)
{
	auto vca = std::unique_ptr<VCA_TransistorCore<dsp::X2Resampler_Order5>>{ new VCA_TransistorCore<dsp::X2Resampler_Order5>(dsp::CreateX2Resampler_Butterworth5) };
	return CallVca(args, std::move(vca));
}
/// arguments = audio numpy array, cv numpy array
/// return numpy array
PyObject* vca_Transistor_cheby7(PyObject*, PyObject* args)
{
	auto vca = std::unique_ptr<VCA_TransistorCore<dsp::X2Resampler_Order7>>{ new VCA_TransistorCore<dsp::X2Resampler_Order7>(dsp::CreateX2Resampler_Chebychev7) };
	return CallVca(args, std::move(vca));
}
/// arguments = audio numpy array, cv numpy array
/// return numpy array
PyObject* vca_Transistor_cheby9(PyObject*, PyObject* args)
{
	auto vca = std::unique_ptr<VCA_TransistorCore<dsp::X2Resampler_Order9>>{ new VCA_TransistorCore<dsp::X2Resampler_Order9>(dsp::CreateX2Resampler_Chebychev9) };
	return CallVca(args, std::move(vca));
}
/// arguments = audio numpy array, cv numpy array
/// return numpy array
PyObject* vca_Transistor_x4_cheby7(PyObject*, PyObject* args)
{
	auto vca = std::unique_ptr<VCA_TransistorCore<dsp::X4Resampler_Order7>>{ new VCA_TransistorCore<dsp::X4Resampler_Order7>(dsp::CreateX4Resampler_Cheby7) };
	return CallVca(args, std::move(vca));
}
PyObject* vdpO(PyObject*, PyObject* args)
{
	VdpOscillator<dsp::X4Resampler_Order7, 3> vdp{dsp::CreateX4Resampler_Cheby7};
	return CallVdpO < dsp::X4Resampler_Order7, 3 > (args, vdp);
}


static PyMethodDef vcv_methods[] = {
	// The first property is the name exposed to Python, fast_tanh, the second is the C++
	// function name that contains the implementation.
	//{ "fast_tanh", (PyCFunction)tanh_impl, METH_O, nullptr },
	{ "vca_OTA_noOversampling", static_cast<PyCFunction>(vca_OTA_noOversampling), METH_VARARGS, nullptr },
	{ "vca_OTA_butterworth5", static_cast<PyCFunction>(vca_OTA_butterworth5), METH_VARARGS, nullptr },
	{ "vca_OTA_cheby7", static_cast<PyCFunction>(vca_OTA_cheby7), METH_VARARGS, nullptr },
	{ "vca_OTA_cheby9", static_cast<PyCFunction>(vca_OTA_cheby9), METH_VARARGS, nullptr },
	{ "vca_OTA_x4_cheby7", static_cast<PyCFunction>(vca_OTA_x4_cheby7), METH_VARARGS, nullptr },

	{ "vca_Transistor_noOversampling", static_cast<PyCFunction>(vca_Transistor_noOversampling), METH_VARARGS, nullptr },
	{ "vca_Transistor_butterworth5", static_cast<PyCFunction>(vca_Transistor_butterworth5), METH_VARARGS, nullptr },
	{ "vca_Transistor_cheby7", static_cast<PyCFunction>(vca_Transistor_cheby7), METH_VARARGS, nullptr },
	{ "vca_Transistor_cheby9", static_cast<PyCFunction>(vca_Transistor_cheby9), METH_VARARGS, nullptr },
	{ "vca_Transistor_x4_cheby7", static_cast<PyCFunction>(vca_Transistor_x4_cheby7), METH_VARARGS, nullptr },
	{ "vdpO", static_cast<PyCFunction>(vdpO), METH_VARARGS, nullptr },

	// Terminate the array with an object containing nulls.
	{ nullptr, nullptr, 0, nullptr }
};
/**
 * \brief Python module to test triggerfish vcv rack models
 */
static PyModuleDef triggerfishvcv_module = {
	PyModuleDef_HEAD_INIT,
	"triggerfishvcv",                        // Module name to use with Python import statements
	"Test models from python",				// Module description
	0,
	vcv_methods								// Structure that defines the methods of the module
};
PyMODINIT_FUNC PyInit_triggerfishvcv() {
	import_array();
	return PyModule_Create(&triggerfishvcv_module);
}