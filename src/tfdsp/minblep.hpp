#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

namespace tfdsp
{

namespace minblep_detail
{

constexpr double Pi = 3.141592653589793238462643383279502884;

template<std::size_t Size>
std::array<std::complex<double>, Size> FourierTransform(
	const std::array<std::complex<double>, Size>& input, bool inverse)
{
	static_assert(Size > 1 && (Size & (Size - 1)) == 0,
		"The minBLEP FFT size must be a power of two");
	auto output = input;
	for (std::size_t index = 1, reversed = 0; index < Size; ++index)
	{
		std::size_t bit = Size >> 1;
		for (; reversed & bit; bit >>= 1)
			reversed ^= bit;
		reversed ^= bit;
		if (index < reversed)
			std::swap(output[index], output[reversed]);
	}

	for (std::size_t length = 2; length <= Size; length <<= 1)
	{
		const double angle = (inverse ? 2.0 : -2.0) * Pi /
			static_cast<double>(length);
		const std::complex<double> rotation(std::cos(angle), std::sin(angle));
		for (std::size_t start = 0; start < Size; start += length)
		{
			std::complex<double> twiddle{1.0, 0.0};
			for (std::size_t offset = 0; offset < length / 2; ++offset)
			{
				const auto even = output[start + offset];
				const auto odd = output[start + offset + length / 2] * twiddle;
				output[start + offset] = even + odd;
				output[start + offset + length / 2] = even - odd;
				twiddle *= rotation;
			}
		}
	}
	if (inverse)
		for (auto& value : output)
			value /= static_cast<double>(Size);
	return output;
}

template<int ZeroCrossings, int TableOversampling>
std::array<double, 2 * ZeroCrossings * TableOversampling + 1>
BuildMinimumPhaseStep()
{
	static_assert(ZeroCrossings >= 2, "A minBLEP needs at least two zero crossings");
	static_assert(TableOversampling >= 2,
		"A minBLEP table needs interpolation oversampling");
	constexpr std::size_t Size = 2 * ZeroCrossings * TableOversampling;
	std::array<std::complex<double>, Size> windowedSinc{};
	for (std::size_t index = 0; index < Size; ++index)
	{
		const double position = -static_cast<double>(ZeroCrossings) +
			2.0 * ZeroCrossings * static_cast<double>(index) /
				static_cast<double>(Size - 1);
		const double sinc = std::abs(position) < 1.0e-14 ? 1.0 :
			std::sin(Pi * position) / (Pi * position);
		const double phase = 2.0 * Pi * static_cast<double>(index) /
			static_cast<double>(Size - 1);
		const double window = 0.35875 - 0.48829 * std::cos(phase) +
			0.14128 * std::cos(2.0 * phase) - 0.01168 * std::cos(3.0 * phase);
		windowedSinc[index] = sinc * window;
	}

	const auto spectrum = FourierTransform(windowedSinc, false);
	std::array<std::complex<double>, Size> logMagnitude{};
	for (std::size_t index = 0; index < Size; ++index)
	{
		const double magnitude = std::max(std::abs(spectrum[index]),
			std::exp(-30.0));
		logMagnitude[index] = std::log(magnitude);
	}
	const auto cepstrum = FourierTransform(logMagnitude, true);
	std::array<std::complex<double>, Size> minimumPhaseCepstrum{};
	minimumPhaseCepstrum[0] = cepstrum[0];
	for (std::size_t index = 1; index < Size / 2; ++index)
		minimumPhaseCepstrum[index] = 2.0 * cepstrum[index];

	const auto minimumPhaseLogSpectrum = FourierTransform(
		minimumPhaseCepstrum, false);
	std::array<std::complex<double>, Size> minimumPhaseSpectrum{};
	for (std::size_t index = 0; index < Size; ++index)
		minimumPhaseSpectrum[index] = std::exp(minimumPhaseLogSpectrum[index]);
	const auto impulse = FourierTransform(minimumPhaseSpectrum, true);

	std::array<double, Size + 1> step{};
	double total = 0.0;
	for (std::size_t index = 0; index < Size; ++index)
	{
		total += impulse[index].real();
		step[index] = total;
	}
	const double normalization = std::abs(total) > 1.0e-15 ? 1.0 / total : 1.0;
	for (std::size_t index = 0; index < Size; ++index)
		step[index] *= normalization;
	step[Size] = 1.0;
	return step;
}

} // namespace minblep_detail

/** Minimum-phase band-limited step correction with fractional event timing.
 *
 * InsertDiscontinuity() accepts the raw signal's signed step magnitude and its
 * position in samples relative to the current frame. Zero means the event is
 * at the current frame; negative values describe a recently observed event.
 * Supporting positions older than one sample is useful when a host-rate edge
 * is corrected inside an oversampled processor.
 */
template<int ZeroCrossings = 8, int TableOversampling = 32,
	typename Sample = double>
class MinBlepGenerator
{
public:
	static constexpr int CorrectionSamples = 2 * ZeroCrossings;
	static constexpr int KernelSamples = 2 * ZeroCrossings * TableOversampling;

	void Reset()
	{
		_buffer.fill(Sample{});
		_position = 0;
	}

	void InsertDiscontinuity(double samplePosition, Sample magnitude)
	{
		if (!std::isfinite(samplePosition) || !std::isfinite(
			static_cast<double>(magnitude)) || samplePosition > 0.0 ||
			samplePosition <= -CorrectionSamples)
			return;

		const auto& kernel = Kernel();
		for (int offset = 0; offset < CorrectionSamples; ++offset)
		{
			const double tablePosition =
				(static_cast<double>(offset) - samplePosition) * TableOversampling;
			if (tablePosition >= KernelSamples)
				break;
			const int tableIndex = std::max(0,
				static_cast<int>(std::floor(tablePosition)));
			const double fraction = tablePosition - tableIndex;
			const double step = kernel[tableIndex] + fraction *
				(kernel[tableIndex + 1] - kernel[tableIndex]);
			const int bufferIndex = (_position + offset) % CorrectionSamples;
			_buffer[bufferIndex] += magnitude * static_cast<Sample>(step - 1.0);
		}
	}

	Sample Process()
	{
		const Sample value = _buffer[_position];
		_buffer[_position] = Sample{};
		_position = (_position + 1) % CorrectionSamples;
		return value;
	}

	static const std::array<double, KernelSamples + 1>& Kernel()
	{
		static const auto kernel = minblep_detail::BuildMinimumPhaseStep<
			ZeroCrossings, TableOversampling>();
		return kernel;
	}

private:
	std::array<Sample, CorrectionSamples> _buffer{};
	int _position{};
};

} // namespace tfdsp
