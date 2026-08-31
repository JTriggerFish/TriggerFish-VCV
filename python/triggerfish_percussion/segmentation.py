"""Explicit perceptual regions for impact calibration."""

from dataclasses import dataclass


@dataclass(frozen=True)
class PerceptualRegions:
    onset_sample: int
    contact_end: int
    bloom_end: int
    early_body_end: int
    signal_end: int

    def validate(self) -> None:
        boundaries = (
            self.onset_sample,
            self.contact_end,
            self.bloom_end,
            self.early_body_end,
            self.signal_end,
        )
        if boundaries[0] < 0 or any(
            right <= left for left, right in zip(boundaries, boundaries[1:])
        ):
            raise ValueError("perceptual-region boundaries must increase")

    def slices(self) -> dict[str, slice]:
        self.validate()
        return {
            "contact": slice(self.onset_sample, self.contact_end),
            "bloom": slice(self.contact_end, self.bloom_end),
            "early_body": slice(self.bloom_end, self.early_body_end),
            "tail": slice(self.early_body_end, self.signal_end),
        }


def nominal_regions(
    onset_sample: int,
    signal_samples: int,
    sample_rate: float,
    contact_seconds: float = 0.015,
    bloom_seconds: float = 0.120,
    early_body_seconds: float = 0.600,
) -> PerceptualRegions:
    if sample_rate <= 0 or onset_sample < 0 or signal_samples <= onset_sample:
        raise ValueError("invalid signal dimensions for percussion segmentation")
    boundaries = [
        onset_sample,
        onset_sample + round(contact_seconds * sample_rate),
        onset_sample + round(bloom_seconds * sample_rate),
        onset_sample + round(early_body_seconds * sample_rate),
        signal_samples,
    ]
    for index in range(1, len(boundaries)):
        boundaries[index] = min(
            signal_samples, max(boundaries[index], boundaries[index - 1] + 1)
        )
    regions = PerceptualRegions(*boundaries)
    regions.validate()
    return regions
