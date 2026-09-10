// Ear-first guidance; implementation equations belong in the architecture docs.
const ControlHelp = {
  field_phase_tilt: "Move the blur towards the bass or treble without moving pitches, spreading the packets or changing their energy. Negative values soften low rings while keeping the highs clearer; positive values do the opposite. The balance pivots at 1 kHz. Zero keeps the original blur profile; it is not equal blur at every frequency. Needs Phase blur above zero.",
  model_level_db: "Overall synth volume, without changing its character or bloom. Reference selection never adjusts it automatically. Double-click to return to the default level.",
  direct_gain: "How much of the initial stick, mallet or brush contact you hear alongside the ringing body. Raise for a closer, more immediate attack; lower for a body-led sound.",
  impact_tone_noise: "Balances pitched stick ping against broadband contact noise.",
  impact_width: "How long the contact lasts. Shorter gives a sharper tap; longer softens and spreads the attack. This changes the contact, not the ringing tail.",
  impact_chirp_pitch: "Raises or lowers the pitch of the initial ping without retuning the ringing body. Most noticeable when the contact mix favours ping over noise.",
  impact_noise_tilt: "Makes the contact noise brighter or darker. Raise for more hiss and crispness; lower for a rounder attack. This does not change the body's tonal balance.",
  impact_micro_density: "How densely the tiny contacts fill a brush gesture. Raise for a smoother texture; lower for more separated strokes. The effect depends on the chosen implement.",
  velocity_brightness: "How much brighter the body becomes as you strike harder. Raise for a bigger contrast between light taps and hard crashes; lower for a more consistent colour across playing strengths.",
  bloom_rate: "How quickly the ringing spreads through the spectrum. Raise it for a faster crash; lower it for a slower developing bloom. With a dark initial strike, energy usually travels into the highs, but it can spread downward too. Zero turns this movement off.",
  bloom_energy_acceleration: "Changes how energy spreads across the spectrum. Near zero, weakly filled regions let it through readily. Higher settings favour concentrated regions and can hold back the developing highs. This follows the relative energy distribution, not overall loudness; use Energy sensitivity for the hard-versus-soft strike response.",
  bloom_energy_sensitivity: "Makes bloom faster while more energy is stored: harder strikes and repeated hits can spread faster, then settle as the tail fades. Zero removes this energy-dependent speed change, not other velocity effects such as brightness. One makes diffusion speed proportional to total energy; two squares that response. Does not add energy or change the T60 damping law.",
  body_excitation: "How strongly each contact drives the ringing body. Raise for more energy and a stronger response from energy-dependent bloom. To change loudness without changing that behaviour, use Model level instead.",
  field_gain: "How loud the ringing body is compared with the direct contact. This is a listening balance: it does not change the energy driving the bloom.",
  body_brightness: "Where the strike starts the body ringing: negative values favour a dark, low start; positive values excite more highs immediately. For a low note that opens into a bright bloom, start darker and let diffusion carry energy upward. This shapes the strike, not the final output EQ.",
  body_excitation_centre: "Where the initial strike's brightness slope begins. Lower it to leave more of the highs for the bloom to develop; raise it to excite a wider range straight away. Use with Initial excitation tilt.",
  body_tune: "Moves all body tones up or down together, keeping their frequency ratios. One is the painted tuning; two is one octave higher. The contact ping has its own pitch control.",
  field_turbulence: "How much surrounding sizzle a ring near 1 kHz has. Lower keeps a clear pitch; higher shifts energy into nearby tones. Slope makes the bass or treble noisier. The old centre control has been folded into this level without changing your noisiness curve. The slider gives fine control near the useful low values.",
  field_turbulence_slope: "Choose which end sizzles most. Positive values keep the lows clearer and make the highs noisier; negative values do the opposite. Noisiness at 1 kHz stays fixed. This also changes the surrounding tones available for bloom and decay.",
  field_packet_spread: "Spreads the surrounding tones farther from each painted mode. Lower gives a focused ring; higher gives a broader metallic shimmer. This changes their spacing, not how many tones there are. It has little effect on a mode with zero noisiness.",
  field_satellite_density: "How many surrounding tones fill the packets. Lower can reveal individual whistles and beats; higher gives a fuller, denser texture. Unlike Packet spread, this changes the number of tones, not the intended width. Zero leaves only the painted centre tones. Changing it rebuilds the sound.",
  field_doublet_split: "How quickly paired tones pulse: 1 Hz is one beat per second. Try 0.3–1 Hz for slow breathing. This is the speed at 125 Hz in both paired layouts; Beat rate tilt can make higher rings pulse faster. Beat depth changes the strength without changing this speed. Other nearby tones can also beat. Zero removes the deliberate split.",
  field_beat_depth: "Strength of the paired pulsation. In Paired ring this shapes the main ring; in Beating doublets it shapes the surrounding shimmer. Try 0.1–0.3 for gentle movement; one gives equal partners and the deepest pulses. Zero removes the weaker partner, but other nearby tones may still beat. The packet keeps the same excitation energy. Double-click restores 0.3.",
  field_beat_rate_tilt: "Let high packets pulse faster than low packets instead of making them all breathe at one speed. Beat rate sets the speed at 125 Hz; +0.25 doubles it over four octaves. Zero uses the same gap everywhere. Works for Paired ring and Beating doublets. It varies speed between packets, not randomly over time; Pitch wander adds irregular movement if wanted.",
  field_phase_bandwidth: "Softens steady metallic ringing into a less tonal wash. Zero preserves clear tones and their natural beating; a little can tame synthetic whistles; too much can sound like hiss. Noisier high-frequency packets blur more strongly. Blur tilt can soften low rings while keeping the highs clear. It changes coherence, not the damping curve.",
  field_wander_hz: "Makes beating less clockwork by gently moving individual pitches independently. This is the maximum deviation in Hz, equally in bass and treble, including clear centre tones. Try 0.3–1 Hz for subtle movement, not more Phase blur. Zero keeps frequencies fixed; large values can sound detuned. Stored energy and the damping curve are unchanged.",
  field_wander_rate: "How frequently each pitch chooses a new random destination and glides smoothly towards it. Try 0.2–1 changes per second for unhurried movement. Different tones start at different times: this is not a shared vibrato. Only matters when Pitch wander is above zero.",
};

export function helpFor(key) {
  if (ControlHelp[key]) return ControlHelp[key];
  if (key.startsWith("body_decay_seconds_")) {
    return "How long ringing near this frequency takes to fade by 60 dB when left to decay. Raise for a longer ring, lower for tighter damping. Bloom can keep feeding this region, so the audible tail may last longer. The curve applies across all modes.";
  }
  if (key.startsWith("body_decay_frequency_")) {
    return "Where this decay point sits in the spectrum. Move it lower or higher to choose which part of the sound gets this decay time.";
  }
  if (key.startsWith("resolved_frequency_")) {
    return "The pitch at the centre of this painted ring and its surrounding tones. Drag left for a lower ring, right for a higher one. Spread and noisiness determine how widely the surrounding tones extend.";
  }
  if (key.startsWith("resolved_level_")) {
    return "How strongly you hear this packet in the finished sound. Lower a ring that sticks out, or raise one you want to hear more clearly. This does not feed more energy into the bloom. At the exact minimum, the handle is switched off.";
  }
  if (key.startsWith("resolved_turbulence_")) {
    return "Keep this particular ring clearer or noisier than the rest. Zero gives a clear ring (still a beating pair in Paired ring); one follows the main noisiness controls. Raise for more surrounding sizzle, or lower to preserve a bell-like pitch.";
  }
  if (key.startsWith("resolved_allocation_")) {
    return "Give this packet a larger share of the available surrounding tones. Raise it for a fuller cluster here; the other packets then get fewer tones. This does not make the packet louder or wider. Zero keeps only its main ring.";
  }
  if (key === "output_eq_enabled") {
    return "Shape the complete contact and body mix with one final EQ. Bypass to hear the unfiltered mix. Neither setting changes the energy moving inside the instrument.";
  }
  if (key.includes("low_cut")) return "Removes low frequencies from what you hear. Raise to reduce rumble or weight; lower for a fuller bottom end. This does not drain energy from the bloom.";
  if (key.includes("high_cut")) return "Removes high frequencies from what you hear. Lower to soften the top end; raise for more brightness. Use the T60 curve to shorten high ringing instead of just making it quieter.";
  if (key.includes("colour_frequency")) return "Where the colour boost or cut sits. Sweep it to find the part of the sound you want to bring forward or soften.";
  if (key.includes("colour_gain")) return "Boosts or softens the region around Colour frequency. Zero leaves that region unchanged.";
  if (key.includes("colour_q")) return "How narrowly the colour boost or cut focuses on its centre frequency. Higher values affect a narrower region.";
  return "Double-click to return to the default setting.";
}
