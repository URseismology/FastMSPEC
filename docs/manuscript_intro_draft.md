# Manuscript introduction: working draft

Not part of any specific stage's deliverable — a running draft of the opening paragraph for the
eventual FastMSPEC manuscript, refined in conversation across the Round 2 hypothesis-evaluation
work (see `docs/round2_hypothesis_evaluation.tex`'s Significance section, which draws on the same
ideas at report scale, and `docs/notebook5_revamp_progress.md`'s 2026-09-04 log for how these
arguments were derived). Kept here so it doesn't only exist in conversation history.

**Style note**: short, clean sentences, no em-dashes — the user's own stated preference, applied
here and worth carrying into the eventual manuscript draft generally, not just this paragraph.

## Technical version

For decades, earth scientists have used earthquake and ambient-noise correlations to image the
lithosphere, monitor groundwater, and probe the deep Earth. Yet how these correlations should
actually be computed has never been chosen on principled grounds. At its core, the problem is
deciding what spectral resolution, or smoothing bandwidth, to use to improve signal recovery. In
this paper we return to the classical theory of time-frequency concentration, developed by
Slepian and Thomson, and extend it to this problem. We derive a distance-aware bandwidth
criterion for a global seismic array spanning heterogeneous station spacing and data quality. The
resulting framework improves resolution without trading away scale, correctness, or
affordability. The cost of this estimator, FastMSPEC, is provably bounded at roughly two dozen
correction tapers across three orders of magnitude of bandwidth (NW = 5 to 1600). This cost is
already decoupled from the resolution it delivers. A first test against a real 380-station-pair
deployment shows the same distance-dependent criterion correctly anticipating which pairs are
structurally unreachable. This is an early, encouraging signal for the fully principled, per-pair
bandwidth selection framework this paper develops.

## Non-technical version

For decades, earth scientists have used the seismic signals recorded at pairs of stations, set
off by earthquakes or picked out of the planet's constant background noise, to map the structure
of the crust, track groundwater, and probe the deep interior. To compare two stations' recordings
properly, the signal has to be smoothed across frequency in a particular way. Too much smoothing
blurs together features that are actually distinct. Too little leaves a result too noisy to
trust. Getting that smoothing setting right has never been done on principled grounds for this
specific problem. Crucially, the right setting is not the same for every pair of stations. It
depends on how far apart they are, which can vary a hundredfold or more across a real network. In
this paper we return to a decades-old mathematical theory, developed by Slepian and Thomson,
about how to best extract this kind of frequency information from a limited stretch of recorded
data, and extend it to account for that distance dependence directly. The payoff is a method that
adjusts its smoothing correctly for every station pair without becoming expensive to run, even
when a very fine or very coarse setting is needed. The extra computation involved stays small
across an enormous, roughly 300-fold range of settings. Tested on a real set of 380 station
pairs, this distance-based rule already correctly identifies which pairs are simply too far apart
for any smoothing setting to work at all. This is an early, encouraging result for the complete
method this paper is building toward.

## Status

Both versions are drafts, not final. Deliberately avoid overclaiming completion: the bias-variance
lower bound (`NW_low`) is still undetermined (`docs/stage5_bandwidth_theory.tex`), and the one
direct empirical test of the window-length recommendation
(`verification/window_length_splice_test/`) came back negative for the single pair tested. The
closing sentences reflect this — "an early, encouraging signal," not "we show."
