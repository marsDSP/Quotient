## Summary
<!-- One or two sentences: what does this PR change and why? -->

## Type of change
<!-- Tick all that apply. -->
- [ ] Bug fix (non-breaking)
- [ ] New feature (non-breaking)
- [ ] Breaking change (parameter rename/removal, state format change, etc.)
- [ ] DSP / algorithmic change
- [ ] UI / editor change
- [ ] Build / CI / tooling
- [ ] Documentation only

## DSP / audio risk
<!--
If this touches the audio path, answer briefly:
- Filters, dynamics, or topology affected?
- Does it change latency, group delay, or phase response?
- Stability across sample rates (44.1 → 192 kHz) and buffer sizes (16 → 4096)?
- Real-time safe (no locks/allocations on the audio thread)?
- Denormal-safe? Reset behaviour on `prepareToPlay` / `releaseResources`?
-->

## Testing
<!--
Describe what you tested and how. Include hosts, formats, sample rates,
buffer sizes, and inputs where relevant.
-->
- [ ] Pluginval passes locally at strictness ≥ 8
- [ ] Manual smoke test in at least one DAW (note which):
- [ ] Standalone build runs

## Screenshots / spectra
<!-- Optional: spectra, oscilloscope captures, A/B with reference, etc. -->

## Checklist
- [ ] Code compiles cleanly on macOS, Linux, and Windows (CI is green)
- [ ] No new compiler warnings introduced
- [ ] Public/parameter API documented if changed
- [ ] State save/restore tested if parameters changed
- [ ] Version bumped if this is a release-worthy change

## Related issues
<!-- e.g. Closes #123, Refs #456 -->
