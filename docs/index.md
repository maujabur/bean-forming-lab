# Bean Forming Lab

*[maujabur/bean-forming-lab](https://github.com/maujabur/bean-forming-lab) · made with Kiro*

Open experimental laboratory for multichannel audio, loudspeaker arrays, beam steering and near-field focusing. Everything here is a **simulation in a free-field, ideal point-source model**: predictions, not measurements. Source code, tests and the command-line tools are in the [repository](https://github.com/maujabur/bean-forming-lab).

## Setup and conventions

- [Speaker identity, geometry and hardware mapping](array-configuration.md)
- [Speaker geometry JSON](grid-geometry.md)

## Studies, in the order they were built

1. [Static acoustic focusing](static-focusing.md): real time delays, fractional-delay synthesis, multichannel WAV.
2. [Virtual microphone](virtual-microphone.md): what an ideal omnidirectional microphone would receive.
3. [Two-dimensional field maps](field-map.md): steady-state pressure maps of a tone.
4. [Parametric beamforming study](parametric-study.md): frequency, spacing and focal direction.
5. [Broadband acoustic field study](broadband-study.md): energy maps for band-limited signals.
6. [Three-dimensional field of the arrays](three-dimensional-field.md): orthogonal cuts through the focus, line versus planar arrays.
7. [Geometric optimization of loudspeaker arrays](array-optimization.md): regular, irregular and optimized layouts with 8, 16 and 32 elements.
