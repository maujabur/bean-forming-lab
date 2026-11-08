"""CLI for generating and optionally playing deterministic multichannel PCM."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import math

import numpy as np

from .analysis import check_full_scale, db_ratio, peak_abs, rms, steady_window
from .array_config import load_array
from .focus import SPEED_OF_SOUND, focus_delays, route_to_pcm, synthesize_speakers
from .field import (DEFAULT_FLOOR_DB, complex_field, difference_db, field_maps, grid_axis,
                    level_db)
from .hardware_config import load_hardware
from .microphone import receive
from .signal import faded_tone, pulse, to_s16le, tone, write_wav


def focus_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="bean-lab focus",
                                     description="Static acoustic focusing of a tone at a 3-D point")
    parser.add_argument("--array", type=Path, required=True, help="geometry JSON")
    parser.add_argument("--hardware", type=Path, default=Path("configs/hardware/cmedia-0d8c-0102.json"))
    parser.add_argument("--target", type=float, nargs=3, metavar=("X", "Y", "Z"), required=True,
                        help="focal point in millimeters (same frame as the array)")
    parser.add_argument("--frequency", type=float, default=1000.0)
    parser.add_argument("--duration", type=float, default=2.0, help="tone seconds, before delay tail")
    parser.add_argument("--gain", type=float, default=0.05, help="per-channel peak amplitude (0, 1]")
    parser.add_argument("--rate", type=int, default=48000)
    parser.add_argument("--speed-of-sound", type=float, default=SPEED_OF_SOUND, help="m/s")
    parser.add_argument("--wav", type=Path, help="write a multichannel WAV instead of playing")
    parser.add_argument("--device", help="override ALSA device (default from hardware JSON)")
    args = parser.parse_args(argv)
    try:
        array = load_array(args.array)
        hw = load_hardware(args.hardware)
        focus = tuple(v / 1000 for v in args.target)
        delays = focus_delays(array, focus, args.speed_of_sound)
        mono = faded_tone(sample_rate=args.rate, duration=args.duration,
                          frequency=args.frequency, gain=args.gain)
        samples = route_to_pcm(synthesize_speakers(mono, delays, args.rate), hw)
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))

    print(f"Focus ({args.target[0]:g}, {args.target[1]:g}, {args.target[2]:g}) mm, "
          f"c = {args.speed_of_sound:g} m/s, {args.frequency:g} Hz, gain {args.gain:g}")
    print("speaker  pcm  distance_mm  delay_us  delay_samples")
    for d in delays:
        print(f"{d.speaker:7d}  {hw.pcm_for(d.speaker):3d}  {d.distance * 1000:11.3f}  "
              f"{d.delay * 1e6:8.2f}  {d.delay * args.rate:13.3f}")
    if args.wav:
        write_wav(args.wav, samples, args.rate)
        print(f"Wrote {args.wav}: {len(samples)} frames, {samples.shape[1]} channels, {args.rate} Hz")
        return 0
    device = args.device or hw.device
    cmd = ["aplay", "-q", "-D", device, "-t", "raw", "-f", "S16_LE",
           "-c", str(samples.shape[1]), "-r", str(args.rate)]
    print(f"Sending focused signal to {device}", file=sys.stderr)
    try:
        subprocess.run(cmd, input=to_s16le(samples), check=True)
    except FileNotFoundError:
        parser.exit(2, "aplay not found (install alsa-utils)\n")
    except subprocess.CalledProcessError as exc:
        parser.exit(2, f"aplay failed with exit code {exc.returncode}\n")
    return 0


MM_TRIPLE = dict(type=float, nargs=3, metavar=("X", "Y", "Z"))


def _add_tone_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--array", type=Path, required=True, help="geometry JSON")
    parser.add_argument("--frequency", type=float, default=1000.0)
    parser.add_argument("--duration", type=float, default=2.0, help="emitted tone seconds")
    parser.add_argument("--gain", type=float, default=0.05,
                        help="per-speaker amplitude at the 1 m reference distance, (0, 1]")
    parser.add_argument("--rate", type=int, default=48000)
    parser.add_argument("--speed-of-sound", type=float, default=SPEED_OF_SOUND, help="m/s")


def _mm_to_m(point) -> tuple[float, float, float]:
    return tuple(v / 1000 for v in point)


def _listen_case(args, array, mono, mic_mm, focused: bool, target_mm=None):
    """Virtual-microphone reception of one experiment (all speakers get the same tone)."""
    emission = None
    if focused:
        delays = focus_delays(array, _mm_to_m(target_mm), args.speed_of_sound)
        emission = {d.speaker: d.delay for d in delays}
    speakers = [e.speaker for e in array.elements]
    return receive({s: mono for s in speakers}, array, _mm_to_m(mic_mm), args.rate,
                   emission_delays=emission, speed_of_sound=args.speed_of_sound)


def _write_mono(path: Path, signal: np.ndarray, rate: int) -> None:
    write_wav(path, signal[:, None], rate)


def listen_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="bean-lab listen",
                                     description="Virtual omnidirectional microphone (free field, 1/r)")
    _add_tone_args(parser)
    parser.add_argument("--mode", choices=("simultaneous", "focused"), default="focused")
    parser.add_argument("--target", help="focal point in mm (required for --mode focused)",
                        **MM_TRIPLE)
    parser.add_argument("--mic", required=True, help="microphone position in mm", **MM_TRIPLE)
    parser.add_argument("--window", type=float, nargs=2, metavar=("T0", "T1"),
                        help="analysis window in seconds (default: own steady-state window)")
    parser.add_argument("--wav", type=Path, help="write the received mono signal as WAV")
    args = parser.parse_args(argv)
    try:
        focused = args.mode == "focused"
        if focused and args.target is None:
            parser.error("--target is required for --mode focused")
        array = load_array(args.array)
        mono = faded_tone(sample_rate=args.rate, duration=args.duration,
                          frequency=args.frequency, gain=args.gain)
        rec = _listen_case(args, array, mono, args.mic, focused, args.target)
        peak = check_full_scale(rec.signal)
        if args.window:
            i0, i1 = (round(t * args.rate) for t in args.window)
            if not 0 <= i0 < i1 <= len(rec.signal):
                raise ValueError("window must lie inside the recording")
        else:
            i0, i1 = steady_window([[p.arrival_time for p in rec.paths]], args.duration,
                                   args.rate, args.frequency)
        level = rms(rec.signal[i0:i1])
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))

    print(f"Mode {args.mode}; mic ({args.mic[0]:g}, {args.mic[1]:g}, {args.mic[2]:g}) mm; "
          f"c = {args.speed_of_sound:g} m/s; {args.frequency:g} Hz; gain {args.gain:g} at 1 m")
    print("speaker  distance_mm  travel_ms  emission_us  arrival_ms  gain(1/r)")
    for p in rec.paths:
        print(f"{p.speaker:7d}  {p.distance * 1000:11.3f}  {p.travel_time * 1e3:9.4f}  "
              f"{p.emission_delay * 1e6:11.2f}  {p.arrival_time * 1e3:10.4f}  {p.gain:9.4f}")
    _print_levels(level, peak, args.gain, i0, i1, args.rate)
    if args.wav:
        _write_mono(args.wav, rec.signal, args.rate)
        print(f"Wrote {args.wav}: {len(rec.signal)} frames, 1 channel, {args.rate} Hz")
    return 0


def _print_levels(level, peak, gain, i0, i1, rate) -> None:
    single = gain / math.sqrt(2)  # RMS of one speaker's sine at 1 m
    print(f"Window: {i0 / rate:.9f} to {i1 / rate:.9f} s (samples {i0}:{i1})")
    print(f"RMS {level:.6f} ({db_ratio(level, 1.0):.2f} dBFS, "
          f"{db_ratio(level, single):+.2f} dB re one speaker at 1 m); "
          f"peak {peak:.6f} ({db_ratio(peak, 1.0):.2f} dBFS)")


def compare_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="bean-lab compare",
        description="Experiments A-D: simultaneous/focused emission, mic on/off target")
    _add_tone_args(parser)
    parser.add_argument("--target", required=True, help="focal point and on-target mic, mm",
                        **MM_TRIPLE)
    parser.add_argument("--off-target", required=True, help="off-target mic, mm", **MM_TRIPLE)
    parser.add_argument("--outdir", type=Path, required=True, help="directory for exp-A..D WAVs")
    args = parser.parse_args(argv)
    cases = {"A": (False, args.target, "simultaneous", "target"),
             "B": (True, args.target, "focused", "target"),
             "C": (False, args.off_target, "simultaneous", "off-target"),
             "D": (True, args.off_target, "focused", "off-target")}
    try:
        array = load_array(args.array)
        mono = faded_tone(sample_rate=args.rate, duration=args.duration,
                          frequency=args.frequency, gain=args.gain)
        recs = {k: _listen_case(args, array, mono, mic, foc, args.target)
                for k, (foc, mic, _, _) in cases.items()}
        peaks = {k: check_full_scale(r.signal) for k, r in recs.items()}  # nothing written if any clips
        i0, i1 = steady_window([[p.arrival_time for p in r.paths] for r in recs.values()],
                               args.duration, args.rate, args.frequency)
        levels = {k: rms(r.signal[i0:i1]) for k, r in recs.items()}
        args.outdir.mkdir(parents=True, exist_ok=True)
        for k, r in recs.items():
            _write_mono(args.outdir / f"exp-{k}.wav", r.signal, args.rate)
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))

    print(f"Common window: {i0 / args.rate:.9f} to {i1 / args.rate:.9f} s (samples {i0}:{i1}); "
          f"{args.frequency:g} Hz, gain {args.gain:g} at 1 m, c = {args.speed_of_sound:g} m/s")
    print("exp  emission      mic         rms       dBFS     peak  file")
    for k, (_, _, emission, where) in cases.items():
        print(f"{k:>3}  {emission:<12} {where:<10} {levels[k]:.6f} {db_ratio(levels[k], 1.0):8.2f} "
              f"{peaks[k]:8.6f}  {args.outdir / f'exp-{k}.wav'}")
    for label, a, b in (("B - A  focusing gain at target", "B", "A"),
                        ("D - C  focusing effect off target", "D", "C"),
                        ("A - C  simultaneous, target vs off", "A", "C"),
                        ("B - D  focused, target vs off", "B", "D")):
        print(f"{label}: {db_ratio(levels[a], levels[b]):+.3f} dB")
    return 0


def field_map_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="bean-lab field-map",
        description="2-D steady-state pressure maps (free field, 1/r): simultaneous, focused, difference")
    parser.add_argument("--array", type=Path, required=True, help="geometry JSON")
    parser.add_argument("--target", "--focus", dest="target", required=True,
                        help="focal point in mm", **MM_TRIPLE)
    parser.add_argument("--frequency", type=float, default=4000.0)
    parser.add_argument("--gain", type=float, default=0.05, help="per-speaker amplitude at 1 m")
    parser.add_argument("--speed-of-sound", type=float, default=SPEED_OF_SOUND, help="m/s")
    parser.add_argument("--x-range", type=float, nargs=2, default=(-500.0, 1000.0), metavar=("MIN", "MAX"))
    parser.add_argument("--y-range", type=float, nargs=2, default=(100.0, 2000.0), metavar=("MIN", "MAX"))
    parser.add_argument("--z", type=float, default=0.0, help="plane height in mm")
    parser.add_argument("--nx", type=int, default=151, help="points along x")
    parser.add_argument("--ny", type=int, default=191, help="points along y")
    parser.add_argument("--range-db", type=float, default=40.0,
                        help="color range below the common maximum (display only)")
    parser.add_argument("--vmax", type=float, help="top of the level color scale, dB (display only)")
    parser.add_argument("--diff-limit", type=float, default=10.0,
                        help="difference map spans +/- this many dB (display only)")
    parser.add_argument("--floor-db", type=float, default=DEFAULT_FLOOR_DB,
                        help="dB values are clamped here so exact zeros stay finite")
    parser.add_argument("--probe", action="append", default=[], nargs=3, type=float,
                        metavar=("X", "Y", "Z"), help="print levels at this point (mm); repeatable")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        array = load_array(args.array)
        x = grid_axis(args.x_range[0] / 1000, args.x_range[1] / 1000, args.nx, "x")
        y = grid_axis(args.y_range[0] / 1000, args.y_range[1] / 1000, args.ny, "y")
        if not args.range_db > 0 or not args.diff_limit > 0:
            raise ValueError("--range-db and --diff-limit must be positive")
        focus = _mm_to_m(args.target)
        maps = field_maps(array, focus, args.frequency, x, y, args.z / 1000,
                          speed_of_sound=args.speed_of_sound, gain=args.gain)
        ref = maps.reference
        lev_s = level_db(maps.simultaneous, ref, args.floor_db)
        lev_f = level_db(maps.focused, ref, args.floor_db)
        diff = difference_db(maps.focused, maps.simultaneous, ref, args.floor_db)
        emission = {d.speaker: d.delay for d in focus_delays(array, focus, args.speed_of_sound)}
        pts = [("focus", args.target)] + [(f"probe {i + 1}", p) for i, p in enumerate(args.probe)]
        probe_xyz = np.array([_mm_to_m(p) for _, p in pts])
        kw = dict(speed_of_sound=args.speed_of_sound, gain=args.gain)
        ps = complex_field(array, probe_xyz, args.frequency, **kw)
        pf = complex_field(array, probe_xyz, args.frequency, emission_delays=emission, **kw)
        vmax = args.vmax if args.vmax is not None else max(float(lev_s.max()), float(lev_f.max()))
        vmin = vmax - args.range_db
        args.outdir.mkdir(parents=True, exist_ok=True)
        from .fieldplot import save_difference_map, save_level_map

        x_mm, y_mm = x * 1000, y * 1000
        sub = (f"{args.frequency:g} Hz, focus ({args.target[0]:g}, {args.target[1]:g}, "
               f"{args.target[2]:g}) mm")
        common = (array, args.target, x_mm, y_mm, args.z)
        save_level_map(args.outdir / "field-simultaneous.png", lev_s, *common,
                       f"Simultaneous emission, {sub}", vmin, vmax, mark_max=True)
        save_level_map(args.outdir / "field-focused.png", lev_f, *common,
                       f"Focused emission, {sub}", vmin, vmax, mark_max=True)
        save_difference_map(args.outdir / "field-difference.png", diff, *common,
                            f"Focused - simultaneous (dB), {sub}", args.diff_limit)
        params = dict(array=array.name, frequency_hz=args.frequency, gain=args.gain,
                      speed_of_sound=args.speed_of_sound, focus_mm=list(args.target),
                      z_mm=args.z, floor_db=args.floor_db, reference_amplitude=ref,
                      phasor="p(t)=|P| sin(2*pi*f*t + arg P), levels in dB re one speaker at 1 m")
        np.savez_compressed(args.outdir / "field-data.npz", x_mm=x_mm, y_mm=y_mm,
                            simultaneous=maps.simultaneous, focused=maps.focused,
                            simultaneous_db=lev_s, focused_db=lev_f, difference_db=diff,
                            parameters=json.dumps(params))
    except (OSError, ValueError, KeyError, ImportError) as exc:
        parser.error(str(exc))

    print(f"Wrote field-simultaneous.png, field-focused.png, field-difference.png, field-data.npz "
          f"to {args.outdir}")
    print(f"Grid {args.nx} x {args.ny} ({x_mm[0]:g}..{x_mm[-1]:g} x {y_mm[0]:g}..{y_mm[-1]:g} mm, "
          f"z = {args.z:g} mm); {args.frequency:g} Hz; gain {args.gain:g} at 1 m; "
          f"c = {args.speed_of_sound:g} m/s")
    print(f"Level color scale: {vmin:.2f} to {vmax:.2f} dB; difference scale +/-{args.diff_limit:g} dB; "
          f"dB floor {args.floor_db:g}")
    for name, lev in (("simultaneous", lev_s), ("focused", lev_f)):
        iy, ix = np.unravel_index(int(np.argmax(lev)), lev.shape)
        print(f"Map maximum, {name}: {lev[iy, ix]:+.2f} dB at ({x_mm[ix]:g}, {y_mm[iy]:g}) mm")
    print("point        x_mm     y_mm     z_mm   simultaneous_dB   focused_dB   diff_dB")
    for (name, p), a, b in zip(pts, ps, pf):
        ls, lf = level_db(a, ref, args.floor_db), level_db(b, ref, args.floor_db)
        print(f"{name:<10} {p[0]:8.1f} {p[1]:8.1f} {p[2]:8.1f} {ls:16.2f} {lf:12.2f} {lf - ls:9.2f}")
    return 0


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if argv and argv[0] == "focus":
        return focus_main(argv[1:])
    if argv and argv[0] == "listen":
        return listen_main(argv[1:])
    if argv and argv[0] == "campaign":
        from .campaign import campaign_main

        return campaign_main(argv[1:])
    if argv and argv[0] == "field-map":
        return field_map_main(argv[1:])
    if argv and argv[0] == "compare":
        return compare_main(argv[1:])
    parser = argparse.ArgumentParser(description="Minimal multichannel ALSA audio laboratory")
    parser.add_argument("mode", choices=("tone", "pulse"))
    selector = parser.add_mutually_exclusive_group(required=True)
    selector.add_argument("--channel", type=int, help="diagnostic raw PCM channel (1-based)")
    selector.add_argument("--speaker", type=int, help="stable physical speaker number")
    parser.add_argument("--hardware", type=Path, default=Path("configs/hardware/cmedia-0d8c-0102.json"))
    parser.add_argument("--array", type=Path, help="optional geometry JSON, to validate speaker membership")
    parser.add_argument("--channels", type=int, default=8, help="raw PCM channel count for --channel mode")
    parser.add_argument("--rate", type=int, default=48000)
    parser.add_argument("--seconds", type=float, default=2.0)
    parser.add_argument("--frequency", type=float, default=1000.0)
    parser.add_argument("--gain", type=float, default=0.15)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--device", help="override ALSA device (default from hardware JSON)")
    args = parser.parse_args(argv)
    try:
        if args.speaker is not None:
            hw = load_hardware(args.hardware)
            if args.array is not None:
                array = load_array(args.array)
                array.position_for(args.speaker)
            channel = hw.pcm_for(args.speaker)
            channels = hw.channels
            device = args.device or hw.device
        else:
            channel = args.channel
            channels = args.channels
            device = args.device or "hw:CARD=ICUSBAUDIO7D,DEV=0"
            if args.array is not None:
                parser.error("--array is only valid with --speaker")
        if channels < 1 or not 1 <= channel <= channels:
            parser.error(f"channel must be between 1 and {channels}")
        common = dict(channels=channels, sample_rate=args.rate, duration=args.seconds, active=(channel - 1,))
        samples = (
            tone(**common, frequency=args.frequency, gain=args.gain)
            if args.mode == "tone"
            else pulse(**common, amplitude=args.gain)
        )
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))

    if args.output:
        write_wav(args.output, samples, args.rate)
        print(f"Wrote {args.output}: {len(samples)} frames, {channels} channels, {args.rate} Hz")
        return 0

    cmd = ["aplay", "-q", "-D", device, "-t", "raw", "-f", "S16_LE",
           "-c", str(channels), "-r", str(args.rate)]
    label = f"speaker {args.speaker} -> " if args.speaker is not None else ""
    print(f"Sending {label}PCM channel {channel}/{channels} to {device}", file=sys.stderr)
    try:
        subprocess.run(cmd, input=to_s16le(samples), check=True)
    except FileNotFoundError:
        parser.exit(2, "aplay not found (install alsa-utils)\n")
    except subprocess.CalledProcessError as exc:
        parser.exit(2, f"aplay failed with exit code {exc.returncode}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
