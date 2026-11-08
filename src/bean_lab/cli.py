"""CLI for generating and optionally playing deterministic multichannel PCM."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

from .array_config import load_array
from .focus import SPEED_OF_SOUND, focus_delays, route_to_pcm, synthesize_speakers
from .hardware_config import load_hardware
from .signal import pulse, to_s16le, tone, write_wav


FADE_SECONDS = 0.005  # raised-cosine fade at both ends of the focus tone


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
        mono = tone(channels=1, sample_rate=args.rate, duration=args.duration,
                    frequency=args.frequency, gain=args.gain)[:, 0]
        fade = min(round(FADE_SECONDS * args.rate), len(mono) // 2)
        ramp = 0.5 - 0.5 * np.cos(np.pi * (np.arange(fade) + 0.5) / fade)
        mono[:fade] *= ramp
        mono[len(mono) - fade:] *= ramp[::-1]
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


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if argv and argv[0] == "focus":
        return focus_main(argv[1:])
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
