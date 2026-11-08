"""CLI for generating and optionally playing deterministic multichannel PCM."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from .array_config import load_array
from .signal import pulse, to_s16le, tone, write_wav


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Minimal 8-channel ALSA audio laboratory")
    parser.add_argument("mode", choices=("tone", "pulse"))
    selector = parser.add_mutually_exclusive_group(required=True)
    selector.add_argument("--channel", type=int, help="1-based raw PCM channel (not speaker-test -s)")
    selector.add_argument("--element", help="speaker element ID from --array JSON")
    parser.add_argument("--array", type=Path, help="speaker array JSON configuration")
    parser.add_argument("--channels", type=int, default=8)
    parser.add_argument("--rate", type=int, default=48000)
    parser.add_argument("--seconds", type=float, default=2.0)
    parser.add_argument("--frequency", type=float, default=1000.0, help="tone frequency in Hz")
    parser.add_argument("--gain", type=float, default=0.15)
    parser.add_argument("--output", type=Path, help="write a WAV file instead of playing")
    parser.add_argument("--device", default="hw:CARD=ICUSBAUDIO7D,DEV=0")
    args = parser.parse_args(argv)
    if args.channels < 1:
        parser.error("--channels must be positive")
    if args.element is not None:
        if args.array is None:
            parser.error("--element requires --array")
        try:
            config = load_array(args.array, max_channels=args.channels)
            channel = config.pcm_for(args.element)
        except (OSError, ValueError, KeyError) as exc:
            parser.error(str(exc))
    else:
        channel = args.channel
        if args.array is not None:
            try:
                load_array(args.array, max_channels=args.channels)
            except (OSError, ValueError) as exc:
                parser.error(str(exc))
    if not 1 <= channel <= args.channels:
        parser.error("--channel must be between 1 and --channels")
    try:
        common = dict(channels=args.channels, sample_rate=args.rate, duration=args.seconds, active=(channel - 1,))
        samples = (
            tone(**common, frequency=args.frequency, gain=args.gain)
            if args.mode == "tone"
            else pulse(**common, amplitude=args.gain)
        )
    except ValueError as exc:
        parser.error(str(exc))
    if args.output:
        write_wav(args.output, samples, args.rate)
        print(f"Wrote {args.output}: {len(samples)} frames, {args.channels} channels, {args.rate} Hz")
        return 0

    cmd = [
        "aplay", "-q", "-D", args.device, "-t", "raw", "-f", "S16_LE",
        "-c", str(args.channels), "-r", str(args.rate),
    ]
    label = f"element {args.element} -> " if args.element is not None else ""
    print(f"Sending {label}PCM channel {channel}/{args.channels} directly to {args.device}", file=sys.stderr)
    try:
        subprocess.run(cmd, input=to_s16le(samples), check=True)
    except FileNotFoundError:
        parser.exit(2, "aplay not found (install alsa-utils)\n")
    except subprocess.CalledProcessError as exc:
        parser.exit(2, f"aplay failed with exit code {exc.returncode}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
