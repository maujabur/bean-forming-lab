"""CLI for generating and optionally playing deterministic multichannel PCM."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from .signal import pulse, to_s16le, tone, write_wav


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Minimal 8-channel ALSA audio laboratory")
    parser.add_argument("mode", choices=("tone", "pulse"))
    parser.add_argument("--channel", type=int, required=True, help="PCM channel index (1-based; not speaker-test -s)")
    parser.add_argument("--channels", type=int, default=8)
    parser.add_argument("--rate", type=int, default=48000)
    parser.add_argument("--seconds", type=float, default=2.0)
    parser.add_argument("--frequency", type=float, default=1000.0, help="tone frequency in Hz")
    parser.add_argument("--gain", type=float, default=0.15)
    parser.add_argument("--output", type=Path, help="write a WAV file instead of playing")
    parser.add_argument("--device", default="hw:CARD=ICUSBAUDIO7D,DEV=0")
    args = parser.parse_args(argv)
    if not 1 <= args.channel <= args.channels:
        parser.error("--channel must be between 1 and --channels")
    try:
        common = dict(channels=args.channels, sample_rate=args.rate, duration=args.seconds, active=(args.channel - 1,))
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
    print(f"Sending PCM channel {args.channel}/{args.channels} directly to {args.device}", file=sys.stderr)
    try:
        subprocess.run(cmd, input=to_s16le(samples), check=True)
    except FileNotFoundError:
        parser.exit(2, "aplay not found (install alsa-utils)\n")
    except subprocess.CalledProcessError as exc:
        parser.exit(2, f"aplay failed with exit code {exc.returncode}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
