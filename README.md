# DaSpeak Formant Synthesizer (Compact Guide)

## What this is

DaSpeak is a simple formant-based speech synthesizer written in C++.

It reads phoneme input and produces a mono WAV speech file.

It also supports **interactive CLI**.

---

## File Names

- Source file: `daspeak.cpp`
- Compacted version: `daspeak_compacted.cpp`

Use the compacted version if you want faster compilation or cleaner code.

---

## How to Compile (g++)

Make sure you have **g++** installed.

### Compile the compacted version (recommended)

```bash
g++ -std=c++17 -O3 -s daspeak_compacted.cpp -o daspeak
```

### Compile the original version

```bash
g++ -std=c++17 -O3 -s daspeak.cpp -o daspeak
```

---

## How to Run

You must provide **either** phoneme mode OR specification mode.

### 1. Phoneme Mode (Easier)

Just list phonemes separated by spaces.

Example:

```bash
./daspeak -phon="HH EH L OW" -o=hello.wav
```

This will generate:

- Default timing
- Default pitch
- Automatic smoothing

---

### 2. Specification Mode (Advanced)

Allows manual control.

Format:

```
PHONEME DURATION OVERLAP PITCH1 PITCH2 ...
```

Example:

```bash
./daspeak -spec="HH 0.12 0.015 95|EH 0.14 0.018 105 110|L 0.10 0.015 100|OW 0.20 0.02 115 100" -o=output.wav
```

Rules:

- Phoneme Specs seperated by | (pipe/vertical line)
- Duration = seconds (0.01 to 2.0)
- Overlap = seconds (0.0 to 0.5)
- You can add up to 8 pitch values

---

## Basic Examples

Interactive mode:

```bash
./daspeak
```

Interactive mode opens when no arguments are provided.

Generate speech:

```bash
./daspeak -phon="S AH N D" -o=sound.wav
```

Change sample rate:

```bash
./daspeak -phon="IY Z IY" -r=96000 -o=high.wav
```

Change volume:

```bash
./daspeak -phon="HH AE P IY" -v=6 -o=loud.wav
```

---

## Output

- Produces **16-bit mono WAV file**
- Default sample rate is **44100 Hz**
- Output file default name is **output.wav**

---

## Help

Run:

```bash
./daspeak -h
```

---

## Notes

- Only use phonemes listed in the code.
- If input is invalid, silence is inserted.
- This is a research-style synthesizer, not production TTS.

---

## License

MIT License
