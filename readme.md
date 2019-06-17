# DigitalComm package 

## Purpose 

This package aims to provide some usefull tools to manipulate digital
communication blocks in Julia. 
Currently, the package support the following elements 
- Bit manipulation 
  * Generation of random binary sequence 
  * Conversion between binary sequences and octal sequences 
- Modulation // demodulation
  * Quadrature Amplitude Modulation (QAM) with 4-QAM (QPSK), 16-QAM, 64-QAM and 256-QAM. 
  * Hard demapper for the x-QAM formats 
  * Max log Soft demapper for the x-QAM formats
- Single carrier pulses shapes 
  * Raised Cosine pulse shape 
  * Square root raised Cosine pulse shape 
- Multicarrier Waveform generation and decoding 
  * Support of multicarrier Waveforms: OFDM, UF-OFDM, WOLA, BF-OFDM 

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add DigitalComm
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("DigitalComm")
```

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **documentation of the most recently tagged version.**